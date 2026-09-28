"""Onshore wind constraint: the strongest onshore wind component between the surface and
`layer_top_ft` at LC-39A and SLC-40, from every model that reaches that layer.

    python -m pipeline.onshore [--only hrrr,gefs]

Writes docs/data/onshore.json. Runs after pipeline.run in the same workflow and keeps its
own cache (cache/onshore.json), so the static-fire output is untouched.

Onshore component (m/s) = -(speed x cos(direction - coast_normal)): negative when the wind
blows from the sea toward land, matching the constraint's sign convention. For each
member-hour the most negative value over the layer is kept, using every model level from
10 m to the layer top plus the wind interpolated to the layer top (log-height between the
levels either side of it). If no level reaches the top, the highest level (at least
`min_top_level_m`) is power-law scaled to it; models whose highest level is lower than that
are left out.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import yaml

from . import grib
from .common import SESSION, Context, PointCache, SourceResult, floor_hour, iso, log, sd_from_uv, uv_from_sd
from .grib import Missing, earth_relative, exists, match_fields, read_idx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HRRR_BASE = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
HRRR_SFC = "hrrr.{ymd}/conus/hrrr.t{hh}z.wrfsfcf{fh:02d}.grib2"
HRRR_PRS = "hrrr.{ymd}/conus/hrrr.t{hh}z.wrfprsf{fh:02d}.grib2"
PLEVS = (1000, 975, 950)          # 950 mb (~500 m) is the first level above 1000 ft

FIELDS = {"u10": r":UGRD:10 m above ground:", "v10": r":VGRD:10 m above ground:",
          "u80": r":UGRD:80 m above ground:", "v80": r":VGRD:80 m above ground:"}
for _p in PLEVS:
    FIELDS[f"u{_p}"] = rf":UGRD:{_p} mb:"
    FIELDS[f"v{_p}"] = rf":VGRD:{_p} mb:"
    FIELDS[f"h{_p}"] = rf":HGT:{_p} mb:"


# ---------------------------------------------------------------- physics
def onshore(u: float, v: float, normal_deg: float) -> float:
    """Negative = onshore. u, v are the wind's east/north components (m/s)."""
    n = math.radians(normal_deg)
    # component of the wind vector pointing from the sea (normal direction) toward land
    toward_land = -(u * math.sin(n) + v * math.cos(n))
    return -toward_land


def layer_min(levels: list[tuple[float, float, float]], pad: dict, top_m: float, min_top: float, alpha: float,
              full: bool = False):
    """levels: [(height_m, u, v)]. Returns (most negative component, height) or (None, None);
    with full=True returns (component, height, u, v) of that level."""
    lv = sorted((z, u, v) for z, u, v in levels if z is not None and z > 1 and u is not None and v is not None)
    none = (None, None, None, None) if full else (None, None)
    if not lv:
        return none
    below = [x for x in lv if x[0] <= top_m]
    above = [x for x in lv if x[0] > top_m]
    pts = list(below)
    if below and above:
        (z1, u1, v1), (z2, u2, v2) = below[-1], above[0]
        w = math.log(top_m / z1) / math.log(z2 / z1)
        pts.append((top_m, u1 + w * (u2 - u1), v1 + w * (v2 - v1)))
    elif below and below[-1][0] >= min_top:
        z, u, v = below[-1]
        k = (top_m / z) ** alpha
        pts.append((top_m, u * k, v * k))
    else:
        return none                            # highest level too low to speak for the layer
    best = min(pts, key=lambda p: onshore(p[1], p[2], pad["coast_normal_deg"]))
    c = onshore(best[1], best[2], pad["coast_normal_deg"])
    return (c, best[0], best[1], best[2]) if full else (c, best[0])


# ---------------------------------------------------------------- GRIB (several points per file)
def _decode(blob, points):
    out = []
    fd, path = tempfile.mkstemp(suffix=".grib2")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(blob)
        with open(path, "rb") as f:
            while True:
                gid = grib.eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                try:
                    miss = grib.eccodes.codes_get(gid, "missingValue")
                    vals = []
                    for lat, lon in points:
                        v = grib._nearest(gid, lat, lon)
                        vals.append(None if abs(v - miss) < 1e-6 or abs(v) > 1e10 else v)
                    out.append((vals, grib._grid_meta(gid)))
                finally:
                    grib.eccodes.codes_release(gid)
    finally:
        os.unlink(path)
    return out


def _multipart(resp) -> list[tuple[int, bytes]]:
    """[(start, bytes)] from a multipart/byteranges (or single-range) response."""
    ctype = resp.headers.get("Content-Type", "")
    if "multipart/byteranges" not in ctype:
        m = re.search(r"bytes (\d+)-", resp.headers.get("Content-Range", ""))
        return [(int(m.group(1)) if m else 0, resp.content)]
    boundary = re.search(r"boundary=\"?([^\";]+)\"?", ctype).group(1).encode()
    parts = []
    for chunk in resp.content.split(b"--" + boundary):
        if b"\r\n\r\n" not in chunk:
            continue
        head, body = chunk.split(b"\r\n\r\n", 1)
        m = re.search(rb"Content-Range:\s*bytes (\d+)-(\d+)", head, re.I)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            parts.append((a, body[: b - a + 1]))
    return parts


def fetch_points(url: str, points: list) -> tuple[dict, dict]:
    if grib.eccodes is None:
        raise RuntimeError("eccodes not installed")
    inv = read_idx(url + ".idx")
    found = sorted(match_fields(inv, FIELDS).items(), key=lambda kv: kv[1][1])
    groups: list[list] = []
    for name, rec in found:
        whole = "." not in rec[0]
        g = groups[-1] if groups else None
        if g and whole and g[-1][2] and g[-1][1][2] is not None and rec[1] == g[-1][1][2] + 1:
            g.append((name, rec, whole))
        else:
            groups.append([(name, rec, whole)])
    def rng(g):
        a, b = g[0][1][1], g[-1][1][2]
        return f"{a}-{b}" if b is not None else f"{a}-"

    blobs = {}
    # NOMADS accepts several byte ranges in one request, so each file costs one hit there.
    if "nomads.ncep.noaa.gov" in url and len(groups) > 1 and all(g[-1][1][2] is not None for g in groups):
        r = grib._get(url, headers={"Range": "bytes=" + ",".join(rng(g) for g in groups)}, timeout=90)
        if r.status_code in (403, 404, 416):
            raise Missing(url)
        r.raise_for_status()
        if r.status_code == 206:
            for a, body in _multipart(r):
                blobs[a] = body
    vals, meta = {}, {}
    for g in groups:
        blob = blobs.get(g[0][1][1])
        if blob is None:
            r = grib._get(url, headers={"Range": "bytes=" + rng(g)}, timeout=60)
            if r.status_code in (403, 404, 416):
                raise Missing(url)
            r.raise_for_status()
            blob = r.content
        dec = _decode(blob, points)
        if len(g) == 1 and not g[0][2]:
            k = int(g[0][1][0].split(".")[1]) - 1
            pairs = [(g[0][0], dec[k] if k < len(dec) else ([None] * len(points), {}))]
        else:
            pairs = [(n, dec[i] if i < len(dec) else ([None] * len(points), {})) for i, (n, _, _) in enumerate(g)]
        for n, (v, m) in pairs:
            vals[n], meta[n] = v, m
    return vals, meta


def _fmt(tmpl, cycle, fh):
    g = time.gmtime(cycle)
    return tmpl.format(ymd=time.strftime("%Y%m%d", g), hh=f"{g.tm_hour:02d}", fh=fh)


def _levels_from_grib(urls, cycle, ctx) -> dict | None:
    """{pad_id: [[z, u, v], ...]} merged from one or more files for the same valid time."""
    key = "|".join(urls)
    cached = ctx.cache.get(key)
    if cached is not None:
        return cached
    if time.time() > ctx.deadline:          # out of time this run; picked up next hour from where it left off
        ctx.skipped += 1
        return None
    points = [(p["lat"], p["lon"]) for p in ctx.pads]
    merged, got_any = {}, False
    for url in urls:
        try:
            v, m = fetch_points(url, points)
        except Missing:
            if url == urls[0]:
                return None
            continue
        got_any = True
        merged.update({k: (v[k], m[k]) for k in v})
    if not got_any:
        return None
    out = {}
    for i, pad in enumerate(ctx.pads):
        lv = []
        def uv(ku, kv):
            if ku not in merged or kv not in merged:
                return None
            u, v = merged[ku][0][i], merged[kv][0][i]
            if u is None or v is None:
                return None
            return earth_relative(u, v, merged[ku][1], pad["lon"])
        for z, ku, kv in ((10.0, "u10", "v10"), (80.0, "u80", "v80")):
            r = uv(ku, kv)
            if r:
                lv.append([z, round(r[0], 2), round(r[1], 2)])
        for p in PLEVS:
            r = uv(f"u{p}", f"v{p}")
            h = merged.get(f"h{p}", ([None] * len(points), {}))[0][i]
            if r and h is not None and h > 15:          # height above sea level ~ above ground at the Cape
                lv.append([round(h, 1), round(r[0], 2), round(r[1], 2)])
        out[pad["id"]] = lv
    ctx.cache.put(key, cycle, out)
    return out


def _records(levels_by_pad, ctx):
    rec = {}
    for pad in ctx.pads:
        lv = levels_by_pad.get(pad["id"]) or []
        c, z, u, v = layer_min([tuple(x) for x in lv], pad, ctx.top_m, ctx.min_top, ctx.alpha, full=True)
        if c is not None:
            _put(rec, pad["id"], c, z, u, v)
    return rec or None


def _put(rec, pid, c, z, u, v):
    """Component, its height, and the wind (speed, direction) at that height for the polar plots."""
    s, d = sd_from_uv(u, v)
    rec[f"c@{pid}"] = round(c, 2)
    rec[f"z@{pid}"] = int(round(z))
    rec[f"s@{pid}"] = round(s, 1)
    rec[f"d@{pid}"] = int(round(d)) % 360


def _pick(bases, tmpl, cycles):
    for base in bases:
        for c in cycles:
            if exists(f"{base}/{_fmt(tmpl, c, 1)}.idx"):
                return base, c
    return None


def _run(tasks, ctx, workers, sid):
    # newest cycles and shortest lead times first, so a run cut short by the time budget
    # still has the most useful members
    tasks = sorted(tasks, key=lambda t: (-t[2], t[3]))

    def run(task):
        mid, urls, c, valid = task
        try:
            lv = _levels_from_grib(urls, c, ctx)
            return task, (_records(lv, ctx) if lv else None)
        except Exception as e:
            log.warning("%s %s: %s", sid, urls[0].rsplit("/", 1)[-1], e)
            return task, None
    members = {}
    with ThreadPoolExecutor(workers) as ex:
        for (mid, _, _, valid), rec in ex.map(run, tasks):
            if rec:
                members.setdefault(mid, {})[valid] = rec
    return members


def src_hrrr(scfg, ctx):
    recent = [floor_hour(ctx.now) - k * 3600 for k in range(10)]
    base = scfg.get("base", HRRR_BASE)
    picked = _pick([base], HRRR_SFC, recent)
    if not picked:
        return SourceResult({}, status="missing", note="no recent cycle found")
    _, latest = picked
    cycles = [latest - k * 3600 for k in range(int(scfg.get("lag_cycles", 6)))]
    cycles += [c for c in (latest - k * 3600 for k in range(30))
               if time.gmtime(c).tm_hour % 6 == 0 and c not in cycles][: int(scfg.get("synoptic_extra", 2))]
    tasks = []
    for c in cycles:
        mx = 48 if time.gmtime(c).tm_hour % 6 == 0 else 18
        for fh in range(0, mx + 1):
            if ctx.in_window(c + fh * 3600):
                tasks.append((time.strftime("%d/%HZ", time.gmtime(c)),
                              [f"{base}/{_fmt(HRRR_SFC, c, fh)}", f"{base}/{_fmt(HRRR_PRS, c, fh)}"], c, c + fh * 3600))
    mem = _run(tasks, ctx, 16, scfg["id"])
    return SourceResult(mem, cycle=iso(latest), note=f"{len(mem)}/{len(cycles)} cycles; 10 m, 80 m and 1000-950 mb",
                        status="ok" if len(mem) == len(cycles) else ("partial" if mem else "missing"))


def src_multi(scfg, ctx):
    tasks, notes, missing = [], [], []
    for comp in scfg["components"]:
        hours = set(comp.get("cycles", [0, 6, 12, 18]))
        cands = [c for c in (floor_hour(ctx.now) - k * 3600 for k in range(60)) if time.gmtime(c).tm_hour in hours]
        base = comp.get("base") or comp["bases"][0]
        tmpl = comp["file"]
        if comp["id"] == "HRRR":                              # HRRR pressure levels live in the prs file
            picked = _pick([base], tmpl, cands)
            extra = [HRRR_PRS]
        else:
            picked = _pick([base], tmpl, cands)
            extra = []
        if not picked:
            missing.append(comp["id"])
            continue
        _, latest = picked
        use = [c for c in cands if c <= latest][: int(comp.get("lag", 2))]
        notes.append(f"{comp['id']} " + ", ".join(time.strftime("%HZ", time.gmtime(c)) for c in use))
        for c in use:
            for fh in range(0, int(comp.get("max_fh", 48)) + 1):
                if ctx.in_window(c + fh * 3600):
                    urls = [f"{base}/{_fmt(tmpl, c, fh)}"] + [f"{base}/{_fmt(x, c, fh)}" for x in extra]
                    tasks.append((f"{comp['id']} {time.strftime('%d/%HZ', time.gmtime(c))}", urls, c, c + fh * 3600))
    mem = _run(tasks, ctx, int(scfg.get("workers", 4)), scfg["id"])
    note = "; ".join(notes + ([f"missing: {', '.join(missing)}"] if missing else []))
    return SourceResult(mem, note=note, status="missing" if not mem else ("partial" if missing else "ok"))


# ---------------------------------------------------------------- Open-Meteo
ENS_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"
DET_URL = "https://api.open-meteo.com/v1/forecast"
HLEVS = (10, 80, 100, 120, 180)
_OMKEY = re.compile(r"^(wind_speed|wind_direction|geopotential_height)_(\d+)(m|hPa)(?:_member(\d+))?$")


def _om_vars():
    v = []
    for z in HLEVS:
        v += [f"wind_speed_{z}m", f"wind_direction_{z}m"]
    for p in PLEVS:
        v += [f"wind_speed_{p}hPa", f"wind_direction_{p}hPa", f"geopotential_height_{p}hPa"]
    return v


def _om_fetch(url, model, ctx):
    variables = _om_vars()
    days = max(1, math.ceil((ctx.t_end - ctx.now) / 86400) + 1)
    r = None
    for _ in range(len(variables)):
        r = SESSION.get(url, timeout=120, params={
            "latitude": ",".join(str(p["lat"]) for p in ctx.pads), "longitude": ",".join(str(p["lon"]) for p in ctx.pads),
            "models": model, "hourly": ",".join(variables), "wind_speed_unit": "ms",
            "timeformat": "unixtime", "timezone": "GMT", "past_days": 1, "forecast_days": min(days, 16)})
        if r.status_code != 400:
            break
        bad = [v for v in variables if v in r.text]
        if not bad:
            break
        variables.remove(bad[0])
    r.raise_for_status()
    js = r.json()
    return js if isinstance(js, list) else [js]


def _om_members(locs, ctx):
    out = {}
    for pad, loc in zip(ctx.pads, locs):
        h = loc["hourly"]
        times = h["time"]
        per = {}
        for key, arr in h.items():
            m = _OMKEY.match(key)
            if not m or not any(x is not None for x in arr):
                continue
            mid = f"m{int(m.group(4)):02d}" if m.group(4) else "m00"
            per.setdefault(mid, {})[(m.group(1), int(m.group(2)), m.group(3))] = arr
        for mid, f in per.items():
            for i, t in enumerate(times):
                if not ctx.in_window(t):
                    continue
                lv = []
                for z in HLEVS:
                    s, d = f.get(("wind_speed", z, "m")), f.get(("wind_direction", z, "m"))
                    if s and d and s[i] is not None and d[i] is not None:
                        lv.append((float(z), *uv_from_sd(s[i], d[i])))
                for p in PLEVS:
                    s, d, gh = (f.get(("wind_speed", p, "hPa")), f.get(("wind_direction", p, "hPa")),
                                f.get(("geopotential_height", p, "hPa")))
                    if s and d and gh and None not in (s[i], d[i], gh[i]) and gh[i] > 15:
                        lv.append((float(gh[i]), *uv_from_sd(s[i], d[i])))
                c, z, u, v = layer_min(lv, pad, ctx.top_m, ctx.min_top, ctx.alpha, full=True)
                if c is not None:
                    _put(out.setdefault(mid, {}).setdefault(int(t), {}), pad["id"], c, z, u, v)
    return out


def src_openmeteo(scfg, ctx, url):
    path = os.path.join(ROOT, "cache", f"onshore_om_{scfg['id']}.json")
    try:
        with open(path) as f:
            old = json.load(f)
        if ctx.now - old["fetched"] < ctx.om_refresh_h * 3600:
            mem = {m: {int(t): r for t, r in s.items() if ctx.in_window(int(t))} for m, s in old["members"].items()}
            return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(old["fetched"])) + " (cached)",
                                note=f"{len(mem)} members via Open-Meteo", status="ok" if mem else "missing")
    except (OSError, ValueError, KeyError):
        pass
    mem = _om_members(_om_fetch(url, scfg["model"], ctx), ctx)
    if mem:
        with open(path, "w") as f:
            json.dump({"fetched": ctx.now, "members": mem}, f, separators=(",", ":"))
    return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(ctx.now)),
                        note=f"{len(mem)} members via Open-Meteo", status="ok" if mem else "missing")


KINDS = {
    "hrrr": src_hrrr,
    "multi_model": src_multi,
    "openmeteo_ens": lambda s, c: src_openmeteo(s, c, ENS_URL),
    "openmeteo_det": lambda s, c: src_openmeteo(s, c, DET_URL),
}


# ---------------------------------------------------------------- build
def build(cfg, only=None):
    oc = cfg["onshore"]
    now = int(time.time())
    t_start = floor_hour(now) - int(cfg["window"]["hours_back"]) * 3600
    t_end = floor_hour(now) + int(cfg["window"]["hours_ahead"]) * 3600
    timeline = list(range(t_start, t_end + 1, 3600))
    idx = {t: i for i, t in enumerate(timeline)}
    cache = PointCache(os.path.join(ROOT, "cache", "onshore.json"))
    cache.prune(now - 4 * 86400)
    ctx = Context(now=now, lat=oc["pads"][0]["lat"], lon=oc["pads"][0]["lon"], target_m=0.0, t_start=t_start, t_end=t_end,
                  cache=cache, alpha=float(cfg.get("vertical", {}).get("power_law_alpha", 0.14)))
    ctx.pads = oc["pads"]
    ctx.top_m = float(oc["layer_top_ft"]) * 0.3048
    ctx.min_top = float(oc.get("min_top_level_m", 80))
    ctx.om_refresh_h = float(oc.get("openmeteo_refresh_hours", 3))
    ctx.deadline = time.time() + 60 * float(oc.get("time_budget_min", 20))
    ctx.skipped = 0

    main = {s["id"]: s for s in cfg["sources"]}
    sources = []
    for sid in oc["sources"]:
        scfg = dict(main.get(sid, {}), id=sid)
        if (only and sid not in only) or not main.get(sid, {}).get("enabled", True) and sid not in oc.get("force", []):
            continue
        kind = scfg.get("kind")
        fn = KINDS.get(kind)
        if fn is None:
            log.info("onshore: %s (%s) has no onshore reader; skipped", sid, kind)
            continue
        t0 = time.time()
        try:
            res = fn(scfg, ctx)
        except Exception as e:
            log.exception("onshore %s failed", sid)
            res = SourceResult({}, status="error", note=f"{type(e).__name__}: {e}"[:200])
        members = []
        for mid, series in sorted(res.members.items()):
            v = {}
            for t, rec in series.items():
                i = idx.get(int(t))
                if i is None:
                    continue
                for k, x in rec.items():
                    v.setdefault(k, [None] * len(timeline))[i] = x
            if v:
                members.append({"id": mid, "v": v})
        skipped, ctx.skipped = ctx.skipped, 0
        if skipped:
            res.note += f"; {skipped} files left for the next run (time budget)"
            if res.status == "ok":
                res.status = "partial"
        el = round(time.time() - t0, 1)
        log.info("onshore %-10s %-8s %3d members %6.1fs  %s", sid, res.status, len(members), el, res.note)
        sources.append({"id": sid, "label": scfg.get("label", sid), "family": scfg.get("family", "global"),
                        "weight": float(scfg.get("weight", 1)), "status": res.status, "cycle": res.cycle,
                        "note": res.note, "seconds": el, "members": members})
    cache.save()
    return {"generated": iso(now), "generated_unix": now, "display_tz": cfg["site"]["display_tz"],
            "layer_top_ft": oc["layer_top_ft"], "pads": oc["pads"], "times": timeline, "sources": sources}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(ROOT, "config.yaml"))
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    with open(a.config) as f:
        cfg = yaml.safe_load(f)
    if "onshore" not in cfg:
        log.info("no onshore section in config; nothing to do")
        return 0
    data = build(cfg, set(filter(None, a.only.split(","))) or None)
    if not any(s["members"] for s in data["sources"]):
        log.error("onshore: no members from any source; keeping previous output")
        return 1
    out = os.path.join(ROOT, "docs", "data", "onshore.json")
    with open(out + ".tmp", "w") as f:
        json.dump(data, f, separators=(",", ":"))
    os.replace(out + ".tmp", out)
    log.info("wrote onshore.json (%.0f kB)", os.path.getsize(out) / 1024)
    return 0


if __name__ == "__main__":
    sys.exit(main())
