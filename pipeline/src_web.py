"""Point sources served as JSON: Open-Meteo (ensembles, NBM), NWS gridpoints (NDFD), meteoblue."""
from __future__ import annotations

import math
import os
import re
from datetime import datetime, timedelta

from .common import SESSION, Context, SourceResult, log, sample_from_levels, uv_from_sd

ENS_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"
DET_URL = "https://api.open-meteo.com/v1/forecast"
LEVELS = (10, 80, 100, 120)
_KEY = re.compile(r"^(wind_speed|wind_direction|wind_gusts)_(\d+)m(?:_member(\d+))?$")


def _om_params(model: str, ctx: Context, levels) -> dict:
    hourly = [f"wind_speed_{z}m" for z in levels] + [f"wind_direction_{z}m" for z in levels]
    hourly.append("wind_gusts_10m")
    days = max(1, math.ceil((ctx.t_end - ctx.now) / 86400) + 1)
    return {
        "latitude": ctx.lat, "longitude": ctx.lon, "models": model,
        "hourly": ",".join(hourly), "wind_speed_unit": "mph",
        "timeformat": "unixtime", "timezone": "GMT",
        "past_days": 1, "forecast_days": min(days, 16),
    }


def _om_fetch(url: str, model: str, ctx: Context) -> dict:
    r = SESSION.get(url, params=_om_params(model, ctx, LEVELS), timeout=60)
    if r.status_code == 400:            # some models reject upper levels: fall back to 10 m
        log.info("%s: retrying with 10 m only (%s)", model, r.text[:120])
        r = SESSION.get(url, params=_om_params(model, ctx, (10,)), timeout=60)
    r.raise_for_status()
    return r.json()


def _om_members(js: dict, ctx: Context) -> dict:
    h = js["hourly"]
    times = h["time"]
    by_member: dict[str, dict] = {}
    for key, arr in h.items():
        m = _KEY.match(key)
        if not m or not any(x is not None for x in arr):
            continue
        var, z, mem = m.group(1), int(m.group(2)), m.group(3)
        mid = f"m{int(mem):02d}" if mem else "m00"
        by_member.setdefault(mid, {})[(var, z)] = arr

    out = {}
    for mid, fields in by_member.items():
        series = {}
        for i, t in enumerate(times):
            if not ctx.in_window(t):
                continue
            levels = {}
            for z in LEVELS:
                s = fields.get(("wind_speed", z))
                d = fields.get(("wind_direction", z))
                if s and d and s[i] is not None and d[i] is not None:
                    levels[float(z)] = uv_from_sd(s[i], d[i])
            if not levels:
                continue
            g = fields.get(("wind_gusts", 10))
            w = fields.get(("wind_speed", 10))
            smp = sample_from_levels(levels, ctx,
                                     gust10=g[i] if g else None,
                                     wind10=w[i] if w else None)
            if smp:
                series[t] = smp
        if series:
            out[mid] = series
    return out


def openmeteo_ens(scfg: dict, ctx: Context) -> SourceResult:
    js = _om_fetch(ENS_URL, scfg["model"], ctx)
    mem = _om_members(js, ctx)
    return SourceResult(mem, cycle="latest", note=f"{len(mem)} members via Open-Meteo",
                        status="ok" if mem else "missing")


def openmeteo_det(scfg: dict, ctx: Context) -> SourceResult:
    js = _om_fetch(DET_URL, scfg["model"], ctx)
    mem = _om_members(js, ctx)
    return SourceResult(mem, cycle="latest", note="via Open-Meteo",
                        status="ok" if mem else "missing")


# ---------------------------------------------------------------- NWS / NDFD
_DUR = re.compile(r"P(?:(\d+)D)?(?:T(?:(\d+)H)?)?")


def _expand(prop: dict) -> tuple[dict, str]:
    out = {}
    for v in prop.get("values", []):
        start, dur = v["validTime"].split("/")
        t0 = datetime.fromisoformat(start)
        m = _DUR.match(dur)
        hours = int(m.group(1) or 0) * 24 + int(m.group(2) or 0) if m else 1
        for k in range(max(hours, 1)):
            out[int((t0 + timedelta(hours=k)).timestamp())] = v["value"]
    return out, prop.get("uom", "")


def _to_mph(x, uom: str):
    if x is None:
        return None
    if "km_h" in uom:
        return x * 0.621371
    if "m_s" in uom:
        return x * 2.2369363
    if "kt" in uom or "knot" in uom:
        return x * 1.15078
    return x


def nws_grid(scfg: dict, ctx: Context) -> SourceResult:
    hdr = {"Accept": "application/geo+json"}
    pt = SESSION.get(f"https://api.weather.gov/points/{ctx.lat:.4f},{ctx.lon:.4f}",
                     headers=hdr, timeout=30)
    pt.raise_for_status()
    grid = SESSION.get(pt.json()["properties"]["forecastGridData"], headers=hdr, timeout=60)
    grid.raise_for_status()
    p = grid.json()["properties"]
    spd, su = _expand(p["windSpeed"])
    dirs, _ = _expand(p["windDirection"])
    gst, gu = _expand(p.get("windGust", {}))
    series = {}
    for t, s in spd.items():
        d = dirs.get(t)
        if not ctx.in_window(t) or s is None or d is None:
            continue
        s = _to_mph(s, su)
        smp = sample_from_levels({10.0: uv_from_sd(s, d)}, ctx,
                                 gust10=_to_mph(gst.get(t), gu), wind10=s)
        if smp:
            series[t] = smp
    upd = p.get("updateTime", "")[:16].replace("T", " ")
    return SourceResult({"m00": series} if series else {}, cycle=upd, note="10 m, scaled to height",
                        status="ok" if series else "missing")


# ---------------------------------------------------------------- meteoblue
def meteoblue(scfg: dict, ctx: Context) -> SourceResult:
    key = os.environ.get("METEOBLUE_API_KEY")
    if not key:
        return SourceResult({}, status="disabled", note="set METEOBLUE_API_KEY to enable")
    r = SESSION.get("https://my.meteoblue.com/packages/basic-1h", timeout=60, params={
        "lat": ctx.lat, "lon": ctx.lon, "apikey": key, "format": "json",
        "windspeed": "mph", "timeformat": "timestamp_utc", "tz": "GMT"})
    r.raise_for_status()
    d1 = r.json()["data_1h"]
    series = {}
    for t, s, d in zip(d1["time"], d1["windspeed"], d1["winddirection"]):
        if isinstance(t, str):
            t = int(datetime.fromisoformat(t.replace(" ", "T") + "+00:00").timestamp())
        if s is None or d is None or not ctx.in_window(t):
            continue
        smp = sample_from_levels({10.0: uv_from_sd(s, d)}, ctx)
        if smp:
            series[t] = smp
    return SourceResult({"m00": series} if series else {}, cycle="latest", note="10 m, scaled to height",
                        status="ok" if series else "missing")
