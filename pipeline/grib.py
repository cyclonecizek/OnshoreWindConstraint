"""Pull single grid-point values out of remote GRIB2 files.

Only the records we need are downloaded (HTTP Range from the .idx inventory),
then eccodes finds the nearest grid point. Lambert-conformal u/v that are
grid-relative (HRRR, RRFS, HREF, REFS) are rotated to earth-relative."""
from __future__ import annotations

import math
import os
import re
import tempfile

from .common import SESSION, log

try:
    import eccodes
    eccodes.codes_grib_multi_support_on()
except ImportError:  # probe.py can still list inventories without it
    eccodes = None

# Standard wind records. First instantaneous match in the inventory wins.
WIND_FIELDS = {
    "u10": r":UGRD:10 m above ground:",
    "v10": r":VGRD:10 m above ground:",
    "u80": r":UGRD:80 m above ground:",
    "v80": r":VGRD:80 m above ground:",
    "w10": r":WIND:10 m above ground:",
    "w80": r":WIND:80 m above ground:",
    "g0": r":GUST:surface:",
}
# Reject time-processed records (hourly max, averages) and probability records.
_SKIP = re.compile(r"(\bmax\b|\bmin\b|\bave\b|\bacc\b|prob|%)", re.I)


class Missing(Exception):
    """File or inventory not published (yet)."""


def exists(url: str) -> bool:
    try:
        r = SESSION.head(url, timeout=15, allow_redirects=True)
        return r.status_code == 200
    except Exception:
        return False


def read_idx(idx_url: str) -> list[tuple[str, int, int | None, str]]:
    """Returns [(record_no, start, end_or_None, description)]."""
    r = SESSION.get(idx_url, timeout=30)
    if r.status_code in (403, 404):
        raise Missing(idx_url)
    r.raise_for_status()
    recs = []
    for line in r.text.splitlines():
        parts = line.split(":")
        if len(parts) < 3:
            continue
        try:
            start = int(parts[1])
        except ValueError:
            continue
        recs.append((parts[0], start, ":" + ":".join(parts[2:])))
    out = []
    for i, (no, start, desc) in enumerate(recs):
        end = None
        for j in range(i + 1, len(recs)):       # sub-messages share an offset
            if recs[j][1] > start:
                end = recs[j][1] - 1
                break
        out.append((no, start, end, desc))
    return out


def match_fields(inv, wanted: dict[str, str]) -> dict[str, tuple]:
    found = {}
    for name, pat in wanted.items():
        rx = re.compile(pat)
        for rec in inv:
            m = rx.search(rec[3])
            if m and not _SKIP.search(rec[3][m.end():]):
                found[name] = rec
                break
    return found


def _nearest(gid, lat, lon):
    for lo in (lon, lon % 360.0):
        try:
            res = eccodes.codes_grib_find_nearest(gid, lat, lo)
            r = res[0]
            val = r["value"] if isinstance(r, dict) else getattr(r, "value")
            return float(val)
        except Exception:
            continue
    raise RuntimeError("nearest-point lookup failed")


def _grid_meta(gid) -> dict:
    meta = {"grid": eccodes.codes_get(gid, "gridType")}
    try:
        meta["rel"] = int(eccodes.codes_get(gid, "uvRelativeToGrid"))
    except Exception:
        meta["rel"] = 0
    if meta["grid"] == "lambert":
        meta["lov"] = float(eccodes.codes_get(gid, "LoVInDegrees"))
        meta["latin1"] = float(eccodes.codes_get(gid, "Latin1InDegrees"))
    return meta


def _decode(blob: bytes, sub_index: int, lat: float, lon: float):
    fd, path = tempfile.mkstemp(suffix=".grib2")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(blob)
        with open(path, "rb") as f:
            k = 0
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                try:
                    k += 1
                    if k == sub_index:
                        val = _nearest(gid, lat, lon)
                        miss = eccodes.codes_get(gid, "missingValue")
                        if abs(val - miss) < 1e-6 or abs(val) > 1e10:
                            return None, _grid_meta(gid)
                        return val, _grid_meta(gid)
                finally:
                    eccodes.codes_release(gid)
    finally:
        os.unlink(path)
    return None, {}


def fetch_point(grib_url: str, wanted: dict[str, str], lat: float, lon: float,
                idx_url: str | None = None) -> tuple[dict, dict]:
    """Returns (values, meta) for each wanted field found. Raises Missing."""
    if eccodes is None:
        raise RuntimeError("eccodes not installed")
    inv = read_idx(idx_url or grib_url + ".idx")
    found = match_fields(inv, wanted)
    vals, meta = {}, {}
    for name, (no, start, end, _desc) in found.items():
        rng = f"bytes={start}-{end}" if end is not None else f"bytes={start}-"
        r = SESSION.get(grib_url, headers={"Range": rng}, timeout=60)
        if r.status_code in (403, 404, 416):
            raise Missing(grib_url)
        r.raise_for_status()
        sub = int(no.split(".")[1]) if "." in no else 1
        v, m = _decode(r.content, sub, lat, lon)
        vals[name], meta[name] = v, m
    return vals, meta


def earth_relative(u: float, v: float, meta: dict, lon: float) -> tuple[float, float] | None:
    """Rotate grid-relative Lambert u/v to earth-relative (NCEP formula)."""
    if not meta or not meta.get("rel"):
        return u, v
    if meta.get("grid") != "lambert":
        log.warning("grid-relative winds on %s grid not handled; dropping", meta.get("grid"))
        return None
    lov = meta["lov"] if meta["lov"] <= 180 else meta["lov"] - 360
    angle = math.sin(math.radians(meta["latin1"])) * math.radians(lon - lov)
    c, s = math.cos(angle), math.sin(angle)
    return c * u + s * v, -s * u + c * v
