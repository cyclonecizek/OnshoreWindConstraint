"""Shared pieces: HTTP, wind math, vertical adjustment, point cache, run context."""
from __future__ import annotations

import json
import logging
import math
import os
import threading
import time
from dataclasses import dataclass, field

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

log = logging.getLogger("plume")

MPH_PER_MS = 2.2369362920544
FT_TO_M = 0.3048
UA = "plume-constraint/1.0 (+https://github.com)"


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=3, backoff_factor=1.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET", "HEAD"])
    adapter = HTTPAdapter(max_retries=retry, pool_connections=32, pool_maxsize=32)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers["User-Agent"] = UA
    return s


SESSION = make_session()


# ---------------------------------------------------------------- wind math
def uv_from_sd(speed: float, direction: float) -> tuple[float, float]:
    """Meteorological convention: direction the wind blows FROM."""
    r = math.radians(direction)
    return -speed * math.sin(r), -speed * math.cos(r)


def sd_from_uv(u: float, v: float) -> tuple[float, float]:
    speed = math.hypot(u, v)
    direction = (math.degrees(math.atan2(-u, -v)) + 360.0) % 360.0
    return speed, direction


def _ok(x) -> bool:
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def vector_to_height(levels: dict[float, tuple[float, float]], target_m: float,
                     alpha: float) -> tuple[float, float] | None:
    """levels: {height_m: (u, v)}. Returns (u, v) at target_m.

    Bracketed: linear in ln(z) per component. One-sided: nearest level,
    power-law scaled (direction held)."""
    lv = {z: uv for z, uv in levels.items() if uv and _ok(uv[0]) and _ok(uv[1])}
    if not lv:
        return None
    if target_m in lv:
        return lv[target_m]
    below = [z for z in lv if z < target_m]
    above = [z for z in lv if z > target_m]
    if below and above:
        z1, z2 = max(below), min(above)
        w = math.log(target_m / z1) / math.log(z2 / z1)
        (u1, v1), (u2, v2) = lv[z1], lv[z2]
        return u1 + w * (u2 - u1), v1 + w * (v2 - v1)
    z = min(lv, key=lambda zz: abs(math.log(zz / target_m)))
    k = (target_m / z) ** alpha
    return lv[z][0] * k, lv[z][1] * k


def scalar_to_height(levels: dict[float, float], target_m: float, alpha: float) -> float | None:
    lv = {z: s for z, s in levels.items() if _ok(s)}
    if not lv:
        return None
    if target_m in lv:
        return lv[target_m]
    below = [z for z in lv if z < target_m]
    above = [z for z in lv if z > target_m]
    if below and above:
        z1, z2 = max(below), min(above)
        w = math.log(target_m / z1) / math.log(z2 / z1)
        return lv[z1] + w * (lv[z2] - lv[z1])
    z = min(lv, key=lambda zz: abs(math.log(zz / target_m)))
    return lv[z] * (target_m / z) ** alpha


def gust_at_height(gust10, wind10, speed_h):
    """Hold the turbulent gust excess measured at 10 m constant with height."""
    if not (_ok(gust10) and _ok(wind10) and _ok(speed_h)):
        return None
    return speed_h + max(0.0, gust10 - wind10)


def sample_from_levels(uv_levels: dict[float, tuple[float, float]], ctx: "Context",
                       gust10=None, wind10=None):
    """Turn per-level u/v (mph) into (speed_h, dir_h, gust_h) at the evaluation height."""
    uv = vector_to_height(uv_levels, ctx.target_m, ctx.alpha)
    if uv is None:
        return None
    s, d = sd_from_uv(*uv)
    if wind10 is None and 10.0 in uv_levels and uv_levels[10.0]:
        u10, v10 = uv_levels[10.0]
        if _ok(u10) and _ok(v10):
            wind10 = math.hypot(u10, v10)
    return s, d, gust_at_height(gust10, wind10, s)


# ---------------------------------------------------------------- cache
class PointCache:
    """Point values already pulled from GRIB, keyed by file URL.

    Committed with the repo so each hourly run only downloads new cycles."""

    def __init__(self, path: str):
        self.path = path
        self.lock = threading.Lock()
        self.data: dict = {}
        if os.path.exists(path):
            try:
                with open(path) as f:
                    self.data = json.load(f)
            except (OSError, ValueError):
                log.warning("cache unreadable, starting fresh")
        self.hits = 0

    def get(self, key: str):
        with self.lock:
            rec = self.data.get(key)
            if rec is not None:
                self.hits += 1
                return rec["v"]
        return None

    def put(self, key: str, cycle_ts: int, value: dict):
        with self.lock:
            self.data[key] = {"c": cycle_ts, "v": value}

    def prune(self, min_cycle_ts: int):
        with self.lock:
            self.data = {k: r for k, r in self.data.items() if r.get("c", 0) >= min_cycle_ts}

    def save(self):
        tmp = self.path + ".tmp"
        with self.lock, open(tmp, "w") as f:
            json.dump(self.data, f, separators=(",", ":"))
        os.replace(tmp, self.path)


# ---------------------------------------------------------------- context
@dataclass
class Context:
    now: int                      # unix seconds
    lat: float
    lon: float
    target_m: float
    alpha: float
    t_start: int                  # first hour on the timeline
    t_end: int                    # last hour on the timeline
    cache: PointCache
    notes: list = field(default_factory=list)

    def in_window(self, t: int) -> bool:
        return self.t_start <= t <= self.t_end


@dataclass
class SourceResult:
    """members: {member_id: {valid_unix: (speed_mph, dir_deg, gust_mph|None)}}"""
    members: dict
    cycle: str = ""
    note: str = ""
    status: str = "ok"


def floor_hour(ts: float) -> int:
    return int(ts // 3600 * 3600)


def iso(ts: int) -> str:
    return time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime(ts))
