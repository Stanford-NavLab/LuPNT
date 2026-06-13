#!/usr/bin/env python3
"""
ANTEXLoader: parse IGS ANTEX (.atx) files and query satellite antenna PCO/PCV.

- Loads and indexes the file at construction.
- Query by epoch, constellation, PRN, and frequency code (or friendly name).
- Returns satellite PCO (N/E/U) in meters and NOAZI PCV in meters.

Notes:
- ANTEX PCO/PCV are in satellite antenna NEU frame (not LOS-projected).
- This class returns the model values; LOS projection requires attitude/yaw model.

Example:
    loader = ANTEXLoader("igs20.atx")
    epoch = "2025-01-15T12:00:00Z"
    pco = loader.get_pco(epoch, "G", 5, freq="L1")
    pcv30 = loader.get_pcv_noazi(epoch, "G", 5, freq="L1", elevation_deg=30.0)
"""

from __future__ import annotations

import datetime as dt
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
import pylupnt as pnt

LABEL_COL = 60  # label starts at col 61 (1-based), typical ANTEX


def _unit(v: np.ndarray, eps: float = 1e-15) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(3)
    n = np.linalg.norm(v)
    if n < eps:
        raise ValueError("Cannot normalize near-zero vector.")
    return v / n


def _get_label(line: str) -> str:
    return line[LABEL_COL:].strip() if len(line) > LABEL_COL else ""


def _get_data(line: str) -> str:
    return line[:LABEL_COL].rstrip("\n")


def _parse_epoch(s: str | dt.datetime) -> dt.datetime:
    """Accept dt.datetime or ISO-like string; return tz-aware UTC datetime."""
    if isinstance(s, dt.datetime):
        t = s
        if t.tzinfo is None:
            return t.replace(tzinfo=dt.timezone.utc)
        return t.astimezone(dt.timezone.utc)

    s = str(s).strip()
    if s.endswith("Z"):
        s = s[:-1].strip()
    fmts = [
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d",
    ]
    for f in fmts:
        try:
            t = dt.datetime.strptime(s, f)
            return t.replace(tzinfo=dt.timezone.utc)
        except ValueError:
            pass
    raise ValueError(f"Could not parse epoch: {s}")


def _parse_antex_datetime(line: str) -> Optional[dt.datetime]:
    """
    ANTEX lines like: ' 2020 01 05 00 00 00.0                 VALID FROM'
    """
    parts = line[:60].split()
    if len(parts) < 3:
        return None
    try:
        y = int(parts[0])
        mo = int(parts[1])
        d = int(parts[2])
        hh = int(parts[3]) if len(parts) > 3 else 0
        mm = int(parts[4]) if len(parts) > 4 else 0
        ss = float(parts[5]) if len(parts) > 5 else 0.0
        sec = int(ss)
        micro = int(round((ss - sec) * 1e6))
        return dt.datetime(y, mo, d, hh, mm, sec, microsecond=micro, tzinfo=dt.timezone.utc)
    except Exception:
        return None


def _in_range(t: dt.datetime, start: Optional[dt.datetime], end: Optional[dt.datetime]) -> bool:
    if start is not None and t < start:
        return False
    if end is not None and t > end:
        return False
    return True


def _sat_id(const: str, prn: int) -> str:
    return f"{const.strip().upper()}{int(prn):02d}"


@dataclass
class FreqPattern:
    freq_code: str  # e.g., "G01", "G02", "E01", "E05", ...
    pco_neu_m: Optional[Tuple[float, float, float]] = None
    # NOAZI PCV values in meters, corresponding to zenith grid
    pcv_noazi_m: List[float] = field(default_factory=list)
    zen1_deg: Optional[float] = None
    zen2_deg: Optional[float] = None
    dzen_deg: Optional[float] = None


@dataclass
class SatAntennaEntry:
    sat_id: str  # e.g., "G05"
    valid_from: Optional[dt.datetime] = None
    valid_until: Optional[dt.datetime] = None
    freqs: Dict[str, FreqPattern] = field(default_factory=dict)


class ANTEXLoader:
    """
    Loader/querier for satellite antenna corrections from an ANTEX file.

    This class focuses on satellite entries (PRN blocks). Receiver antenna
    blocks are ignored by default.
    """

    def __init__(self, atx_path: str):
        self.atx_path = atx_path
        # index: sat_id -> list of entries (different validity intervals)
        self._sat_index: Dict[str, List[SatAntennaEntry]] = {}
        self._parse_antex(atx_path)

    # ---------------- public API ----------------

    def get_pco(
        self,
        epoch: str | dt.datetime,
        constellation: str,
        prn: int,
        freq: str = "L1",
        freqcode: Optional[str] = None,
    ) -> Tuple[float, float, float]:
        """
        Return satellite PCO (North, East, Up) in meters for requested epoch.
        """
        entry, pat = self._get_entry_and_pattern(epoch, constellation, prn, freq, freqcode)
        if pat.pco_neu_m is None:
            raise KeyError(f"PCO not found for {entry.sat_id} freq={pat.freq_code}")
        return pat.pco_neu_m

    def get_pcv_noazi_table(
        self,
        epoch: str | dt.datetime,
        constellation: str,
        prn: int,
        freq: str = "L1",
        freqcode: Optional[str] = None,
    ) -> Tuple[List[float], Tuple[Optional[float], Optional[float], Optional[float]]]:
        """
        Return (pcv_values_m, (zen1, zen2, dzen)) for NOAZI PCV table.
        """
        _, pat = self._get_entry_and_pattern(epoch, constellation, prn, freq, freqcode)
        return pat.pcv_noazi_m.copy(), (pat.zen1_deg, pat.zen2_deg, pat.dzen_deg)

    def get_pcv_noazi(
        self,
        epoch: str | dt.datetime,
        constellation: str,
        prn: int,
        freq: str = "L1",
        freqcode: Optional[str] = None,
        elevation_deg: Optional[float] = None,
    ) -> Optional[float]:
        """
        If elevation_deg is provided, return interpolated NOAZI PCV (meters).
        If elevation_deg is None, return None (use get_pcv_noazi_table instead).
        """
        if elevation_deg is None:
            return None
        _, pat = self._get_entry_and_pattern(epoch, constellation, prn, freq, freqcode)
        return self._interp_noazi(pat, float(elevation_deg))

    def available_freqcodes(
        self, constellation: str, prn: int, epoch: str | dt.datetime
    ) -> List[str]:
        """
        List available ANTEX frequency codes (e.g., G01,G02,...) for this satellite at epoch.
        """
        t = _parse_epoch(epoch)
        sid = _sat_id(constellation, prn)
        entry = self._select_entry(sid, t)
        return sorted(entry.freqs.keys())

    def rot_neu_from_ecef(self, t_tai, r_sat_ecef):
        """
        Return C_ecef_neu (3x3) such that:
            v_ecef = C_ecef_neu @ v_neu

        NEU basis vectors (as columns) expressed in ECEF:
            C = [n_ecef, e_ecef, u_ecef]
        """
        r = np.asarray(r_sat_ecef, dtype=float).reshape(3)

        # z-axis: radial outward (geocentric up)
        kvec = -_unit(r, eps=1e-12)
        # x-axis: projection of ECEF x-axis onto plane perpendicular to z
        t_tdb = pnt.convert_time(t_tai, pnt.TAI, pnt.TDB)
        r_sun = pnt.get_body_pos_vel(t_tdb, pnt.EARTH, pnt.SUN, pnt.ECEF)[:3]
        jvec = _unit(r_sun - r_sat_ecef, eps=1e-12)  # approximate Earth-Sun direction
        # y-axis: cross product
        ivec = np.cross(jvec, kvec)

        # Rotation vector
        Cijk = np.column_stack((ivec, jvec, kvec))  # columns are IJK axes in ECEF

        return Cijk

    # ---------------- internals ----------------
    @staticmethod
    def _freq_to_antex_code(constellation: str, freq: str) -> str:
        """
        Map friendly frequency names to ANTEX codes (common cases).
        If your file differs, pass freqcode explicitly.
        """
        const = constellation.strip().upper()
        f = freq.strip().upper()

        gps = {"L1": "G01", "L2": "G02", "L5": "G05"}
        gal = {"E1": "E01", "E5A": "E05", "E5B": "E07", "E5": "E08", "E6": "E06"}
        glo = {"G1": "R01", "G2": "R02", "G3": "R03"}
        qzs = {"L1": "J01", "L2": "J02", "L5": "J05", "L6": "J06"}
        bds = {"B1I": "C02", "B1C": "C01", "B2A": "C05", "B2I": "C07", "B3I": "C06"}

        if const == "G":
            return gps.get(f, f)
        if const == "E":
            return gal.get(f, f)
        if const == "R":
            return glo.get(f, f)
        if const == "J":
            return qzs.get(f, f)
        if const == "C":
            return bds.get(f, f)
        return f  # fallback: already a code?

    def _get_entry_and_pattern(
        self,
        epoch: str | dt.datetime,
        constellation: str,
        prn: int,
        freq: str,
        freqcode: Optional[str],
    ) -> Tuple[SatAntennaEntry, FreqPattern]:
        t = _parse_epoch(epoch)
        sid = _sat_id(constellation, prn)
        entry = self._select_entry(sid, t)

        code = (
            freqcode.strip().upper() if freqcode else self._freq_to_antex_code(constellation, freq)
        )
        if code not in entry.freqs:
            avail = ", ".join(sorted(entry.freqs.keys()))
            raise KeyError(
                f"Frequency {code} not found for {sid} at {t.isoformat()}. Available: {avail}"
            )
        return entry, entry.freqs[code]

    def _select_entry(self, sat_id: str, t: dt.datetime) -> SatAntennaEntry:
        if sat_id not in self._sat_index:
            raise KeyError(f"Satellite {sat_id} not found in ANTEX: {self.atx_path}")
        candidates = [
            e for e in self._sat_index[sat_id] if _in_range(t, e.valid_from, e.valid_until)
        ]
        if not candidates:
            # helpful error with available validity spans
            spans = []
            for e in self._sat_index[sat_id]:
                spans.append(
                    (
                        e.valid_from.isoformat() if e.valid_from else "N/A",
                        e.valid_until.isoformat() if e.valid_until else "N/A",
                    )
                )
            raise KeyError(
                f"No valid entry for {sat_id} at {t.isoformat()}. Available spans: {spans}"
            )
        candidates.sort(
            key=lambda e: e.valid_from or dt.datetime.min.replace(tzinfo=dt.timezone.utc),
            reverse=True,
        )
        return candidates[0]

    @staticmethod
    def _interp_noazi(pat: FreqPattern, elevation_deg: float) -> Optional[float]:
        """
        Interpolate NOAZI PCV in meters given elevation angle in degrees.
        ANTEX grid is in zenith angle: zen = 90 - elevation.
        """
        if not pat.pcv_noazi_m:
            return None
        if pat.zen1_deg is None or pat.dzen_deg is None:
            return None

        zen = 90.0 - elevation_deg
        z1 = pat.zen1_deg
        dz = pat.dzen_deg
        vals = pat.pcv_noazi_m
        n = len(vals)

        idx = (zen - z1) / dz
        if idx <= 0:
            return vals[0]
        if idx >= n - 1:
            return vals[-1]

        i0 = int(math.floor(idx))
        i1 = i0 + 1
        a = idx - i0
        return (1.0 - a) * vals[i0] + a * vals[i1]

    def _parse_antex(self, path: str) -> None:
        """
        Parse ANTEX and populate self._sat_index.
        Only satellite blocks are indexed (those with sat_id like 'G05').
        """
        cur_entry: Optional[SatAntennaEntry] = None
        cur_freq: Optional[str] = None
        zen_def: Optional[Tuple[float, float, float]] = None
        in_noazi = False

        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                line = raw.rstrip("\n")
                label = _get_label(line)
                data = _get_data(line)

                if label == "START OF ANTENNA":
                    cur_entry = SatAntennaEntry(sat_id="")  # fill later
                    cur_freq = None
                    zen_def = None
                    in_noazi = False
                    continue

                if label == "END OF ANTENNA":
                    if cur_entry and cur_entry.sat_id:
                        self._sat_index.setdefault(cur_entry.sat_id, []).append(cur_entry)
                    cur_entry = None
                    cur_freq = None
                    zen_def = None
                    in_noazi = False
                    continue

                if cur_entry is None:
                    continue

                if label == "TYPE / SERIAL NO":
                    # Detect satellite PRN token like "G05", "E11", etc in the data fields
                    parts = data.split()
                    prn = None
                    for p in parts:
                        if re.fullmatch(r"[A-Z]\d{2}", p):
                            prn = p
                            break
                    # If no PRN token, it's a receiver antenna block -> ignore by clearing sat_id
                    if prn is not None:
                        cur_entry.sat_id = prn
                    else:
                        # mark as non-satellite; we will not index it
                        cur_entry.sat_id = ""
                    continue

                if label == "VALID FROM":
                    cur_entry.valid_from = _parse_antex_datetime(line)
                    continue

                if label == "VALID UNTIL":
                    cur_entry.valid_until = _parse_antex_datetime(line)
                    continue

                if label == "ZEN1 / ZEN2 / DZEN":
                    parts = data.split()
                    if len(parts) >= 3:
                        try:
                            zen_def = (float(parts[0]), float(parts[1]), float(parts[2]))
                        except Exception:
                            zen_def = None
                    continue

                if label == "START OF FREQUENCY":
                    if not cur_entry.sat_id:
                        # receiver block or unknown; ignore its frequency content
                        cur_freq = None
                        continue
                    tok = data.split()
                    if not tok:
                        continue
                    cur_freq = tok[0].strip()
                    cur_entry.freqs.setdefault(cur_freq, FreqPattern(freq_code=cur_freq))
                    if zen_def is not None:
                        pat = cur_entry.freqs[cur_freq]
                        if pat.zen1_deg is None:
                            pat.zen1_deg, pat.zen2_deg, pat.dzen_deg = zen_def
                    in_noazi = False
                    continue

                if label == "END OF FREQUENCY":
                    cur_freq = None
                    in_noazi = False
                    continue

                if not cur_entry.sat_id or cur_freq is None:
                    continue

                if label == "NORTH / EAST / UP":
                    parts = data.split()
                    if len(parts) >= 3:
                        try:
                            # ANTEX stores PCO in mm
                            n = float(parts[0]) / 1000.0
                            e = float(parts[1]) / 1000.0
                            u = float(parts[2]) / 1000.0
                            cur_entry.freqs[cur_freq].pco_neu_m = (n, e, u)
                        except Exception:
                            pass
                    continue

                if label == "NOAZI":
                    in_noazi = True
                    vals = data.split()
                    for v in vals:
                        try:
                            cur_entry.freqs[cur_freq].pcv_noazi_m.append(float(v) / 1000.0)
                        except Exception:
                            pass
                    continue

                # NOAZI can continue on lines with empty label (rare but seen)
                if in_noazi and label == "":
                    vals = data.split()
                    for v in vals:
                        try:
                            cur_entry.freqs[cur_freq].pcv_noazi_m.append(float(v) / 1000.0)
                        except Exception:
                            pass
