"""
Convert a CCSDS OEM (Orbit Ephemeris Message) file with Cartesian state vectors
to sky coordinates (RA/Dec and Alt/Az) for telescope pointing.

Input OEM assumptions (matches your file):
- TIME_SYSTEM = UTC
- REF_FRAME = EME2000 (treated as close to GCRS/J2000 inertial for pointing purposes)
- Rows are: epoch_iso_utc X Y Z VX VY VZ  (km, km/s)

What this script outputs for a requested time:
- Interpolated state (km, km/s)
- ICRS RA/Dec (J2000-like)
- Topocentric Alt/Az for a given observer location
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, GCRS, ICRS, SkyCoord
from astropy.time import Time


@dataclass(frozen=True)
class OemRow:
    t: Time  # UTC time
    r_km: np.ndarray  # shape (3,)
    v_km_s: np.ndarray  # shape (3,)


def _parse_object_name(lines: Iterable[str]) -> str:
    """
    Extract OBJECT_NAME from OEM metadata for naming outputs.
    Returns empty string if not found.
    """
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("COMMENT"):
            continue
        if line.startswith("OBJECT_NAME"):
            # OEM uses "KEY = VALUE"
            parts = line.split("=", 1)
            if len(parts) == 2:
                return parts[1].strip()
    return ""


def _safe_filename_component(s: str) -> str:
    # Keep alphanumerics and a few safe separators.
    s = s.strip()
    if not s:
        return "OBJECT"
    out = []
    for ch in s:
        if ch.isalnum() or ch in ("-", "_", "."):
            out.append(ch)
        elif ch.isspace():
            out.append("_")
        # else drop
    cleaned = "".join(out).strip("._-")
    return cleaned or "OBJECT"


def _parse_oem_rows(lines: Iterable[str]) -> List[OemRow]:
    """
    Parse OEM rows from a CCSDS OEM text file.
    We ignore metadata and read numeric rows of 7 columns.
    """
    rows: List[OemRow] = []
    in_data = False
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("META_STOP"):
            in_data = True
            continue
        if not in_data:
            continue
        if line.startswith("COMMENT"):
            continue

        parts = line.split()
        if len(parts) != 7:
            # Some OEMs can include covariance blocks etc; skip anything that's not a state row.
            continue

        epoch = parts[0]
        try:
            x, y, z, vx, vy, vz = map(float, parts[1:])
        except ValueError:
            continue

        t = Time(epoch, format="isot", scale="utc")
        r = np.array([x, y, z], dtype=float)
        v = np.array([vx, vy, vz], dtype=float)
        rows.append(OemRow(t=t, r_km=r, v_km_s=v))

    if len(rows) < 2:
        raise ValueError("Parsed fewer than 2 state rows. Is this a CCSDS OEM with data after META_STOP?")

    rows.sort(key=lambda rr: rr.t.utc.jd)
    return rows


def _hermite_interpolate(
    t0: Time,
    r0: np.ndarray,
    v0: np.ndarray,
    t1: Time,
    r1: np.ndarray,
    v1: np.ndarray,
    t: Time,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Cubic Hermite interpolation using position and velocity at endpoints.
    Returns (r_km, v_km_s) at time t.
    """
    dt_total = (t1.utc.jd - t0.utc.jd) * 86400.0  # seconds
    if dt_total == 0.0:
        return r0.copy(), v0.copy()

    tau = ((t.utc.jd - t0.utc.jd) * 86400.0) / dt_total  # [0,1]
    tau = float(tau)

    # Basis polynomials
    h00 = 2 * tau**3 - 3 * tau**2 + 1
    h10 = tau**3 - 2 * tau**2 + tau
    h01 = -2 * tau**3 + 3 * tau**2
    h11 = tau**3 - tau**2

    r = h00 * r0 + h10 * (dt_total * v0) + h01 * r1 + h11 * (dt_total * v1)

    # Derivatives for velocity
    dh00 = 6 * tau**2 - 6 * tau
    dh10 = 3 * tau**2 - 4 * tau + 1
    dh01 = -6 * tau**2 + 6 * tau
    dh11 = 3 * tau**2 - 2 * tau
    v = (dh00 * r0 + dh10 * (dt_total * v0) + dh01 * r1 + dh11 * (dt_total * v1)) / dt_total

    return r, v


def interpolate_state(rows: List[OemRow], t: Time) -> Tuple[np.ndarray, np.ndarray]:
    """
    Interpolate a state at time t using cubic Hermite interpolation between bracketing rows.
    """
    jds = np.array([rr.t.utc.jd for rr in rows], dtype=float)
    tjd = float(t.utc.jd)

    if tjd < jds[0] or tjd > jds[-1]:
        raise ValueError(
            f"Requested time {t.utc.isot} is outside OEM range {rows[0].t.utc.isot} .. {rows[-1].t.utc.isot}"
        )

    i1 = int(np.searchsorted(jds, tjd, side="right"))
    if i1 == 0:
        i1 = 1
    if i1 >= len(rows):
        i1 = len(rows) - 1
    i0 = i1 - 1

    r, v = _hermite_interpolate(
        rows[i0].t, rows[i0].r_km, rows[i0].v_km_s, rows[i1].t, rows[i1].r_km, rows[i1].v_km_s, t
    )
    return r, v


def compute_pointing(
    r_km: np.ndarray,
    t: Time,
    location: EarthLocation,
) -> Tuple[SkyCoord, SkyCoord]:
    """
    Return (icrs, altaz) SkyCoord for the target given geocentric inertial position vector.
    """
    # Treat EME2000 as (approximately) GCRS for pointing use.
    gcrs = GCRS(
        obstime=t,
        obsgeoloc=location.get_gcrs_posvel(t)[0],
        obsgeovel=location.get_gcrs_posvel(t)[1],
    )

    # Represent target by its geocentric position in the same inertial frame, with distance.
    target_gcrs = SkyCoord(
        x=r_km[0] * u.km,
        y=r_km[1] * u.km,
        z=r_km[2] * u.km,
        frame=GCRS(obstime=t),
        representation_type="cartesian",
    )

    icrs = target_gcrs.transform_to(ICRS())
    altaz = target_gcrs.transform_to(AltAz(obstime=t, location=location))
    return icrs, altaz


def _parse_time_utc(s: str) -> Time:
    # Accept "YYYY-MM-DDTHH:MM:SS[.sss]" (OEM-style) or "YYYY-MM-DD HH:MM:SS"
    s2 = s.strip().replace(" ", "T")
    # Force 'Z' handling if user includes it
    if s2.endswith("Z"):
        s2 = s2[:-1]
    return Time(s2, format="isot", scale="utc")


def main() -> int:
    p = argparse.ArgumentParser(description="Convert OEM states to RA/Dec and Alt/Az for telescope pointing.")
    p.add_argument("--oem", required=True, help="Path to CCSDS OEM file (.asc)")
    p.add_argument("--time", required=True, help="UTC time, e.g. 2026-04-02T07:10:00.000")
    p.add_argument("--lat", type=float, required=False, help="Observer latitude in degrees (+N). Optional if only RA/Dec needed.")
    p.add_argument("--lon", type=float, required=False, help="Observer longitude in degrees (+E). Optional if only RA/Dec needed.")
    p.add_argument("--height-m", type=float, default=0.0, help="Observer height above ellipsoid in meters (optional)")
    p.add_argument("--interval", type=float, default=0.0, help="If >0, minutes between positions (time series mode)")
    p.add_argument("--positions", type=int, default=1, help="Number of positions to output (time series mode)")
    p.add_argument(
        "--output-to-csv",
        dest="output_to_csv",
        default="",
        help="If set, write positions to a CSV file (one position per line). Provide an output directory.",
    )
    p.add_argument(
        "--no-altaz",
        action="store_true",
        help="Suppress Alt/Az. Note: Alt/Az is only computed when --lat and --lon are provided.",
    )
    args = p.parse_args()

    oem_path = Path(args.oem)
    text = oem_path.read_text(encoding="utf-8", errors="replace").splitlines()
    rows = _parse_oem_rows(text)
    object_name = _parse_object_name(text)
    object_name_safe = _safe_filename_component(object_name)

    have_location = args.lat is not None and args.lon is not None
    loc: EarthLocation | None = None
    if have_location:
        loc = EarthLocation.from_geodetic(lon=args.lon * u.deg, lat=args.lat * u.deg, height=args.height_m * u.m)
    else:
        # If no location is provided, implicitly suppress Alt/Az.
        args.no_altaz = True
    t0 = _parse_time_utc(args.time)

    # Time series parameters
    n_positions = max(1, int(args.positions))
    step_min = max(0.0, float(args.interval))
    step_sec = step_min * 60.0

    csv_writer: csv.DictWriter | None = None
    csv_file = None
    if args.output_to_csv:
        requested = Path(args.output_to_csv)
        # Treat the argument as an output directory if it is a directory or doesn't end with ".csv".
        out_dir = requested if requested.exists() and requested.is_dir() else requested
        if out_dir.suffix.lower() == ".csv":
            out_dir = out_dir.parent
        out_dir.mkdir(parents=True, exist_ok=True)

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        csv_path = out_dir / f"{object_name_safe}_{stamp}.csv"
        csv_file = csv_path.open("w", newline="", encoding="utf-8")
        fieldnames = ["time_utc", "ra_icrs", "dec_icrs"]
        if have_location and not args.no_altaz:
            fieldnames += ["alt_deg", "az_deg"]
        csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        csv_writer.writeheader()

    wrote_any = False
    for i in range(n_positions):
        ti = Time(t0.utc.jd + (i * step_sec) / 86400.0, format="jd", scale="utc")
        r_km, v_km_s = interpolate_state(rows, ti)
        # Build geocentric inertial coordinate and convert to ICRS (RA/Dec)
        target_gcrs = SkyCoord(
            x=r_km[0] * u.km,
            y=r_km[1] * u.km,
            z=r_km[2] * u.km,
            frame=GCRS(obstime=ti),
            representation_type="cartesian",
        )
        icrs = target_gcrs.transform_to(ICRS())
        ra = icrs.ra.to(u.hourangle)
        dec = icrs.dec.to(u.deg)

        ra_s = ra.to_string(sep=":", precision=2, pad=True)
        dec_s = dec.to_string(sep=":", precision=1, alwayssign=True, pad=True)

        if csv_writer is not None:
            row = {"time_utc": ti.utc.isot, "ra_icrs": ra_s, "dec_icrs": dec_s}
            if have_location and not args.no_altaz:
                altaz = target_gcrs.transform_to(AltAz(obstime=ti, location=loc))
                row["alt_deg"] = f"{altaz.alt.to_value(u.deg):.6f}"
                row["az_deg"] = f"{altaz.az.to_value(u.deg):.6f}"
            csv_writer.writerow(row)
            wrote_any = True
        else:
            print(f"Time (UTC): {ti.utc.isot}")
            print(f"  RA  (ICRS)  = {ra_s}")
            print(f"  Dec (ICRS)  = {dec_s}")
            if have_location and not args.no_altaz:
                altaz = target_gcrs.transform_to(AltAz(obstime=ti, location=loc))
                az = altaz.az.to(u.deg)
                alt = altaz.alt.to(u.deg)
                print(f"  Az          = {az.to_string(sep=':', precision=1, pad=True)}")
                print(f"  Alt         = {alt.to_string(sep=':', precision=1, alwayssign=True, pad=True)}")
            print(f"  r [km]      = {r_km[0]: .6f} {r_km[1]: .6f} {r_km[2]: .6f}")
            print(f"  v [km/s]    = {v_km_s[0]: .9f} {v_km_s[1]: .9f} {v_km_s[2]: .9f}")
            if i != n_positions - 1:
                print()

    if csv_file is not None:
        csv_file.close()
        if wrote_any:
            print(f"Wrote CSV: {csv_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

