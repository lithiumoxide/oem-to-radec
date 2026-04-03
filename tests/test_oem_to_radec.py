import re
from pathlib import Path

import numpy as np
import pytest
from astropy.time import Time

import oem_to_radec as mod


SAMPLE_OEM = """CCSDS_OEM_VERS = 2.0
COMMENT Orion/Planning
CREATION_DATE = 2026-04-02T14:06:23
ORIGINATOR = NASA/JSC/FOD/FDO

META_START
OBJECT_NAME = EM2
OBJECT_ID = 24
CENTER_NAME = EARTH
REF_FRAME = EME2000
TIME_SYSTEM = UTC
START_TIME = 2026-04-02T03:07:49.583
STOP_TIME = 2026-04-02T03:09:34.583
META_STOP

2026-04-02T03:07:49.583 1 2 3 0.1 0.2 0.3
2026-04-02T03:09:34.583 4 6 8 0.4 0.6 0.8
"""


def test_parse_object_name():
    name = mod._parse_object_name(SAMPLE_OEM.splitlines())
    assert name == "EM2"


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("EM2", "EM2"),
        (" Artemis II ", "Artemis_II"),
        ("A/B\\C:D*E?F\"G<H>I|J", "ABCDEFGHIJ"),
        ("", "OBJECT"),
        ("---", "OBJECT"),
    ],
)
def test_safe_filename_component(raw, expected):
    assert mod._safe_filename_component(raw) == expected


def test_parse_oem_rows_reads_two_rows():
    rows = mod._parse_oem_rows(SAMPLE_OEM.splitlines())
    assert len(rows) == 2
    assert rows[0].t.utc.isot == "2026-04-02T03:07:49.583"
    assert np.allclose(rows[0].r_km, [1.0, 2.0, 3.0])
    assert np.allclose(rows[0].v_km_s, [0.1, 0.2, 0.3])


def test_interpolate_state_at_endpoints_matches_exact_rows():
    rows = mod._parse_oem_rows(SAMPLE_OEM.splitlines())
    t0 = rows[0].t
    t1 = rows[1].t

    r0, v0 = mod.interpolate_state(rows, t0)
    r1, v1 = mod.interpolate_state(rows, t1)

    assert np.allclose(r0, rows[0].r_km)
    assert np.allclose(v0, rows[0].v_km_s)
    assert np.allclose(r1, rows[1].r_km)
    assert np.allclose(v1, rows[1].v_km_s)


def test_interpolate_state_out_of_range_raises():
    rows = mod._parse_oem_rows(SAMPLE_OEM.splitlines())
    before = Time("2026-04-02T03:00:00.000", format="isot", scale="utc")
    with pytest.raises(ValueError, match="outside OEM range"):
        mod.interpolate_state(rows, before)


def test_ra_dec_computation_returns_finite_values():
    # Use a real-ish state vector and time; only test that transform runs and returns finite angles.
    rows = mod._parse_oem_rows(SAMPLE_OEM.splitlines())
    t = rows[0].t
    r_km, _ = mod.interpolate_state(rows, t)

    target_gcrs = mod.SkyCoord(
        x=r_km[0] * mod.u.km,
        y=r_km[1] * mod.u.km,
        z=r_km[2] * mod.u.km,
        frame=mod.GCRS(obstime=t),
        representation_type="cartesian",
    )
    icrs = target_gcrs.transform_to(mod.ICRS())
    ra = icrs.ra.to_value(mod.u.deg)
    dec = icrs.dec.to_value(mod.u.deg)
    assert np.isfinite(ra)
    assert np.isfinite(dec)
    assert 0.0 <= ra < 360.0
    assert -90.0 <= dec <= 90.0

