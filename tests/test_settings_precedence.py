import argparse

import oem_to_radec as mod


def make_args(lat=None, lon=None, height_m=None, visibility_deg=None):
    ns = argparse.Namespace()
    ns.lat = lat
    ns.lon = lon
    ns.height_m = height_m
    ns.visibility_deg = visibility_deg
    return ns


def test_settings_used_when_cli_missing():
    args = make_args(lat=None, lon=None, height_m=None, visibility_deg=None)
    settings = {"lat": 10.0, "lon": 20.0, "height_m": 100.0, "visibility_deg": 25.0}
    lat, lon, h, vis = mod._resolve_effective_inputs(args, settings)
    assert lat == 10.0 and lon == 20.0
    assert h == 100.0
    assert vis == 25.0


def test_cli_overrides_settings():
    args = make_args(lat=1.0, lon=2.0, height_m=5.0, visibility_deg=35.0)
    settings = {"lat": 10.0, "lon": 20.0, "height_m": 100.0, "visibility_deg": 25.0}
    # visibility_deg override is applied in main before calling resolver; emulate that:
    settings["visibility_deg"] = args.visibility_deg
    lat, lon, h, vis = mod._resolve_effective_inputs(args, settings)
    assert lat == 1.0 and lon == 2.0
    assert h == 5.0
    assert vis == 35.0


def test_defaults_when_neither_provided():
    args = make_args(lat=None, lon=None, height_m=None, visibility_deg=None)
    settings = {}
    lat, lon, h, vis = mod._resolve_effective_inputs(args, settings)
    assert lat == 0.0 and lon == 0.0
    assert h == 0.0
    assert vis == 20.0


def test_visibility_function_uses_threshold():
    # At 15 deg: with threshold 20 -> close_to_horizon; with threshold 10 -> visible
    above_a, label_a = mod._visibility_from_alt_deg(15.0, threshold_visible_deg=20.0)
    above_b, label_b = mod._visibility_from_alt_deg(15.0, threshold_visible_deg=10.0)
    assert (above_a, label_a) == (True, "close_to_horizon")
    assert (above_b, label_b) == (True, "visible")
