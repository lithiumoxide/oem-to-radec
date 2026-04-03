import oem_to_radec as mod


def test_visibility_below_horizon():
    above, label = mod._visibility_from_alt_deg(-5.0)
    assert above is False
    assert label == ""


def test_visibility_close_to_horizon_range():
    # Just above 0 and below 20 should be "close_to_horizon"
    above, label = mod._visibility_from_alt_deg(0.1)
    assert above is True
    assert label == "close_to_horizon"

    above2, label2 = mod._visibility_from_alt_deg(19.999)
    assert above2 is True
    assert label2 == "close_to_horizon"


def test_visibility_visible_at_or_above_20():
    above, label = mod._visibility_from_alt_deg(20.0)
    assert above is True
    assert label == "visible"

    above2, label2 = mod._visibility_from_alt_deg(45.0)
    assert above2 is True
    assert label2 == "visible"
