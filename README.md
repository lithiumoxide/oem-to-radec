## Ephemeris to RA/Dec (and Alt/Az) for CCSDS Orbit Ephemeris Message files

This repo contains a small script to convert a **CCSDS OEM (Orbit Ephemeris Message)** file with Earth-centered Cartesian state vectors into telescope-friendly **RA/Dec** (ICRS/J2000-like), and optionally **Alt/Az** for a given observer location.

This may be useful for converting OEM files issued by NASA and other space agencies into RA/Dec coordinates, such as the Artemis II ephemeris data [available here](https://www.nasa.gov/missions/artemis/artemis-2/track-nasas-artemis-ii-mission-in-real-time/).

### Files
- `oem_to_radec.py`: main CLI
- `requirements.txt`: runtime deps
- `requirements-dev.txt`: dev deps (tests)
- `tests/`: unit tests

### Install

```bash
python -m pip install -r "c:\repos\oem_to_radec\requirements.txt"
```

### Usage

#### RA/Dec only

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675"
```

Notes:
- By default, the observer location is `(lat, lon, height_m) = (0, 0, 0)` and Alt/Az will be shown.
- Use `--no-altaz` to suppress Alt/Az output.
- You can set persistent defaults in `settings.ini` (see below).

#### Time series (RA/Dec only)

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675" ^
  --interval 10 --positions 20
```

#### RA/Dec + Alt/Az (provide location or use defaults)

Longitude is positive-east, latitude is positive-north. If you do not provide a location on the CLI or in `settings.ini`, the default `(0, 0, 0 m)` is used and Alt/Az will still be computed.

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675" ^
  --lat 53.8616 --lon -6.5378 --height-m 30 ^
  --interval 10 --positions 20
```

When a location is provided, the output also indicates whether the target is above the local horizon and a visibility label:

- **Above horizon**: printed as "Yes"/"No"
- **Visibility labels** (based on altitude):
  - `< 20°`: "close to horizon"
  - `≥ 20°`: "visible"

### Settings via settings.ini

You can provide default inputs in a `settings.ini` file. The tool looks for it in the current working directory first, then next to `oem_to_radec.py`. CLI flags always override the file.

Section and keys:

```ini
[defaults]
# Observer location (optional). If lat/lon are omitted here and not provided on CLI, Alt/Az is suppressed.
lat = 33.8812        # degrees (+N)
lon = -9.5301        # degrees (+E)
height_m = 37        # meters above ellipsoid

# Visibility threshold in degrees for the "visible" label (default 20 if not set)
visibility_deg = 20
```

If `lat`, `lon`, or `height_m` are omitted in both CLI and `settings.ini`, they default to `0, 0, 0`.

CLI overrides example:

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675" ^
  --visibility-deg 25
```

### CSV output

Use `--output-to-csv <output_dir>`.

- The generated filename is: `<OBJECT_NAME>_<YYYYMMDDTHHMMSSZ>.csv`
- Columns always include: `time_utc, ra_icrs, dec_icrs`
- If Alt/Az is enabled (default unless `--no-altaz`), the CSV also includes: `alt_deg, az_deg, above_horizon, visibility`
- `above_horizon` is `true` or `false`. `visibility` is included only when it matches a label ("close_to_horizon" or "visible").

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675" ^
  --lat 53.8616 --lon -6.5378 --height-m 30 ^
  --interval 10 --positions 20 ^
  --output-to-csv "c:\repos\oem_to_radec\csv_out"
```

### Run tests

```bash
python -m pip install -r "c:\repos\oem_to_radec\requirements-dev.txt"
python -m pytest -q "c:\repos\oem_to_radec\tests"
```

