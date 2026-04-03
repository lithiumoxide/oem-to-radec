## Ephemeris to RA/Dec (and Alt/Az) for Artemis II OEM

This repo contains a small script to convert a **CCSDS OEM (Orbit Ephemeris Message)** file with Earth-centered Cartesian state vectors into telescope-friendly **RA/Dec** (ICRS/J2000-like), and optionally **Alt/Az** for a given observer location.

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

#### RA/Dec only (no location required)

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675"
```

#### Time series (RA/Dec only)

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675" ^
  --interval 10 --positions 20
```

#### RA/Dec + Alt/Az (provide location)

Longitude is positive-east (so Ireland is negative), latitude is positive-north.

```bash
python "c:\repos\oem_to_radec\oem_to_radec.py" ^
  --oem "c:\repos\oem_to_radec\Artemis_II_OEM_2026_04_02_to_EI_v3.asc" ^
  --time "2026-04-02T07:10:56.675" ^
  --lat 53.8616 --lon -6.5378 --height-m 30 ^
  --interval 10 --positions 20
```

### CSV output

Use `--output-to-csv <output_dir>`.

- The generated filename is: `<OBJECT_NAME>_<YYYYMMDDTHHMMSSZ>.csv`
- Columns always include: `time_utc, ra_icrs, dec_icrs`
- If a location is provided (and Alt/Az is enabled), the CSV also includes: `alt_deg, az_deg`

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

