# Contributing

Thanks for helping improve this workflow.

## Development setup

```bash
conda env create -f environment.yml
conda activate s2-landcover
python -m pip install pytest
```

Or install `requirements.txt` into an environment that already has GDAL/rasterio.

## What to change

- **Scene / AOI behaviour** belongs in `classify_sentinel2_landcover.py` (CLI flags, not hardcoded paths).
- **Method details** belong in `docs/methodology.md` and a short update to the README methods section.
- **Published products** (`outputs/*.jpg`, `outputs/*.tif`) should only be regenerated when the default case study changes. Commit the sidecars (`metadata.json`, `legend.csv`) with them.

## Tests

Unit tests cover reflectance conversion and the cluster-labelling rules. They do not download imagery.

```bash
pytest -q
```

A full end-to-end run hits the Planetary Computer STAC API and pulls several hundred MB of COGs:

```bash
python classify_sentinel2_landcover.py --help
```

Do not add that download to CI.

## Style

- Python 3.10+, four-space indent, no unused imports.
- Keep the default `--item-id` pinned so the published Mandalay map stays reproducible.
- Do not commit `__pycache__/`, virtualenvs, or `copilot-session-*.md`.

## Pull requests

1. Describe the scientific or engineering change in a few sentences.
2. If class colours or labelling rules change, regenerate the map and update the README table.
3. Do not add API keys or signed Planetary Computer URLs. Signing happens at runtime.
