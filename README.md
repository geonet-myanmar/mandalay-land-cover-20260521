# Sentinel-2 unsupervised land-cover mapping

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

Reproducible workflow that downloads **Sentinel-2 MSI Level-2A** imagery from the [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/), runs an unsupervised **7-class** spectral classification, and exports a high-resolution, publication-ready land-cover map.

The default case study covers the **Mandalay region, Myanmar** (WGS84 bbox `95.964203, 21.725697, 96.330872, 22.104726`). Later monsoon acquisitions over this AOI are almost fully cloudy, so the published map uses the latest *usable* scene: **21 May 2026** (MGRS tile `46QHK`, ~5% cloud inside the AOI after SCL masking).

<p align="center">
  <img src="outputs/sentinel2_landcover_7class.jpg" alt="7-class unsupervised land-cover map of the Mandalay region from Sentinel-2 L2A on 21 May 2026" width="920">
</p>

<p align="center"><em>300&nbsp;dpi JPEG · 10&nbsp;m pixels · UTM 46N (EPSG:32646)</em></p>

## Features

- STAC search against Planetary Computer (`sentinel-2-l2a`), with a pinned scene for reproducibility or `--latest` to pick the newest low-cloud item
- 10 m stack of B02–B12 (20 m bands resampled on load) plus NDVI, NDWI, NDBI, and visible/NIR brightness
- Scene Classification Layer (SCL) mask for cloud, cloud shadow, cirrus, saturation, and no-data
- Mini-batch k-means with a fixed seed; class names assigned from cluster-mean indices
- GeoTIFF with class colormap and nodata, JPEG map (legend, scale bar, north arrow, true-colour inset), plus `metadata.json` and `legend.csv`

## Results (published scene)

| Value | Class | Valid pixels | Mean NDVI | Interpretation cue |
|------:|-------|-------------:|----------:|--------------------|
| 0 | Water | 2.9% | −0.07 | Highest NDWI, negative NDVI |
| 1 | Dense vegetation | 18.1% | +0.80 | Highest NDVI |
| 2 | Cropland | 24.1% | +0.72 | High NDVI, lower than canopy |
| 3 | Grassland / shrub | 24.4% | +0.54 | Intermediate NDVI |
| 4 | Sparse vegetation | 10.3% | +0.39 | Lower NDVI, mixed soil |
| 5 | Bare soil | 16.0% | +0.28 | Bright, residual vegetation |
| 6 | Built-up | 4.2% | +0.18 | High brightness + NDBI |
| — | Cloud / no data | 5.1% of AOI | — | SCL 3, 8, 9, 10 |

Class names are **heuristic labels**, not a validated land-cover legend. Built-up and bright bare surfaces (for example river sand) are spectrally similar and will mix.

## Quick start

Python 3.10+ with GDAL/rasterio. Conda-forge is the easiest way to get the geospatial stack.

```bash
git clone https://github.com/mrtinkooo/experiment-hub.git
cd experiment-hub

# conda (recommended)
conda env create -f environment.yml
conda activate s2-landcover

# or pip, if GDAL/rasterio are already installed
python -m pip install -r requirements.txt

# reproduce the published Mandalay map
python classify_sentinel2_landcover.py
```

Internet access is required. Planetary Computer asset URLs are signed automatically; no API key is needed for this public collection.

### Other useful commands

```bash
# newest scene with tile cloud cover < 20%
python classify_sentinel2_landcover.py --latest --max-cloud 20

# different AOI
python classify_sentinel2_landcover.py \
  --bbox 96.0 21.8 96.2 22.0 \
  --latest \
  --region "Custom AOI" \
  --out-dir outputs/custom

python classify_sentinel2_landcover.py --help
```

## Outputs

Written to `outputs/` (see [`outputs/README.md`](outputs/README.md)):

| File | Description |
|------|-------------|
| [`sentinel2_landcover_7class.jpg`](outputs/sentinel2_landcover_7class.jpg) | Publication map, 300 dpi, JPEG quality 95 |
| [`sentinel2_landcover_7class.tif`](outputs/sentinel2_landcover_7class.tif) | Class raster, UInt8, LZW, nodata = 255, embedded colormap |
| [`legend.csv`](outputs/legend.csv) | Class values, names, hex colours, and index stats |
| [`metadata.json`](outputs/metadata.json) | Scene id, bbox, CRS, algorithm settings, class table |

## Methods (short)

1. Query the Planetary Computer STAC API for `sentinel-2-l2a` over the AOI.
2. Load B02, B03, B04, B05, B06, B07, B08, B8A, B11, B12, and SCL at 10 m in the native MGRS CRS with [`odc-stac`](https://odc-stac.readthedocs.io/).
3. Convert digital numbers to surface reflectance using the processing-baseline ≥ 04.00 offset/scale: `(DN − 1000) / 10000`.
4. Mask SCL values `{0, 1, 3, 8, 9, 10}`.
5. Append NDVI `(B08−B04)/(B08+B04)`, NDWI `(B03−B08)/(B03+B08)`, NDBI `(B11−B08)/(B11+B08)`, and mean brightness of B02/B03/B04/B08.
6. Standardise features and fit `MiniBatchKMeans` (`k=7`, `random_state=42`) on up to 400 000 random valid pixels; predict every valid pixel.
7. Name clusters from mean indices (water → dense vegetation → built-up → bare soil → remaining vegetation by NDVI).
8. Apply a 3×3 median sieve and render the map.

A longer write-up is in [`docs/methodology.md`](docs/methodology.md).

## Limitations

- Labels are unsupervised interpretations. Do not treat them as ESA WorldCover / Dynamic World equivalents.
- Single date: cropping calendars and phenology are not modelled.
- “Latest available” is not the same as “latest calendar date”. June–August 2026 over this AOI is 80–100% cloud.
- No ground-truth accuracy assessment is included.
- The AOI sits on the UTM 46/47 boundary. The raster stays in the native tile CRS (EPSG:32646) to avoid an extra warp.

## Citation

If you use this repository, please cite the code and the upstream data:

```text
Sentinel-2 unsupervised land-cover mapping (2026).
github.com/geonet-myanmar/mandalay-land-cover-20260521
```

A machine-readable record is in [`CITATION.cff`](CITATION.cff).

**Data.** Contains modified Copernicus Sentinel data (2026) accessed through the Microsoft Planetary Computer STAC API.

- [Copernicus Sentinel data terms](https://sentinels.copernicus.eu/documents/247904/690755/Sentinel_Data_Legal_Notice)
- [Planetary Computer](https://planetarycomputer.microsoft.com/)
- Scene: `S2A_MSIL2A_20260521T041211_R047_T46QHK_20260521T075059`

## License

Code in this repository is released under the [MIT License](LICENSE). Sentinel-2 imagery remains subject to the Copernicus Sentinel data terms. Microsoft Planetary Computer terms apply to access via their API.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md).
