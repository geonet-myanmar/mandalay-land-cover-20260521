# Methodology

This document describes the default Mandalay case study. CLI flags can change the AOI, scene, class count, and resolution; the processing steps stay the same.

## Area of interest

| Item | Value |
|------|--------|
| WGS84 bbox | `95.964203, 21.725697, 96.330872, 22.104726` (min lon, min lat, max lon, max lat) |
| Place | Mandalay region, Myanmar (Ayeyarwady / Irrawaddy floodplain and eastern foothills) |
| Extent | ≈ 38.8 × 42.8 km |
| Native raster CRS | EPSG:32646 (UTM zone 46N), taken from MGRS tile `46QHK` |

The AOI sits on the UTM 46/47 meridian (96°E). Tile `46QHK` covers 100% of the bbox. Neighbouring tiles `47QKE` (~64%) and `46QGK` (~9%) are not required. The raster is kept in EPSG:32646 so 10 m Sentinel-2 pixels are not warped across the zone boundary.

## Scene selection

Sentinel-2 L2A items are discovered through the [Microsoft Planetary Computer STAC API](https://planetarycomputer.microsoft.com/api/stac/v1) (`collection=sentinel-2-l2a`).

“Latest available” is defined as the **most recent acquisition that is usable for land-cover mapping**, not the newest file on the archive. Over this AOI, June–August 2026 tiles are 78–100% cloud (southwest monsoon). Tile-level cloud on `46QHK` and SCL-based cloud inside the bbox were both checked:

| Date | `46QHK` tile cloud | SCL cloud/shadow in AOI | Notes |
|------|-------------------:|------------------------:|-------|
| 2026-08-12 | 99.5% | — | Unusable |
| 2026-06-10 | 46.8% | — | Still too cloudy |
| **2026-05-21** | **21.1%** | **5.1%** | **Published scene** |
| 2026-04-24 | 11.2% | 1.9% | Cleaner, older |
| 2026-04-14 | 2.8% | 0.7% | Cleanest spring scene |

Published item:

```text
S2A_MSIL2A_20260521T041211_R047_T46QHK_20260521T075059
2026-05-21T04:12:11Z
```

`python classify_sentinel2_landcover.py` pins this id. `python classify_sentinel2_landcover.py --latest` searches the last `--lookback-days` (default 200) and returns the newest item with tile cloud cover below `--max-cloud` (default 25%) that covers ≥ 90% of the AOI.

## Preprocessing

1. **Load.** `odc.stac.load` clips to the bbox, reprojects to the item CRS, and aligns every band to `--resolution` metres (default 10). Twenty-metre bands (B05, B06, B07, B8A, B11, B12, SCL) are resampled onto that grid.
2. **Reflectance.** Processing baseline ≥ 04.00 stores L2A as `DN = 10000 × ρ + 1000`. Conversion is `(DN − 1000) / 10000`, then clipped to `[0, 1]`.
3. **Cloud mask.** Pixels whose SCL value is in `{0, 1, 3, 8, 9, 10}` (no data, saturated, cloud shadow, medium cloud, high cloud, thin cirrus) are excluded from clustering and drawn as “Cloud / no data”. Vegetation (4), not-vegetated (5), water (6), and unclassified (7) remain.

## Features

| Feature | Formula / source |
|---------|------------------|
| Spectral | B02, B03, B04, B05, B06, B07, B08, B8A, B11, B12 |
| NDVI | `(B08 − B04) / (B08 + B04)` |
| NDWI (McFeeters) | `(B03 − B08) / (B03 + B08)` |
| NDBI | `(B11 − B08) / (B11 + B08)` |
| Brightness | mean of B02, B03, B04, B08 |

Each feature is standardised (zero mean, unit variance) on valid pixels before clustering.

## Clustering and labelling

`sklearn.cluster.MiniBatchKMeans` is fit on a reproducible random sample of at most 400 000 valid pixels (`random_state=42`, `n_init=20`, `batch_size=8192`) and then applied to every valid pixel.

Semantic names are **not** learned from training data. They are assigned from cluster-mean indices, in order:

1. **Water** — maximise `NDWI − 0.4 × NDVI`
2. **Dense vegetation** — maximise NDVI among the rest
3. **Built-up** — maximise `NDBI + 0.3 × brightness`
4. **Bare soil** — maximise `brightness − NDVI`
5. Remaining clusters, high to low NDVI — **Cropland**, **Grassland / shrub**, **Sparse vegetation**

A 3×3 median filter is applied only on valid pixels to reduce salt-and-pepper noise. Class ids in the GeoTIFF match the legend (0 = Water … 6 = Built-up). Nodata is 255 (uint8) so GIS software can display the embedded colour table.

This labelling is a convenience for cartography. Spectrally similar surfaces (bright sand vs. rooftops, vigorous crops vs. woodland) can be swapped or mixed. Treat the names as an index-based interpretation.

## Cartography

The JPEG is 13.5 × 12.2 in at 300 dpi (quality 95). It includes a classified map in the native UTM grid, a 5 km scale bar, a north arrow, a class legend with area shares, a true-colour (B04-B03-B02, 2–98% stretch) thumbnail, methods notes, and Copernicus / Planetary Computer attribution.

## Why not a supervised product?

The brief asked for an unsupervised classification from Sentinel-2 bands. No in-situ labels, cadastral layer, or global land-cover prior is used, so no overall accuracy or kappa is reported. For an accuracy-assessed map, sample the GeoTIFF against independent reference data or compare it to ESA WorldCover / Dynamic World as a separate study.

## Software

Developed with Python 3.12, `pystac-client`, `planetary-computer`, `odc-stac`, `rioxarray`/`rasterio`, `scikit-learn`, `scipy`, and `matplotlib`. Exact install targets are in `environment.yml` and `requirements.txt`.
