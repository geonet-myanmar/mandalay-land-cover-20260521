# Output products

These files are the published Mandalay case study (Sentinel-2 L2A, 21 May 2026). Re-running `python classify_sentinel2_landcover.py` from the repository root overwrites them.

| File | Role |
|------|------|
| `sentinel2_landcover_7class.jpg` | Publication-ready map (300 dpi, JPEG quality 95) |
| `sentinel2_landcover_7class.tif` | Class raster, UInt8, LZW compressed, CRS EPSG:32646, nodata = 255, GDAL colour table |
| `legend.csv` | Class value, name, hex colour, area share, and mean indices |
| `metadata.json` | Full processing record (scene id, bbox, algorithm, class table) |
| `sentinel2_landcover_7class.qml` | Optional QGIS style for the GeoTIFF |

## GeoTIFF class values

| Value | Name | Colour |
|------:|------|--------|
| 0 | Water | `#1d4e89` |
| 1 | Dense vegetation | `#14532d` |
| 2 | Cropland | `#65a30d` |
| 3 | Grassland / shrub | `#a3b18a` |
| 4 | Sparse vegetation | `#d4a373` |
| 5 | Bare soil | `#e9d8a6` |
| 6 | Built-up | `#9b2226` |
| 255 | Cloud / no data | `#d9d9d9` |

QGIS and GDAL should pick up the embedded colour table. If they do not, load `legend.csv`.

## Provenance

- ESA Copernicus Sentinel-2 MSI Level-2A
- Distributed by the Microsoft Planetary Computer (`https://planetarycomputer.microsoft.com/api/stac/v1`)
- Item `S2A_MSIL2A_20260521T041211_R047_T46QHK_20260521T075059`
- Contains modified Copernicus Sentinel data (2026)
