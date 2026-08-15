#!/usr/bin/env python3
"""Unsupervised Sentinel-2 L2A land-cover classification via Microsoft Planetary Computer.

Downloads Sentinel-2 MSI Level-2A surface reflectance for an AOI, runs mini-batch
k-means on spectral bands plus NDVI/NDWI/NDBI, assigns semantic labels from
cluster-mean indices, and writes a GeoTIFF plus a publication-ready JPEG map.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import odc.stac
import planetary_computer
import pystac_client
import rioxarray  # noqa: F401
import xarray as xr
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle
from matplotlib_scalebar.scalebar import ScaleBar
from scipy.ndimage import median_filter
from shapely.geometry import box, shape
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "sentinel-2-l2a"

# Published case study (Mandalay region, Myanmar). Override with --latest or --item-id.
DEFAULT_BBOX = (95.964203, 21.725697, 96.330872, 22.104726)
DEFAULT_ITEM_ID = "S2A_MSIL2A_20260521T041211_R047_T46QHK_20260521T075059"
DEFAULT_REGION = "Mandalay region, Myanmar"

SPECTRAL_BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
ALL_BANDS = SPECTRAL_BANDS + ["SCL"]
SCL_MASK = {0, 1, 3, 8, 9, 10}

CLASS_SPEC = [
    ("Water", "#1d4e89"),
    ("Dense vegetation", "#14532d"),
    ("Cropland", "#65a30d"),
    ("Grassland / shrub", "#a3b18a"),
    ("Sparse vegetation", "#d4a373"),
    ("Bare soil", "#e9d8a6"),
    ("Built-up", "#9b2226"),
]

# Sentinel-2 processing baseline >= 04.00
REFLECTANCE_OFFSET = 1000.0
REFLECTANCE_SCALE = 10000.0

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_OUT_DIR = REPO_ROOT / "outputs"


@dataclass(frozen=True)
class Config:
    bbox: tuple[float, float, float, float]
    item_id: str | None
    latest: bool
    max_cloud: float
    lookback_days: int
    n_classes: int
    resolution: float
    region: str
    out_dir: Path
    seed: int
    sample_pixels: int


def parse_args(argv: list[str] | None = None) -> Config:
    p = argparse.ArgumentParser(
        description="Classify Sentinel-2 L2A into unsupervised land-cover classes "
        "and export a publication-ready JPEG map.",
    )
    p.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        metavar=("MINX", "MINY", "MAXX", "MAXY"),
        default=list(DEFAULT_BBOX),
        help="AOI as WGS84 lon/lat (default: published Mandalay case study).",
    )
    p.add_argument(
        "--item-id",
        default=DEFAULT_ITEM_ID,
        help="Pin a Planetary Computer Sentinel-2 L2A item id (default: published scene).",
    )
    p.add_argument(
        "--latest",
        action="store_true",
        help="Ignore --item-id and pick the newest scene below --max-cloud.",
    )
    p.add_argument("--max-cloud", type=float, default=25.0, help="Max tile cloud cover %% when using --latest.")
    p.add_argument("--lookback-days", type=int, default=200, help="Search window in days when using --latest.")
    p.add_argument("--n-classes", type=int, default=7, help="Number of k-means classes.")
    p.add_argument("--resolution", type=float, default=10.0, help="Output pixel size in metres.")
    p.add_argument("--region", default=DEFAULT_REGION, help="Place name used in the map title.")
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Directory for JPEG, GeoTIFF, and sidecars.")
    p.add_argument("--seed", type=int, default=42, help="Random seed for sampling and k-means.")
    p.add_argument("--sample-pixels", type=int, default=400_000, help="Max pixels used to fit k-means.")
    args = p.parse_args(argv)
    item_id = None if args.latest else args.item_id
    return Config(
        bbox=tuple(args.bbox),
        item_id=item_id,
        latest=args.latest,
        max_cloud=args.max_cloud,
        lookback_days=args.lookback_days,
        n_classes=args.n_classes,
        resolution=args.resolution,
        region=args.region,
        out_dir=args.out_dir,
        seed=args.seed,
        sample_pixels=args.sample_pixels,
    )


def open_catalog() -> pystac_client.Client:
    return pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)


def item_crs(item) -> str:
    code = item.properties.get("proj:code") or item.properties.get("proj:epsg")
    if isinstance(code, str) and code.upper().startswith("EPSG:"):
        return code
    if isinstance(code, int):
        return f"EPSG:{code}"
    return "EPSG:32646"


def item_tile(item) -> str:
    return str(item.properties.get("s2:mgrs_tile") or item.properties.get("grid:code") or "unknown")


def select_item(catalog: pystac_client.Client, cfg: Config):
    if cfg.item_id:
        item = next(catalog.search(collections=[COLLECTION], ids=[cfg.item_id]).items())
        return item

    end = datetime.now(timezone.utc)
    start = end - timedelta(days=cfg.lookback_days)
    search = catalog.search(
        collections=[COLLECTION],
        bbox=list(cfg.bbox),
        datetime=f"{start.date()}/{end.date()}",
        query={"eo:cloud_cover": {"lt": cfg.max_cloud}},
        sortby=[{"field": "properties.datetime", "direction": "desc"}],
    )
    aoi = box(*cfg.bbox)
    best = None
    best_cov = 0.0
    for item in search.items():
        geom = shape(item.geometry)
        cov = aoi.intersection(geom).area / aoi.area if aoi.area else 0.0
        if cov >= 0.90:
            return item
        if cov > best_cov:
            best, best_cov = item, cov
    if best is None:
        raise RuntimeError(
            f"No Sentinel-2 L2A items found for bbox={cfg.bbox} "
            f"in the last {cfg.lookback_days} days with cloud < {cfg.max_cloud}%."
        )
    return best


def load_scene(item, cfg: Config):
    ds = odc.stac.load(
        [item],
        bands=ALL_BANDS,
        bbox=list(cfg.bbox),
        crs=item_crs(item),
        resolution=cfg.resolution,
        groupby="solar_day",
        fail_on_error=False,
    )
    return ds.squeeze(drop=True).compute()


def to_reflectance(arr: np.ndarray) -> np.ndarray:
    """Convert Sentinel-2 L2A DNs (processing baseline >= 04.00) to reflectance."""
    sr = (arr.astype(np.float32) - REFLECTANCE_OFFSET) / REFLECTANCE_SCALE
    return np.clip(sr, 0.0, 1.0)


def build_features(ds):
    scl = ds.SCL.values.astype(np.uint8)
    valid = ~np.isin(scl, list(SCL_MASK))

    stack = np.stack([to_reflectance(ds[b].values) for b in SPECTRAL_BANDS], axis=0)
    _b02, b03, b04, b08, b11 = (stack[i] for i in (0, 1, 2, 6, 8))

    eps = 1e-6
    ndvi = (b08 - b04) / (b08 + b04 + eps)
    ndwi = (b03 - b08) / (b03 + b08 + eps)
    ndbi = (b11 - b08) / (b11 + b08 + eps)
    brightness = np.mean(stack[[0, 1, 2, 6]], axis=0)

    features = np.concatenate(
        [stack, ndvi[None], ndwi[None], ndbi[None], brightness[None]],
        axis=0,
    )
    return features, valid, stack


def classify(features, valid, cfg: Config):
    n_bands, rows, cols = features.shape
    flat = features.reshape(n_bands, -1).T
    mask = valid.ravel() & np.isfinite(flat).all(axis=1)
    X = flat[mask]
    if X.shape[0] < cfg.n_classes:
        raise RuntimeError("Not enough valid pixels to cluster.")

    scaler = StandardScaler()
    xs = scaler.fit_transform(X)

    rng = np.random.default_rng(cfg.seed)
    sample_n = min(cfg.sample_pixels, xs.shape[0])
    sample_idx = rng.choice(xs.shape[0], size=sample_n, replace=False)

    model = MiniBatchKMeans(
        n_clusters=cfg.n_classes,
        random_state=cfg.seed,
        batch_size=8192,
        n_init=20,
        max_iter=200,
        reassignment_ratio=0.01,
    )
    model.fit(xs[sample_idx])
    pred = model.predict(xs)

    labels = np.full(flat.shape[0], -1, dtype=np.int16)
    labels[mask] = pred
    return labels.reshape(rows, cols), pred, X


def interpret_clusters(pred: np.ndarray, X: np.ndarray, n_classes: int) -> tuple[dict[int, str], dict[int, dict]]:
    """Assign semantic names from cluster-mean spectral indices.

    Feature layout: 10 bands + NDVI, NDWI, NDBI, brightness.
    """
    ndvi = X[:, 10]
    ndwi = X[:, 11]
    ndbi = X[:, 12]
    brightness = X[:, 13]

    stats: dict[int, dict] = {}
    for cluster in range(n_classes):
        m = pred == cluster
        if not np.any(m):
            raise RuntimeError(f"Empty cluster {cluster}; try a different seed or fewer classes.")
        stats[cluster] = {
            "ndvi": float(ndvi[m].mean()),
            "ndwi": float(ndwi[m].mean()),
            "ndbi": float(ndbi[m].mean()),
            "brightness": float(brightness[m].mean()),
            "n": int(m.sum()),
        }

    remaining = set(range(n_classes))
    names: dict[int, str] = {}

    water = max(remaining, key=lambda c: stats[c]["ndwi"] - 0.4 * stats[c]["ndvi"])
    names[water] = "Water"
    remaining.remove(water)

    if n_classes == 1:
        return names, stats

    dense = max(remaining, key=lambda c: stats[c]["ndvi"])
    names[dense] = "Dense vegetation"
    remaining.remove(dense)

    if remaining:
        urban = max(remaining, key=lambda c: stats[c]["ndbi"] + 0.3 * stats[c]["brightness"])
        names[urban] = "Built-up"
        remaining.remove(urban)

    if remaining:
        bare = max(remaining, key=lambda c: stats[c]["brightness"] - stats[c]["ndvi"])
        names[bare] = "Bare soil"
        remaining.remove(bare)

    vegetation_names = ["Cropland", "Grassland / shrub", "Sparse vegetation"]
    extra = [f"Vegetation {i + 4}" for i in range(max(0, len(remaining) - len(vegetation_names)))]
    ordered = sorted(remaining, key=lambda c: stats[c]["ndvi"], reverse=True)
    for cluster, name in zip(ordered, vegetation_names + extra):
        names[cluster] = name

    return names, stats


def remap_classes(class_map: np.ndarray, names: dict[int, str], n_classes: int) -> np.ndarray:
    palette = CLASS_SPEC if n_classes == len(CLASS_SPEC) else [
        (names[i], CLASS_SPEC[i % len(CLASS_SPEC)][1]) for i in range(n_classes)
    ]
    name_to_id = {name: i for i, (name, _) in enumerate(palette)}
    remapped = np.full_like(class_map, -1)
    for src, name in names.items():
        remapped[class_map == src] = name_to_id.get(name, src)
    valid = remapped >= 0
    smoothed = median_filter(np.where(valid, remapped, 0), size=3)
    return np.where(valid, smoothed, -1).astype(np.int16)


def rgb_composite(stack: np.ndarray) -> np.ndarray:
    """True-color RGB from B04, B03, B02 with a 2–98% stretch."""
    rgb = np.stack([stack[2], stack[1], stack[0]], axis=-1)
    out = np.zeros_like(rgb)
    for i in range(3):
        band = rgb[..., i]
        lo, hi = np.nanpercentile(band, (2, 98))
        out[..., i] = np.clip((band - lo) / (hi - lo + 1e-6), 0, 1)
    return out


def class_palette(n_classes: int) -> list[tuple[str, str]]:
    if n_classes == len(CLASS_SPEC):
        return list(CLASS_SPEC)
    return [(CLASS_SPEC[i % len(CLASS_SPEC)][0], CLASS_SPEC[i % len(CLASS_SPEC)][1]) for i in range(n_classes)]


def hex_to_rgba(color: str) -> tuple[int, int, int, int]:
    value = color.lstrip("#")
    r, g, b = (int(value[i : i + 2], 16) for i in (0, 2, 4))
    return r, g, b, 255


def save_geotiff(class_map: np.ndarray, ds, crs: str, path: Path, n_classes: int) -> None:
    encoded = np.where(class_map < 0, 255, class_map).astype(np.uint8)
    da = xr.DataArray(
        encoded,
        dims=("y", "x"),
        coords={"y": ds.y, "x": ds.x},
        name="landcover",
        attrs={
            "long_name": "Unsupervised land-cover class",
            "description": "K-means class ids; 255 = cloud / no data",
        },
    )
    da = da.rio.write_crs(crs).rio.write_nodata(255)
    da.rio.to_raster(path, compress="lzw", dtype="uint8")

    import rasterio

    colormap = {i: hex_to_rgba(color) for i, (_, color) in enumerate(class_palette(n_classes))}
    colormap[255] = (217, 217, 217, 255)
    with rasterio.open(path, "r+") as dst:
        dst.write_colormap(1, colormap)


def add_north_arrow(ax) -> None:
    x0, y0 = 0.935, 0.84
    needle = Polygon(
        [(x0, y0 + 0.085), (x0 - 0.018, y0), (x0, y0 + 0.022), (x0 + 0.018, y0)],
        closed=True,
        facecolor="#1a1a1a",
        edgecolor="#1a1a1a",
        lw=0.4,
        transform=ax.transAxes,
        zorder=10,
        clip_on=False,
    )
    ax.add_patch(needle)
    ax.text(
        x0,
        y0 + 0.102,
        "N",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        zorder=10,
    )


def crs_axis_labels(crs: str) -> tuple[str, str]:
    epsg = crs.upper().replace("EPSG:", "")
    hemi = "N" if epsg.startswith("326") else "S" if epsg.startswith("327") else ""
    zone = epsg[-2:] if len(epsg) == 5 else epsg
    suffix = f"UTM {zone}{hemi}" if hemi else crs
    return f"Easting (m, {suffix})", f"Northing (m, {suffix})"


def render_map(class_map, rgb, ds, item, cfg: Config, jpeg_path: Path, cloud_tile: float) -> None:
    palette = class_palette(cfg.n_classes)
    cmap = ListedColormap([color for _, color in palette])
    display = class_map.astype(float)
    display[class_map < 0] = np.nan

    x = ds.x.values
    y = ds.y.values
    extent = [float(x.min()), float(x.max()), float(y.min()), float(y.max())]
    width_m = float(x.max() - x.min())
    height_m = float(y.max() - y.min())
    crs = item_crs(item)
    xlabel, ylabel = crs_axis_labels(crs)

    date = item.properties["datetime"][:10]
    scene_id = item.id
    tile = item_tile(item)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 0.6,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )

    fig = plt.figure(figsize=(13.5, 12.2), facecolor="white")
    ax = fig.add_axes([0.07, 0.10, 0.68, 0.76])
    ax.set_facecolor("#d9d9d9")
    ax.imshow(
        display,
        cmap=cmap,
        vmin=0,
        vmax=max(cfg.n_classes - 1, 1),
        interpolation="nearest",
        extent=extent,
        origin="upper",
    )
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel, fontsize=9, labelpad=6)
    ax.set_ylabel(ylabel, fontsize=9, labelpad=6)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(length=3.5, width=0.6)
    for spine in ax.spines.values():
        spine.set_color("#333333")
    ax.grid(True, color="white", alpha=0.28, lw=0.4, ls="--")

    ax.add_artist(
        ScaleBar(
            1,
            units="m",
            location="lower left",
            box_alpha=0.88,
            box_color="white",
            color="#1a1a1a",
            font_properties={"size": 8},
            pad=0.3,
            border_pad=0.5,
            scale_loc="bottom",
            length_fraction=0.22,
        )
    )
    add_north_arrow(ax)

    fig.text(0.07, 0.955, f"Unsupervised land-cover classification ({cfg.n_classes} classes)", fontsize=16, fontweight="bold", color="#111111", ha="left", va="top")
    fig.text(0.07, 0.922, f"{cfg.region}  •  Sentinel-2 MSI Level-2A  •  Microsoft Planetary Computer", fontsize=10, color="#333333", ha="left", va="top")
    fig.text(
        0.07,
        0.898,
        f"Acquisition {date}  •  Scene {scene_id}  •  Tile {tile}  •  "
        f"AOI {cfg.bbox[0]:.4f}–{cfg.bbox[2]:.4f}°E, {cfg.bbox[1]:.4f}–{cfg.bbox[3]:.4f}°N",
        fontsize=8.2,
        color="#555555",
        ha="left",
        va="top",
    )

    lax = fig.add_axes([0.775, 0.10, 0.205, 0.76])
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    lax.add_patch(
        FancyBboxPatch(
            (0.02, 0.02),
            0.96,
            0.96,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            facecolor="#f7f7f5",
            edgecolor="#cfcfc8",
            linewidth=0.8,
            transform=lax.transAxes,
        )
    )
    lax.text(0.10, 0.93, "Legend", fontsize=12, fontweight="bold", transform=lax.transAxes)
    lax.text(0.10, 0.895, "K-means clusters, labelled by\nspectral indices", fontsize=7.4, color="#444444", transform=lax.transAxes)

    counts = [(class_map == i).sum() for i in range(cfg.n_classes)]
    valid_n = max(int((class_map >= 0).sum()), 1)
    y0 = 0.82
    step = 0.075 if cfg.n_classes <= 7 else 0.055
    for i, (name, color) in enumerate(palette):
        y = y0 - i * step
        lax.add_patch(Rectangle((0.10, y), 0.10, 0.038, facecolor=color, edgecolor="#222222", linewidth=0.4, transform=lax.transAxes))
        pct = 100.0 * counts[i] / valid_n
        lax.text(0.24, y + 0.018, name, fontsize=8.2, va="center", fontweight="medium", transform=lax.transAxes)
        lax.text(0.24, y - 0.004, f"{pct:.1f}% of valid pixels", fontsize=6.8, color="#666666", va="center", transform=lax.transAxes)

    cloud_pct = 100.0 * (class_map < 0).mean()
    lax.add_patch(Rectangle((0.10, 0.275), 0.10, 0.038, facecolor="#d9d9d9", edgecolor="#888888", linewidth=0.4, transform=lax.transAxes))
    lax.text(0.24, 0.293, "Cloud / no data", fontsize=8.2, va="center", transform=lax.transAxes)
    lax.text(0.24, 0.271, f"{cloud_pct:.1f}% of AOI", fontsize=6.8, color="#666666", va="center", transform=lax.transAxes)

    lax.text(0.10, 0.238, "True-color reference", fontsize=9.5, fontweight="bold", transform=lax.transAxes)
    rgb_ax = lax.inset_axes([0.10, 0.118, 0.80, 0.108])
    rgb_ax.imshow(rgb, origin="upper")
    rgb_ax.set_xticks([])
    rgb_ax.set_yticks([])
    for spine in rgb_ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#444444")

    lax.text(0.10, 0.098, "Methods", fontsize=9.5, fontweight="bold", transform=lax.transAxes)
    methods = (
        f"• {', '.join(SPECTRAL_BANDS)}\n"
        "• NDVI, NDWI, NDBI\n"
        f"• Mini-batch k-means, k = {cfg.n_classes}\n"
        f"• {cfg.resolution:.0f} m pixels; SCL 3, 8–10\n"
        f"• Tile cloud cover {cloud_tile:.1f}%"
    )
    lax.text(0.10, 0.082, methods, fontsize=6.5, va="top", color="#333333", linespacing=1.35, transform=lax.transAxes)

    fig.text(
        0.07,
        0.045,
        "Data: Copernicus Sentinel-2 (ESA) via Microsoft Planetary Computer  •  "
        "Unsupervised spectral classification; class names are index-based interpretations, not a supervised legend.",
        fontsize=7.2,
        color="#555555",
        ha="left",
        va="center",
    )
    fig.text(
        0.07,
        0.025,
        f"Map extent ≈ {width_m / 1000:.1f} × {height_m / 1000:.1f} km  •  CRS {crs}  •  "
        f"Contains modified Copernicus Sentinel data ({date[:4]}).",
        fontsize=7.2,
        color="#555555",
        ha="left",
        va="center",
    )

    jpeg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        jpeg_path,
        dpi=300,
        format="jpeg",
        pil_kwargs={"quality": 95, "optimize": True, "progressive": True},
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.18,
    )
    plt.close(fig)
    print(f"Wrote {jpeg_path}  ({jpeg_path.stat().st_size / 1e6:.1f} MB)")


def write_sidecars(class_map, names, stats, item, cfg: Config, crs: str, jpeg_path: Path, tif_path: Path) -> None:
    palette = class_palette(cfg.n_classes)
    valid_n = max(int((class_map >= 0).sum()), 1)
    classes = []
    legend_lines = ["value,name,color,percent_valid,ndvi,ndwi,ndbi,brightness,pixels"]
    for i, (name, color) in enumerate(palette):
        src = next((c for c, n in names.items() if n == name), None)
        s = stats.get(src, {}) if src is not None else {}
        pct = 100.0 * float((class_map == i).sum()) / valid_n
        classes.append(
            {
                "value": i,
                "name": name,
                "color": color,
                "percent_valid": round(pct, 2),
                "ndvi": s.get("ndvi"),
                "ndwi": s.get("ndwi"),
                "ndbi": s.get("ndbi"),
                "brightness": s.get("brightness"),
                "pixels": int((class_map == i).sum()),
            }
        )
        legend_lines.append(
            f"{i},{name},{color},{pct:.2f},{s.get('ndvi', '')},{s.get('ndwi', '')},"
            f"{s.get('ndbi', '')},{s.get('brightness', '')},{int((class_map == i).sum())}"
        )

    metadata = {
        "title": f"Unsupervised land-cover classification ({cfg.n_classes} classes)",
        "region": cfg.region,
        "item_id": item.id,
        "datetime": item.properties.get("datetime"),
        "mgrs_tile": item_tile(item),
        "bbox_wgs84": list(cfg.bbox),
        "crs": crs,
        "resolution_m": cfg.resolution,
        "n_classes": cfg.n_classes,
        "tile_cloud_cover_percent": item.properties.get("eo:cloud_cover"),
        "nodata": 255,
        "aoi_cloud_or_nodata_percent": round(100.0 * float((class_map < 0).mean()), 2),
        "bands": SPECTRAL_BANDS,
        "indices": ["NDVI", "NDWI", "NDBI", "brightness"],
        "algorithm": "MiniBatchKMeans",
        "random_state": cfg.seed,
        "sample_pixels": cfg.sample_pixels,
        "scl_masked_values": sorted(SCL_MASK),
        "reflectance": {"offset": REFLECTANCE_OFFSET, "scale": REFLECTANCE_SCALE},
        "classes": classes,
        "outputs": {
            "jpeg": jpeg_path.name,
            "geotiff": tif_path.name,
            "legend": "legend.csv",
        },
        "stac_catalog": STAC_URL,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "disclaimer": "Class names are heuristic interpretations of cluster spectra, not a validated land-cover legend.",
    }
    (cfg.out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (cfg.out_dir / "legend.csv").write_text("\n".join(legend_lines) + "\n", encoding="utf-8")
    print(f"Wrote {cfg.out_dir / 'metadata.json'}")
    print(f"Wrote {cfg.out_dir / 'legend.csv'}")


def main(argv: list[str] | None = None) -> int:
    cfg = parse_args(argv)
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    jpeg_path = cfg.out_dir / "sentinel2_landcover_7class.jpg"
    tif_path = cfg.out_dir / "sentinel2_landcover_7class.tif"

    print("Searching Planetary Computer STAC…")
    catalog = open_catalog()
    item = select_item(catalog, cfg)
    crs = item_crs(item)
    print(
        f"Using {item.id}  datetime={item.properties.get('datetime')}  "
        f"cloud={item.properties.get('eo:cloud_cover')}  crs={crs}"
    )

    print("Loading Sentinel-2 L2A spectral bands…")
    ds = load_scene(item, cfg)
    print(ds)

    print("Building feature stack…")
    features, valid, stack = build_features(ds)
    print(f"Valid pixels: {valid.mean() * 100:.2f}%   shape={features.shape}")

    print(f"Running MiniBatchKMeans (k={cfg.n_classes})…")
    class_map, pred, X = classify(features, valid, cfg)
    names, stats = interpret_clusters(pred, X, cfg.n_classes)
    print("Cluster interpretation:")
    for cluster, name in sorted(names.items(), key=lambda kv: kv[1]):
        s = stats[cluster]
        print(
            f"  {cluster}: {name:22s}  n={s['n']:8d}  NDVI={s['ndvi']:+.3f}  "
            f"NDWI={s['ndwi']:+.3f}  NDBI={s['ndbi']:+.3f}  bright={s['brightness']:.3f}"
        )

    remapped = remap_classes(class_map, names, cfg.n_classes)

    print("Saving GeoTIFF…")
    save_geotiff(remapped, ds, crs, tif_path, cfg.n_classes)

    print("Rendering publication map…")
    rgb = rgb_composite(stack)
    rgb[~valid] = 1.0
    cloud_tile = float(item.properties.get("eo:cloud_cover") or float("nan"))
    render_map(remapped, rgb, ds, item, cfg, jpeg_path, cloud_tile)
    write_sidecars(remapped, names, stats, item, cfg, crs, jpeg_path, tif_path)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
