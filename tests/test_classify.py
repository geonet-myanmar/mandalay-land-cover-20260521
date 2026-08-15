"""Unit tests that do not download imagery."""

from __future__ import annotations

import numpy as np

from classify_sentinel2_landcover import (
    CLASS_SPEC,
    hex_to_rgba,
    interpret_clusters,
    parse_args,
    to_reflectance,
)


def test_to_reflectance_baseline_04():
    dn = np.array([1000.0, 2000.0, 11000.0, 0.0], dtype=np.float32)
    sr = to_reflectance(dn)
    np.testing.assert_allclose(sr, [0.0, 0.1, 1.0, 0.0], atol=1e-6)


def test_hex_to_rgba():
    assert hex_to_rgba("#1d4e89") == (29, 78, 137, 255)


def test_interpret_clusters_seven_classes():
    # 7 well-separated synthetic clusters in the index slots (cols 10–13).
    rows = []
    means = [
        (-0.10, 0.40, 0.00, 0.08),   # water
        (0.85, -0.70, -0.30, 0.12),  # dense veg
        (0.70, -0.60, -0.20, 0.11),  # cropland
        (0.50, -0.45, -0.05, 0.13),  # grass
        (0.35, -0.30, 0.00, 0.10),   # sparse
        (0.15, -0.20, 0.05, 0.16),   # bare
        (0.10, -0.15, 0.20, 0.25),   # built-up
    ]
    pred_parts = []
    for i, (ndvi, ndwi, ndbi, bright) in enumerate(means):
        block = np.zeros((50, 14), dtype=np.float32)
        block[:, 10] = ndvi
        block[:, 11] = ndwi
        block[:, 12] = ndbi
        block[:, 13] = bright
        rows.append(block)
        pred_parts.append(np.full(50, i, dtype=np.int32))
    X = np.vstack(rows)
    pred = np.concatenate(pred_parts)

    names, stats = interpret_clusters(pred, X, n_classes=7)
    assert set(names.values()) == {name for name, _ in CLASS_SPEC}
    assert names[0] == "Water"
    assert names[1] == "Dense vegetation"
    assert names[6] == "Built-up"
    assert stats[1]["ndvi"] > stats[4]["ndvi"]


def test_parse_args_defaults():
    cfg = parse_args([])
    assert cfg.n_classes == 7
    assert cfg.resolution == 10.0
    assert cfg.latest is False
    assert cfg.item_id is not None


def test_parse_args_latest_clears_item():
    cfg = parse_args(["--latest", "--max-cloud", "15"])
    assert cfg.latest is True
    assert cfg.item_id is None
    assert cfg.max_cloud == 15.0
