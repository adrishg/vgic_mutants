"""Every panel must give one experimental structure the same marker.

A reader identifies a reference by its shape before its color, so a PDB that
is a square in one panel and a circle in the next reads as two structures.
These tests pin the call sites to the canonical table in shared.plotting
rather than to a marker order that happens to line up today.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.markers import MarkerStyle

from shared.plotting import NAV15_EXPERIMENTAL_STYLES, RMSD_REFERENCE_STYLES
from shared.nav15_pore_shape import (
    NAV15_EXPERIMENTAL_GATE_MAPS,
    NAV15_EXPERIMENTAL_STATES,
    plot_shape_landscape,
    plot_shape_distributions,
)


def _canonical_path(marker: str):
    return MarkerStyle(marker).get_path().transformed(
        MarkerStyle(marker).get_transform()
    ).vertices


def _experimental_metrics() -> pd.DataFrame:
    rows = []
    for index, pdb_id in enumerate(NAV15_EXPERIMENTAL_GATE_MAPS):
        rows.append({
            "dataset": pdb_id,
            "state": NAV15_EXPERIMENTAL_STATES[pdb_id],
            "model": pdb_id,
            "mean_diagonal_A": 12.0 + index,
            "side_aspect_ratio": 1.0 + index / 50,
        })
    return pd.DataFrame(rows)


def _ensemble_metrics() -> pd.DataFrame:
    generator = np.random.default_rng(0)
    return pd.DataFrame({
        "dataset": ["WT masked"] * 40,
        "mean_diagonal_A": generator.normal(13.0, 0.4, 40),
        "side_aspect_ratio": generator.normal(1.05, 0.02, 40),
    })


def _markers_by_label(axis) -> dict[str, np.ndarray]:
    found = {}
    for collection in axis.collections:
        label = collection.get_label()
        pdb_id = str(label).split(" | ")[0]
        if pdb_id in NAV15_EXPERIMENTAL_STYLES:
            found[pdb_id] = collection.get_paths()[0].vertices
    return found


def test_pore_shape_landscape_uses_the_canonical_marker_for_every_reference():
    figure, axis = plot_shape_landscape(
        _ensemble_metrics(), _experimental_metrics(), {"WT masked": "#855094"},
    )
    drawn = _markers_by_label(axis)
    plt.close(figure)
    assert set(drawn) == set(NAV15_EXPERIMENTAL_STYLES), (
        "every experimental reference should be drawn and labeled"
    )
    for pdb_id, vertices in drawn.items():
        expected = _canonical_path(NAV15_EXPERIMENTAL_STYLES[pdb_id]["marker"])
        assert vertices.shape == expected.shape and np.allclose(vertices, expected), (
            f"{pdb_id} is drawn with a marker other than "
            f"{NAV15_EXPERIMENTAL_STYLES[pdb_id]['marker']!r}"
        )


def test_pore_shape_panels_agree_with_each_other():
    """The landscape and the distribution panel sit in one figure."""
    ensembles, experimentals = _ensemble_metrics(), _experimental_metrics()
    palette = {"WT masked": "#855094"}

    landscape_figure, landscape_axis = plot_shape_landscape(ensembles, experimentals, palette)
    landscape = _markers_by_label(landscape_axis)
    plt.close(landscape_figure)

    distribution_figure, distribution_axis = plot_shape_distributions(
        ensembles, experimentals, palette,
    )
    handles, _ = distribution_axis.get_legend_handles_labels()
    distribution = {
        str(handle.get_label()).split(" | ")[0]: handle.get_marker()
        for handle in handles
        if str(handle.get_label()).split(" | ")[0] in NAV15_EXPERIMENTAL_STYLES
    }
    plt.close(distribution_figure)

    for pdb_id, marker in distribution.items():
        expected = _canonical_path(marker)
        assert np.allclose(landscape[pdb_id], expected), (
            f"{pdb_id} uses {marker!r} in the distribution panel but a "
            f"different shape in the landscape panel of the same figure"
        )


@pytest.mark.parametrize("table", [RMSD_REFERENCE_STYLES, NAV15_EXPERIMENTAL_STYLES])
def test_every_canonical_marker_is_a_real_matplotlib_marker(table):
    for pdb_id, style in table.items():
        MarkerStyle(style["marker"])
