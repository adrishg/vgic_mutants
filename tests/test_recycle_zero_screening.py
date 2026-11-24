"""Recycle 0 must reach the topology screen.

Our convergence step picks each trajectory's earliest converged recycle as its representative, and
for a large share of trajectories that is recycle 0 (over half of Kv2.1 L403A). A topology screen
that starts at recycle 1 therefore does not make the ensemble cleaner - it lets exactly those
representatives through unjudged, and they are the ones the rest of the analysis reports.

These tests pin the two halves of that: the chain keeps an r0 representative through the topology
step, and the dataset configs ask the screen to judge recycle 0.
"""
from __future__ import annotations

import sys
import tomllib
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "filter_stats"))

import filter_chain as fc  # noqa: E402


# =========================================================================== the chain keeps r0

def _cohort() -> pd.DataFrame:
    """Two Kv2.1 trajectories of two recycles each, every measured distance well inside its cutoff."""
    rows = []
    for seed in (1, 2):
        for recycle in (0, 1):
            row = {"pdb_file": f"out/model_1_seed_{seed}.r{recycle}.pdb"}
            for column in _RING_COLUMNS:
                row[column] = 12.0
            for column in _INTERFACE_COLUMNS:
                row[column] = 9.0
            rows.append(row)
    return pd.DataFrame(rows)


_RING_COLUMNS = [f"CA_CA_{a}_GLY377_CA-{b}_GLY377_CA"
                 for a, b in (("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "D"))]

# The 24 columns the published pore-VSD interface QC reads: S6's cytoplasmic end against the
# residues before S1, within each chain.
_INTERFACE_COLUMNS = [f"CA_CA_{chain}_{first}_CA-{chain}_{second}_CA"
                      for chain in "ABCD"
                      for first in ("LYS429", "GLU425", "LYS422")
                      for second in ("ASN181", "VAL184")]


def _convergence() -> pd.DataFrame:
    """Both trajectories converge at recycle 0, so both representatives are r0 models."""
    rows = []
    for seed in (1, 2):
        for recycle in (0, 1):
            rows.append({"pdb_file": f"model_1_seed_{seed}.r{recycle}.pdb",
                         "all_ok": "True", "earliest_converged_selected": str(recycle == 0)})
    return pd.DataFrame(rows)


def test_topology_judges_and_keeps_recycle_zero_representatives():
    """An r0 representative that passes the screen survives the whole chain."""
    verdicts = pd.DataFrame({"model_file": [f"model_1_seed_{s}.r0.pdb" for s in (1, 2)],
                             "topology_ok": ["True", "True"], "error": ["", ""]})
    result = fc.run_chain(_cohort(), _convergence(), channel="Kv2.1", topology_verdicts=verdicts)
    survivors = sorted(result.representatives["model_file"])
    assert survivors == ["model_1_seed_1.r0.pdb", "model_1_seed_2.r0.pdb"]
    assert result.counts["topology"]["representatives_in"] == 2


def test_topology_rejects_a_failing_recycle_zero_representative():
    """And an r0 representative is judged, not waved through: a failing one is removed."""
    verdicts = pd.DataFrame({"model_file": [f"model_1_seed_{s}.r0.pdb" for s in (1, 2)],
                             "topology_ok": ["False", "True"], "error": ["", ""]})
    result = fc.run_chain(_cohort(), _convergence(), channel="Kv2.1", topology_verdicts=verdicts)
    assert sorted(result.representatives["model_file"]) == ["model_1_seed_2.r0.pdb"]


# =========================================================================== the configs ask for it

CONFIGS = sorted((REPO / "datasets").glob("*/*.toml"))


@pytest.mark.parametrize("path", CONFIGS, ids=[p.stem for p in CONFIGS])
def test_dataset_config_screens_recycle_zero(path):
    """first_recycle = 0 in every dataset config, for the reason in this module's docstring."""
    screen = tomllib.loads(path.read_text())["screen"]
    assert screen["first_recycle"] == 0, (
        f"{path.name} would leave recycle 0 unjudged, and our convergence step selects r0 "
        f"representatives")


def test_config_generator_writes_first_recycle_zero():
    """The configs are generated, so the generator is where the value has to be right."""
    source = (REPO / "scripts" / "filtering" / "configs" / "make_configs.py").read_text()
    assert "'first_recycle = 0'" in source
    assert "'first_recycle = 1'" not in source


# =========================================================================== downstream, which assumed r0 absent

def test_region_dv_takes_the_lowest_screened_recycle_from_the_config():
    """Step 3 used to hardcode recycle >= 1 in two places: the check that the verdicts cover step 1's
    models, and the per-group skip that leaves unscreened recycles alone. With first_recycle = 0 the
    first would reject every screen as out of date and the second would silently keep removed models,
    so both have to read the config instead."""
    source = (REPO / "filter_stats" / "pipeline" / "region_dv.py").read_text()
    assert "conf.screen.first_recycle" in source
    assert "dv.recycle >= 1" not in source
    assert "if recycle == 0:" not in source
