"""The per-prediction verdict table: every prediction accounted for, with the step that took it.

`run_chain` answers what survives. `chain_verdicts` answers what happened to each prediction, which
is what someone outside the project needs: a row per prediction, the measurement each step made, and
one reason when it was removed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "filter_stats"))

import filter_chain as fc  # noqa: E402

RING = [f"CA_CA_{a}_GLY377_CA-{b}_GLY377_CA"
        for a, b in (("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "D"))]
INTERFACE = [f"CA_CA_{chain}_{first}_CA-{chain}_{second}_CA"
             for chain in "ABCD"
             for first in ("LYS429", "GLU425", "LYS422")
             for second in ("ASN181", "VAL184")]
SEEDS = (1, 2, 3, 4)


def _cohort() -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        for recycle in (0, 1):
            row = {"pdb_file": f"out/model_1_seed_{seed}.r{recycle}.pdb"}
            row.update({c: 12.0 for c in RING})
            row.update({c: 9.0 for c in INTERFACE})
            rows.append(row)
    return pd.DataFrame(rows)


def _convergence(fails=()) -> pd.DataFrame:
    """Every trajectory converges at recycle 0 except the ones named in `fails`."""
    rows = []
    for seed in SEEDS:
        for recycle in (0, 1):
            ok = seed not in fails
            rows.append({"pdb_file": f"model_1_seed_{seed}.r{recycle}.pdb",
                         "all_ok": str(ok),
                         "earliest_converged_selected": str(ok and recycle == 0)})
    return pd.DataFrame(rows)


def _topology(fails=()) -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        ok = seed not in fails
        row = {"model_file": f"model_1_seed_{seed}.r0.pdb", "topology_ok": str(ok), "error": ""}
        row.update({f"{c}_ok": str(ok if c == "wiring" else True) for c in fc.TOPOLOGY_CHECKS})
        rows.append(row)
    return pd.DataFrame(rows)


def _measures(detached=()) -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        for recycle in (0, 1):
            rows.append({"pdb_file": f"model_1_seed_{seed}.r{recycle}.pdb",
                         "s6_linker_max_A": 30.0 if seed in detached and recycle == 1 else 6.0})
    return pd.DataFrame(rows)


def _run(**kwargs) -> pd.DataFrame:
    return fc.chain_verdicts(_cohort(), kwargs.pop("convergence", _convergence()),
                             channel="Kv2.1", topology_verdicts=kwargs.pop("topology", _topology()),
                             interface_measures=kwargs.pop("measures", _measures()),
                             dataset="demo", condition="WT", **kwargs)


def test_every_prediction_gets_exactly_one_row():
    table = _run()
    assert len(table) == len(SEEDS)
    assert table["prediction"].is_unique
    assert set(table["dataset"]) == {"demo"} and set(table["condition"]) == {"WT"}


def test_a_clean_cohort_passes_whole():
    table = _run()
    assert set(table["status"]) == {"passed"}
    assert set(table["removed_by"]) == {""}
    assert set(table["removed_because"]) == {""}


def test_each_step_claims_its_own_losses():
    table = _run(convergence=_convergence(fails=(1,)), measures=_measures(detached=(2,)),
                 topology=_topology(fails=(3,)))
    by_seed = table.set_index("seed")["removed_by"].to_dict()
    assert by_seed[1] == "convergence"
    assert by_seed[2] == "pore_vsd_interface_s6_linker"
    assert by_seed[3] == "topology"
    assert by_seed[4] == ""


def test_a_removed_prediction_carries_a_sentence_a_reader_can_act_on():
    table = _run(topology=_topology(fails=(3,)))
    row = table.set_index("seed").loc[3]
    assert row["removed_because"] == fc.STEP_REASON["topology"]
    assert row["topology_result"] == "fail"
    assert row["topology_wiring_result"] == "fail"
    assert row["topology_box_corners_result"] == "pass"


def test_measurements_are_reported_for_every_step_that_ran():
    table = _run(measures=_measures(detached=(2,)))
    row = table.set_index("seed").loc[2]
    assert row["selectivity_filter_ring_max_A"] == pytest.approx(12.0)
    assert row["pore_vsd_s6_linker_max_A"] == pytest.approx(30.0)
    assert row["pore_vsd_published_max_A"] == pytest.approx(9.0)
    assert row["pore_vsd_s6_linker_result"] == "fail"
    assert row["pore_vsd_published_result"] == "pass"


def test_a_check_that_does_not_exist_for_the_channel_says_so_rather_than_passing():
    """Nav1.5 has no published pore-VSD cutoff. A blank column would read as 'it passed'."""
    ring = [f"CA_{a}_CA-{b}_CA" for a, b in (
        ("ASP373", "GLU704"), ("ASP373", "LYS1103"), ("ASP373", "ALA1395"),
        ("GLU704", "LYS1103"), ("GLU704", "ALA1395"), ("LYS1103", "ALA1395"))]
    rows = []
    for seed in SEEDS:
        for recycle in (0, 1):
            row = {"pdb_file": f"out/model_1_seed_{seed}.r{recycle}.pdb"}
            row.update({c: 10.0 for c in ring})
            rows.append(row)
    table = fc.chain_verdicts(pd.DataFrame(rows), _convergence(), channel="Nav1.5",
                              topology_verdicts=_topology(), interface_measures=_measures(),
                              dataset="nav15_vanilla", condition="WT")
    assert set(table["pore_vsd_published_result"]) == {"not_applicable"}
    assert set(table["status"]) == {"passed"}


def test_a_representative_with_no_verdict_stops_the_run():
    """Silently treating an unjudged representative as a failure would misattribute the loss."""
    short = _topology().iloc[:2]
    with pytest.raises(KeyError, match="no verdict"):
        _run(topology=short)
