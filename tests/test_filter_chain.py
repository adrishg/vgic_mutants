"""The four-step Kv2.1 chain: the step order, and the column matching that differs by cohort."""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "filter_stats"))

import filter_chain as fc  # noqa: E402


def _cohort(residue_405: str = "LEU", g377: float = 12.0, interface: float = 20.0) -> pd.DataFrame:
    """Two trajectories of two recycles each, everything well inside every cutoff."""
    rows = []
    for seed in (1, 2):
        for recycle in (0, 1):
            row = {"pdb_file": f"out/model_1_seed_{seed}.r{recycle}.pdb"}
            for a, b in itertools.combinations("ABCD", 2):
                row[f"CA_CA_{a}_GLY377_CA-{b}_GLY377_CA"] = g377
                row[f"CA_CA_{a}_{residue_405}405_CA-{b}_{residue_405}405_CA"] = 16.0
            for chain in "ABCD":
                for first, second in (("LYS429", "ASN181"), ("LYS429", "VAL184"),
                                      ("GLU425", "ASN181"), ("GLU425", "VAL184"),
                                      ("LYS422", "ASN181"), ("LYS422", "VAL184")):
                    row[f"CA_CA_{chain}_{first}_CA-{chain}_{second}_CA"] = interface
            rows.append(row)
    return pd.DataFrame(rows)


def _convergence(representative_recycle: int = 0) -> pd.DataFrame:
    rows = []
    for seed in (1, 2):
        for recycle in (0, 1):
            rows.append({"pdb_file": f"model_1_seed_{seed}.r{recycle}.pdb",
                         "all_ok": str(recycle >= representative_recycle),
                         "earliest_converged_selected": str(recycle == representative_recycle)})
    return pd.DataFrame(rows)


@pytest.mark.parametrize("residue", ["LEU", "ALA"])
def test_ring_columns_match_by_position_not_residue_name(residue):
    """Position 405 is LEU405 in the wild type and ALA405 in L403A; both cohorts must resolve."""
    found = fc.ring_columns(_cohort(residue), 405)
    assert len(found) == 6
    assert all(f"{residue}405" in column for column in found)


def test_ring_columns_refuses_a_position_that_is_not_there():
    with pytest.raises(KeyError):
        fc.ring_columns(_cohort(), 999)


def test_recycle_zero_can_be_the_representative():
    """Over half of L403A's trajectories converge at r0, so a chain that cannot keep r0 is useless."""
    kept, counts = fc.apply_convergence(_cohort(), _convergence(representative_recycle=0))
    assert counts["trajectories_out"] == 2
    assert kept.loc[kept.is_representative, "model_file"].str.endswith(".r0.pdb").all()


def test_g377_ring_rejects_the_whole_trajectory():
    """One recycle past the cutoff takes its trajectory with it, as the published rule has it."""
    frame = _cohort()
    frame.loc[0, "CA_CA_A_GLY377_CA-B_GLY377_CA"] = 25.0          # seed 1, recycle 0 only
    kept, _ = fc.apply_convergence(frame, _convergence())
    kept, counts = fc.apply_g377_ring(kept)
    assert counts["trajectories_rejected"] == 1
    assert set(kept.trajectory) == {"model_1_seed_2.pdb"}


def test_the_interface_qc_sees_the_whole_trajectory():
    """The step order exists so this works: one bad recycle rejects the trajectory, including the
    recycles that look fine. A topology screen run first would have deleted the bad one."""
    frame = _cohort()
    frame.loc[1, "CA_CA_A_LYS429_CA-A_ASN181_CA"] = 31.0          # seed 1, recycle 1 only
    result = fc.run_chain(
        frame, _convergence(),
        topology_verdicts=pd.DataFrame({
            "model_file": [f"model_1_seed_{s}.r0.pdb" for s in (1, 2)],
            "topology_ok": ["True", "True"], "error": ["", ""]}),
        label="synthetic")
    assert result.counts["interface_qc"]["trajectories_rejected"] == 1
    assert list(result.representatives.trajectory) == ["model_1_seed_2.pdb"]


def test_topology_judges_the_representative_only():
    verdicts = pd.DataFrame({"model_file": ["model_1_seed_1.r0.pdb", "model_1_seed_2.r0.pdb"],
                             "topology_ok": ["False", "True"], "error": ["", ""]})
    result = fc.run_chain(_cohort(), _convergence(), topology_verdicts=verdicts, label="synthetic")
    assert result.counts["topology"] == {"representatives_in": 2, "representatives_out": 1,
                                         "could_not_judge": 0}
    assert result.counts["surviving_trajectories"] == 1
    assert set(result.retained.trajectory) == {"model_1_seed_2.pdb"}


def test_a_representative_with_no_verdict_is_an_error_not_a_silent_pass():
    verdicts = pd.DataFrame({"model_file": ["model_1_seed_1.r0.pdb"], "topology_ok": ["True"],
                             "error": [""]})
    with pytest.raises(KeyError):
        fc.run_chain(_cohort(), _convergence(), topology_verdicts=verdicts, label="synthetic")


def test_the_experimental_structures_pass_the_topology_screen():
    """8SD3 and 8SDA both pass: the screen is a misfold detector, not a WT/mutant discriminator."""
    folder = REPO_ROOT / "kv21" / "experimental"
    paths = {p.name: p for p in sorted(folder.glob("8SD*.pdb"))}
    if len(paths) < 2:
        pytest.skip("the experimental structures are not available here")
    verdicts = fc.judge_topology(paths, segments=fc.SEGMENTS_EXPERIMENTAL)
    assert list(verdicts.error) == ["", ""]
    assert verdicts.topology_ok.all()
