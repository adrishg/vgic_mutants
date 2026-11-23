"""The Nav1.5 IFM and Cav1.2 R406 contact columns come from one atom definition."""

import pandas as pd
import pytest

from shared import side_chain_supplements as supplements


@pytest.mark.parametrize(
    "name, expected",
    [
        ("26-07-25_Cav1.2_g406r_vanillaAF2_distances_all_ok_rmsd_3A.csv",
         "cav12/dataDistances/side_chain_20261123/Cav12_g406r_vanilla_R406_side_chain.csv"),
        ("26-11-10_Cav12_g406r_masked_v3_50pct_AF2_distances_all_ok_rmsd_3A.csv",
         "cav12/dataDistances/side_chain_20261123/Cav12_g406r_masked_v3_50pct_R406_side_chain.csv"),
        ("26-07-27_Nav15_wt_vanillaAF2_distances_all_ok_rmsd_3A.csv",
         "nav15/dataDistances/side_chain_20261123/Nav15_wt_vanilla_IFM_side_chain.csv"),
        ("26-11-10_Nav15_qqq_masked_v2_50pct_AF2_distances_all_ok_rmsd_3A.csv",
         "nav15/dataDistances/side_chain_20261123/Nav15_qqq_masked_v2_50pct_IFM_side_chain.csv"),
        ("26-02-10_Cav12_g402s_vanillaAF2_distances_all_ok_rmsd_3A.csv", None),
        ("Kv2.1_wt_vanilla_distances.csv", None),
    ],
)
def test_supplement_path(name, expected):
    assert supplements.supplement_path(name) == expected


def _patch(monkeypatch, tmp_path, table):
    path = tmp_path / "supplement.csv"
    table.to_csv(path, index=False)
    monkeypatch.setattr(supplements.data_access, "resolve", lambda relative: path)


def test_apply_replaces_and_adds_columns(monkeypatch, tmp_path):
    _patch(monkeypatch, tmp_path, pd.DataFrame({
        "pdb_file": ["a.r1.pdb", "b.r1.pdb", "c.r1.pdb"],
        "shortest_ARG406-ASP1528": [3.1, 4.2, 5.3],
        "shortest_ARG406-GLU411": [6.0, 7.0, 8.0],
    }))
    frame = pd.DataFrame({
        "pdb_file": ["/models/b.r1.pdb", "/models/a.r1.pdb"],
        "shortest_ARG406-ASP1528": [9.9, 9.9],
        "CA_x": [1.0, 2.0],
    }, index=[7, 3])
    name = "26-07-25_Cav1.2_g406r_vanillaAF2_distances_all_ok_rmsd_3A.csv"
    out = supplements.apply_side_chain_supplement(frame, name)
    assert out["shortest_ARG406-ASP1528"].tolist() == [4.2, 3.1]
    assert out["shortest_ARG406-GLU411"].tolist() == [7.0, 6.0]
    assert out["CA_x"].tolist() == [1.0, 2.0] and list(out.index) == [7, 3]
    kept = supplements.apply_side_chain_supplement(frame, name, add_missing=False)
    assert "shortest_ARG406-GLU411" not in kept.columns


def test_apply_requires_every_retained_model(monkeypatch, tmp_path):
    _patch(monkeypatch, tmp_path, pd.DataFrame({
        "pdb_file": ["a.r1.pdb"], "shortest_ARG406-ASP1528": [3.1],
    }))
    frame = pd.DataFrame({"pdb_file": ["a.r1.pdb", "zz.r1.pdb"]})
    with pytest.raises(KeyError):
        supplements.apply_side_chain_supplement(
            frame, "26-07-25_Cav1.2_g406r_vanillaAF2_distances_all_ok_rmsd_3A.csv"
        )


def test_tables_without_a_supplement_pass_through():
    frame = pd.DataFrame({"pdb_file": ["a.pdb"], "x": [1.0]})
    assert supplements.apply_side_chain_supplement(frame, "Kv2.1_wt_vanilla.csv") is frame


def test_side_chain_atoms_exclude_backbone_and_hydrogen(tmp_path):
    def atom(serial, name, resname, resid, x, element):
        return (f"ATOM  {serial:5d} {name:<4s} {resname} A{resid:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n")
    pdb = tmp_path / "toy.pdb"
    pdb.write_text(
        atom(1, " CA", "ARG", 406, 0.0, "C") + atom(2, " NH1", "ARG", 406, 1.0, "N")
        + atom(3, " HH1", "ARG", 406, 1.5, "H") + atom(4, " O", "ASP", 1528, 2.0, "O")
        + atom(5, " OD1", "ASP", 1528, 6.0, "O") + atom(6, " CB", "ASP", 1528, 8.0, "C")
    )
    residues = supplements.side_chain_atoms(pdb)
    assert len(residues[("ARG", 406)]) == 1 and len(residues[("ASP", 1528)]) == 2
    distance = supplements.shortest_side_chain_distance(residues, ("ARG", 406), ("ASP", 1528))
    assert distance == pytest.approx(5.0)
