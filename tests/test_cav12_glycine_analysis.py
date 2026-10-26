import matplotlib
matplotlib.use("Agg")

import pandas as pd
import pytest

from shared.cav12_glycine_analysis import (
    matched_ca_table,
    plot_g402s_paired_revision,
    plot_g406r_screen_and_contacts,
)
from shared import data_access
from shared.dataset_selection import load_selected_distance_csv


def test_backbone_comparison_ignores_added_sidechain_distance():
    wt = pd.DataFrame({"pdb_file": ["wt.pdb"], "CA_GLY402_CA-ILE1523_CA": [5.0],
                       "shortest_GLY402-ILE1523": [4.0]})
    mutant = pd.DataFrame({"pdb_file": ["mutant.pdb"], "CA_SER402_CA-ILE1523_CA": [5.0],
                           "shortest_SER402-ILE1523": [2.0]})
    table = matched_ca_table({"vanilla": (wt, mutant)}, 402, "SER402", ["ILE1523"])
    assert table["Cα distance (Å)"].tolist() == [5.0, 5.0]
    assert table["Source column"].tolist() == [
        "CA_GLY402_CA-ILE1523_CA", "CA_SER402_CA-ILE1523_CA"]
    assert table.pdb_file.tolist() == ["wt.pdb", "mutant.pdb"]


def test_backbone_comparison_fails_instead_of_using_legacy_or_dropping_rows():
    wt = pd.DataFrame({"pdb_file": ["wt.pdb"], "CA_GLY406_CA-ASP1528_CA": [5.0]})
    mutant = pd.DataFrame({"pdb_file": ["mutant.pdb"], "shortest_ARG406-ASP1528": [2.0]})
    with pytest.raises(KeyError, match="atom-matched"):
        matched_ca_table({"vanilla": (wt, mutant)}, 406, "ARG406", ["ASP1528"])
    mutant["CA_ARG406_CA-ASP1528_CA"] = float("nan")
    with pytest.raises(ValueError, match="Nonfinite"):
        matched_ca_table({"vanilla": (wt, mutant)}, 406, "ARG406", ["ASP1528"])


def test_cached_qc_table_wins_over_unfiltered_fallback(tmp_path, monkeypatch):
    root, cache = tmp_path / "repo", tmp_path / "cache"
    root.mkdir()
    cache.mkdir()
    (root / "all.csv").write_text("pdb_file,value\nraw.pdb,99\n")
    (cache / "qc.csv").write_text("pdb_file,value\nqc.pdb,5\n")
    monkeypatch.setattr(data_access, "REPO_ROOT", root)
    monkeypatch.setattr(data_access, "CACHE_ROOT", cache)
    monkeypatch.setattr(data_access, "EXTRA_ROOTS", [])
    monkeypatch.setenv("VGIC_OFFLINE", "1")
    frame = load_selected_distance_csv("test", {"all": root / "all.csv",
                                               "all_ok_3": root / "qc.csv"},
                                       "all_ok_3", fallback_to_all=False)
    assert frame.pdb_file.tolist() == ["qc.pdb"]
    assert frame.attrs["actual_selection"] == "all_ok_3"
    with pytest.raises(data_access.DataUnavailable):
        load_selected_distance_csv("test", {"all": root / "all.csv",
                                             "all_ok_3": root / "missing.csv"},
                                   "all_ok_3", fallback_to_all=False)


def test_g402s_figure_uses_all_prespecified_matched_ca_pairs_without_categories():
    root = data_access.REPO_ROOT
    effects, source, figure = plot_g402s_paired_revision(root)
    partners = [
        "PHE1519", "VAL1520", "ALA1521", "VAL1522", "ILE1523", "MET1524",
        "ASP1525", "ASN1526", "PHE1527", "ASP1528", "TYR1529", "LEU1530",
        "THR1531", "ARG1532", "ASP1533", "TRP1534", "SER1535",
    ]
    assert sorted(effects.partner.unique()) == sorted(partners)
    assert sorted(source.partner.unique()) == sorted(partners)
    assert set(effects.contrast) == {
        "vanilla G402S - vanilla WT", "G402S masked - G402S vanilla"
    }
    assert not any(column.startswith("shortest_") for column in source.source_column)
    figure.clear()


def test_r406_screen_is_named_as_screen_and_shown_with_conditional_contacts():
    summary, figure = plot_g406r_screen_and_contacts(data_access.REPO_ROOT)
    assert set(summary.partner_residue) == {"screen_only", "ASP1528", "ASP1533"}
    assert set(summary.atom_selection_rule) == {"side_chain_non_hydrogen_atoms"}
    screen = summary[summary.partner_residue.eq("screen_only")].set_index("protocol")
    assert screen.loc["vanilla", "unconditional_pass_fraction"] == pytest.approx(2419 / 5000)
    assert screen.loc["masked", "unconditional_pass_fraction"] == pytest.approx(835 / 4997)
    contacts = summary[summary.partner_residue.eq("ASP1528")].set_index("protocol")
    assert contacts.loc["vanilla", "conditional_contact_denominator"] == 2419
    assert "does not validate the entire pocket" in " ".join(
        text.get_text() for text in figure.texts
    )
    figure.clear()
