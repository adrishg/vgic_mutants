"""Explicit selection of original or RMSD-filtered distance CSV variants."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


VALID_SELECTIONS = (
    "all", "all_ok", "all_ok_3", "all_ok_3_structural_qc",
    "all_ok_3_structural_interface_qc",
    "all_ok_3_structural_interface_alignment_qc", "all_ok_3p5",
    "first_converged", "first_100_generated",
)

KV21_ACTIVE_MASKED_MANIFEST = "analysis/kv21_active_masked_cohort.tsv"
ACTIVE_MASKED_MANIFEST = "analysis/active_masked_cohorts.tsv"


def _active_masked_rows(repo_root: str | Path) -> pd.DataFrame:
    return pd.read_csv(Path(repo_root) / ACTIVE_MASKED_MANIFEST, sep="\t")


def active_masked_distance_path(
    repo_root: str | Path,
    channel: str,
    condition: str,
    protocol: str = "masked",
) -> Path:
    """Return the registered current masked distance path.

    Legacy masked paths are intentionally absent from this registry.  Missing
    corrected cohorts raise explicitly so a notebook cannot silently analyze
    an old masked dataset.
    """
    channel_aliases = {"kv21": "Kv2.1", "nav15": "NaV1.5", "cav12": "CaV1.2"}
    condition_aliases = {
        "wt": "WT", "l403a": "L403A", "f412l": "F412L",
        "qqq": "QQQ", "g402s": "G402S", "g406r": "G406R",
    }
    channel_key = channel_aliases.get(str(channel).lower(), str(channel))
    condition_key = condition_aliases.get(str(condition).lower(), str(condition))
    if channel_key == "Kv2.1" and protocol == "masked":
        protocol = "masked_v3_50pct"
    if channel_key == "CaV1.2" and protocol == "masked":
        protocol = "masked_v3_50pct"
    table = _active_masked_rows(repo_root)
    rows = table.loc[
        (table["channel"] == channel_key)
        & (table["condition"] == condition_key)
        & (table["protocol"] == protocol)
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one active masked row for {channel_key}/{condition_key}/{protocol}; "
            f"found {len(rows)} in {Path(repo_root) / ACTIVE_MASKED_MANIFEST}"
        )
    row = rows.iloc[0]
    if row["status"] not in {"active", "active_new", "active_corrected", "candidate_new"} or not str(row["distance_relpath"]).strip():
        raise FileNotFoundError(
            f"No corrected masked cohort is available for "
            f"{channel_key}/{condition_key}/{protocol}; status={row['status']}"
        )
    return Path(str(row["distance_relpath"]))


def nav15_active_masked_distance_path(
    repo_root: str | Path, condition: str, variant: str = "v1"
) -> Path:
    """Return corrected NaV1.5 masked v1 or v2 distance input."""
    variant_key = str(variant).lower().replace("-", "_")
    if variant_key in {"v2_noifm", "masked_v2_noifm"}:
        protocol = "masked_v2_noIFM"
    elif variant_key in {"v2", "masked_v2"}:
        protocol = "masked_v2"
    else:
        protocol = "masked_v1"
    return active_masked_distance_path(repo_root, "NaV1.5", condition, protocol)


def cav12_active_masked_distance_path(
    repo_root: str | Path, condition: str = "WT"
) -> Path:
    """Return the registered corrected CaV1.2 50% masked distance path."""
    return active_masked_distance_path(repo_root, "CaV1.2", condition, "masked_v3_50pct")


def kv21_active_masked_distance_path(
    repo_root: str | Path, condition: str
) -> Path:
    """Return the registered current Kv2.1 masked distance-table path.

    The small manifest is deliberately versioned in the repository; the large
    CSV itself is resolved by ``shared.data_access`` from the external data
    package. This keeps every notebook on the same complete corrected 50% v3
    cohort without copying the data into Git. Legacy masked cohorts and the
    alternate legacy cohorts are not active here.
    """
    aliases = {"wt": "WT", "l403a": "L403A", "f412l": "F412L"}
    condition_key = aliases.get(str(condition).lower(), str(condition))
    manifest = Path(repo_root) / KV21_ACTIVE_MASKED_MANIFEST
    table = pd.read_csv(manifest, sep="\t")
    rows = table.loc[
        (table["channel"] == "Kv2.1")
        & (table["condition"] == condition_key)
        & (table["protocol"] == "masked_v3_50pct")
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one active Kv2.1 masked row for {condition_key!r}; "
            f"found {len(rows)} in {manifest}"
        )
    return Path(str(rows.iloc[0]["distance_relpath"]))


def _kv21_active_manifest_row(repo_root: str | Path, condition: str) -> pd.Series:
    aliases = {"wt": "WT", "l403a": "L403A", "f412l": "F412L"}
    condition_key = aliases.get(str(condition).lower(), str(condition))
    manifest = Path(repo_root) / KV21_ACTIVE_MASKED_MANIFEST
    table = pd.read_csv(manifest, sep="\t")
    rows = table.loc[
        (table["channel"] == "Kv2.1")
        & (table["condition"] == condition_key)
        & (table["protocol"] == "masked_v3_50pct")
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one active Kv2.1 masked row for {condition_key!r}; "
            f"found {len(rows)} in {manifest}"
        )
    return rows.iloc[0]


def kv21_active_rmsd_path(repo_root: str | Path, condition: str = "WT") -> Path:
    """Return the registered combined RMSD table for the corrected cohort."""
    return Path(str(_kv21_active_manifest_row(repo_root, condition)["rmsd_relpath"]))


def kv21_active_contact_path(repo_root: str | Path, condition: str = "WT") -> Path:
    """Return the registered contact table for the corrected cohort."""
    return Path(str(_kv21_active_manifest_row(repo_root, condition)["contact_relpath"]))


def _read_csv_resolving_local_lfs(path: Path) -> pd.DataFrame:
    """Read a CSV from the working tree, the local cache, or Hugging Face."""
    from . import data_access
    from .data_access import DataUnavailable

    return pd.read_csv(data_access.resolve(path))


def _companion_extra_path(path: Path) -> Path | None:
    """Return the established companion-distance filename for a CSV variant."""
    name = path.name
    if not name.endswith(".csv"):
        return None
    stem = name[:-4]
    if "Nav15" in name and stem.endswith("_distances"):
        return path.with_name(stem + "_extra_ifm.csv")
    if "_all_ok_rmsd_3A" in stem:
        return path.with_name(stem.replace("_all_ok_rmsd_3A", "_extra_all_ok_rmsd_3A", 1) + ".csv")
    if stem.endswith("_distances"):
        return path.with_name(stem + "_extra.csv")
    return None


def _merge_companion_extra(path: Path, frame: pd.DataFrame) -> pd.DataFrame:
    """Join the established Kv2.1/CaV1.2 companion distance table when present.

    The main and ``_extra`` exports are complementary distance sets with the
    same one-row-per-model key. This keeps the large data external while
    making notebook inputs match the reconciliation workflow.
    """
    if "pdb_file" not in frame.columns:
        return frame
    companion = _companion_extra_path(path)
    if companion is None:
        return frame
    from . import data_access
    from .data_access import DataUnavailable
    try:
        companion_resolved = data_access.resolve(companion)
    except (DataUnavailable, FileNotFoundError):
        return frame
    extra = pd.read_csv(companion_resolved, low_memory=False)
    if "pdb_file" not in extra.columns:
        return frame
    if extra["pdb_file"].duplicated().any():
        raise ValueError(f"Companion distance table has duplicate pdb_file keys: {companion}")
    duplicate = [column for column in extra.columns if column in frame.columns and column != "pdb_file"]
    extra = extra.drop(columns=duplicate)
    merged = frame.merge(extra, on="pdb_file", how="left", validate="one_to_one")
    print(f"Merged companion distance table: {companion}")
    return merged


def select_manifest_rows(path: str | Path, subset: str) -> pd.DataFrame:
    """Select recycle rows from an RMSD convergence manifest."""
    frame = pd.read_csv(path)
    if subset == "all":
        return frame[frame["parse_ok"].fillna(False)].copy()
    if subset == "all_ok":
        return frame[frame["all_ok"].fillna(False)].copy()
    if subset == "all_ok_3":
        selected = []
        work = frame[~frame["is_base"].fillna(False)].copy()
        for _, trajectory in work.groupby(["seed", "model_number"], sort=False):
            trajectory = trajectory.sort_values("recycle_number")
            passed = (
                pd.to_numeric(
                    trajectory["rmsd_to_previous_available"], errors="coerce"
                ).le(3.0)
                & pd.to_numeric(
                    trajectory["aligned_coverage_to_previous"], errors="coerce"
                ).ge(0.9)
            )
            starts = [
                index
                for index in range(len(trajectory))
                if bool(passed.iloc[index:].all())
            ]
            if starts:
                selected.append(trajectory.iloc[starts[0]:])
        return (
            pd.concat(selected, ignore_index=True)
            if selected
            else work.iloc[0:0].copy()
        )
    if subset == "first_converged":
        return frame[frame["earliest_converged_selected"].fillna(False)].copy()
    if subset == "first_100_generated":
        valid = frame[frame["all_ok"].fillna(False)].copy()
        trajectories = (
            valid[["seed", "model_number"]]
            .drop_duplicates()
            .sort_values(["seed", "model_number"])
            .head(100)
        )
        return valid.merge(trajectories, on=["seed", "model_number"], how="inner")
    raise ValueError(f"Unknown subset: {subset}")


def distance_csv_options(
    repo_root: str | Path,
    original_path: str | Path,
    channel: str,
    condition: str,
    protocol: str,
) -> dict[str, Path]:
    """Return all three explicit paths without reading any CSV."""
    root = Path(repo_root)
    original = Path(original_path)
    base = original.stem
    # Some distance exports use ``_all_all`` to distinguish the complete
    # recalculated table. Dataset variants should still receive one clean suffix.
    while base.endswith("_all"):
        base = base[:-4]
    channel_dir = {
        "cav12": "cav12",
        "kv21": "kv21",
        "nav15": "nav15",
    }.get(str(channel).lower())
    if channel_dir is None:
        raise ValueError(f"Unknown channel directory for {channel!r}")
    filtered = root / channel_dir / "rmsd_filtered_distances" / condition / protocol
    threshold_3p5 = (
        root / channel_dir / "rmsd_threshold_sensitivity"
        / "3p5A" / condition / protocol
    )
    return {
        "all": original,
        "all_ok": filtered / f"{base}_all_ok_5.csv",
        # Keep the publication-facing 3 Å tables beside their source CSVs so
        # notebooks cloned from GitHub can switch datasets without recreating
        # the repository's intermediate threshold-sensitivity directory.
        "all_ok_3": original.parent / f"{base}_all_ok_rmsd_3A.csv",
        "all_ok_3_structural_qc": original.parent / f"{base}_all_ok_rmsd_3A_structural_qc.csv",
        "all_ok_3_structural_interface_qc": original.parent / f"{base}_all_ok_rmsd_3A_structural_interface_qc.csv",
        "all_ok_3_structural_interface_alignment_qc": original.parent / f"{base}_all_ok_rmsd_3A_structural_interface_alignment_qc.csv",
        "all_ok_3p5": threshold_3p5 / f"{base}_all_ok_rmsd_3p5A.csv",
        "first_converged": filtered / f"{base}_earliest_converged.csv",
        "first_100_generated": filtered / f"{base}_first_100_generated.csv",
    }


def load_selected_distance_csv(
    label: str,
    options: dict[str, Path],
    selection: str,
    *,
    fallback_to_all: bool = True,
) -> pd.DataFrame:
    """Load the requested variant, visibly falling back only when authorized."""
    if selection not in VALID_SELECTIONS:
        raise ValueError(f"selection must be one of {VALID_SELECTIONS}; got {selection!r}")
    selected = options[selection]
    actual = selection
    from .data_access import DataUnavailable

    try:
        frame = _read_csv_resolving_local_lfs(selected)
        frame = _merge_companion_extra(selected, frame)
    except (DataUnavailable, FileNotFoundError) as error:
        if not fallback_to_all:
            if isinstance(error, FileNotFoundError):
                raise FileNotFoundError(
                    f"{label}: requested {selection!r} file does not exist: {selected}"
                ) from error
            raise
        if selection == "all":
            raise
        selected = options["all"]
        actual = "all"
        print(f"WARNING: {label}: {selection!r} is unavailable; using original 'all': {selected}")
        frame = _read_csv_resolving_local_lfs(selected)
        frame = _merge_companion_extra(selected, frame)
    print(f"{label}: requested={selection}; actual={actual}; path={selected}")
    frame.attrs["dataset_label"] = label
    frame.attrs["requested_selection"] = selection
    frame.attrs["actual_selection"] = actual
    frame.attrs["source_path"] = str(selected)
    return frame


def apply_kv21_interface_qc(
    frame: pd.DataFrame,
    *,
    threshold_A: float = 27.0,
    trajectory_level: bool = True,
) -> pd.DataFrame:
    """Remove Kv2.1 trajectories with a detached pore–VSD interface.

    The existing structural-QC CSVs verify the G377 model-numbered
    selectivity-filter ring. A converged model can still retain that ring while
    separating the S6/C-terminal end of the pore from the neighboring
    voltage-sensor region. The six K429/E425/K422-to-N181/V184 contacts in each
    chain provide a direct check for this failure.

    Grossly failed trajectories begin above approximately 32.9 Å. One
    additional L403A vanilla recycle reaches 27.35 Å even though every other
    recycle in that trajectory remains near 20–21 Å. The default 27 Å cutoff
    removes this isolated excursion while remaining above the largest
    coordinate-derived 8SD3/8SDA reference distance (25.51 Å).
    """
    residue_pairs = (
        ("LYS429", "ASN181"), ("LYS429", "VAL184"),
        ("GLU425", "ASN181"), ("GLU425", "VAL184"),
        ("LYS422", "ASN181"), ("LYS422", "VAL184"),
    )
    columns = [
        f"CA_CA_{chain}_{first}_CA-{chain}_{second}_CA"
        for chain in "ABCD"
        for first, second in residue_pairs
    ]
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(f"Kv2.1 interface QC requires 24 distance columns; missing {missing}")
    if "pdb_file" not in frame.columns:
        raise KeyError("Kv2.1 interface QC requires the 'pdb_file' column")

    interface_max = frame[columns].apply(pd.to_numeric, errors="coerce").max(axis=1)
    directly_failed = interface_max.gt(threshold_A)
    trajectory = frame["pdb_file"].astype(str).str.replace(
        r"\.r\d+\.pdb$", ".pdb", regex=True
    )
    failed_trajectories = set(trajectory[directly_failed])
    rejected = trajectory.isin(failed_trajectories) if trajectory_level else directly_failed

    result = frame.loc[~rejected].copy()
    result.attrs.update(frame.attrs)
    result.attrs["kv21_interface_qc_threshold_A"] = float(threshold_A)
    result.attrs["kv21_interface_qc_trajectory_level"] = bool(trajectory_level)
    result.attrs["kv21_interface_qc_rejected_rows"] = int(rejected.sum())
    result.attrs["kv21_interface_qc_rejected_trajectories"] = len(failed_trajectories)
    result.attrs["kv21_interface_qc_max_retained_A"] = float(
        result[columns].apply(pd.to_numeric, errors="coerce").max(axis=1).max()
    )

    label = frame.attrs.get("dataset_label", "Kv2.1 ensemble")
    print(
        f"{label}: pore–VSD interface QC retained {len(result)}/{len(frame)} rows; "
        f"rejected {int(rejected.sum())} rows from {len(failed_trajectories)} trajectories; "
        f"maximum retained interface distance={result.attrs['kv21_interface_qc_max_retained_A']:.2f} Å"
    )
    return result
