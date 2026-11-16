#!/usr/bin/env python3
"""Validate and publish the current Kv2.1 masked 50% derived analyses.

The coordinate shards are already published separately.  This script packages
the derived distance tables, convergence manifests, and publication-facing
analysis tables without modifying or deleting any existing dataset files.

By default it performs a read-only validation and prints the upload plan.  Use
``--stage-dir`` to materialize the generated ``all_ok_rmsd_3A`` tables, and
``--upload`` only after reviewing the plan.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
import sys
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from shared.dataset_selection import select_manifest_rows


DEFAULT_SOURCE = Path("kv21_50pct_reanalysis_2026-11-12")
DEFAULT_REPO = "adrishgz/vgic-mutant-structural-ensembles"
RELEASE_TAG = "50pct_reanalysis_2026-11-12"

COHORTS = {
    "wt": {"label": "WT", "source_key": "wt_masked_v3_50pct"},
    "l403a": {"label": "L403A", "source_key": "l403a_masked_v3_50pct"},
    "f412l": {"label": "F412L", "source_key": "f412l_masked_v3_50pct"},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def distance_paths(source_root: Path, key: str) -> tuple[Path, Path]:
    directory = source_root / "distances" / COHORTS[key]["source_key"]
    main = directory / f"Kv21_{key}_masked_v3_50pct_distances_all.csv"
    extra = directory / f"Kv21_{key}_masked_v3_50pct_distances_extra.csv"
    return main, extra


def manifest_path(source_root: Path, key: str) -> Path:
    return source_root / "convergence" / f"Kv21_rmsd_convergence_{key}_masked" / "all_models_manifest.csv"


def selected_distance_table(main: Path, manifest: Path, output: Path) -> dict[str, int | str]:
    selected = select_manifest_rows(manifest, "all_ok_3")
    if "pdb_basename" not in selected.columns:
        raise ValueError(f"{manifest} lacks pdb_basename")
    selected_keys = selected["pdb_basename"].astype(str).map(lambda value: Path(value).name)
    selected_keys = set(selected_keys)

    header = pd.read_csv(main, nrows=0)
    if "pdb_file" not in header.columns:
        raise ValueError(f"{main} lacks pdb_file")
    frame = pd.read_csv(main, low_memory=False)
    frame_keys = frame["pdb_file"].astype(str).map(lambda value: Path(value).name)
    filtered = frame.loc[frame_keys.isin(selected_keys)].copy()
    if len(filtered) != len(selected):
        missing = selected_keys.difference(frame_keys)
        raise ValueError(
            f"Manifest-to-distance join failed for {main.name}: "
            f"selected={len(selected)} joined={len(filtered)} missing={len(missing)}"
        )
    filtered.to_csv(output, index=False)
    return {
        "manifest_rows": int(len(selected)),
        "distance_rows": int(len(frame)),
        "all_ok_rmsd_3A_rows": int(len(filtered)),
        "sha256": sha256(output),
    }


def check_source(source_root: Path) -> dict[str, dict[str, int | str]]:
    summary: dict[str, dict[str, int | str]] = {}
    for key in COHORTS:
        main, extra = distance_paths(source_root, key)
        manifest = manifest_path(source_root, key)
        for required in (main, extra, manifest):
            if not required.is_file():
                raise FileNotFoundError(required)
        header = pd.read_csv(main, nrows=0)
        required_columns = {
            "pdb_file",
            "CA_CA_A_GLU425_CA-A_ASN181_CA",
            "CA_CA_B_GLU425_CA-B_ASN181_CA",
            "CA_CA_C_GLU425_CA-C_ASN181_CA",
            "CA_CA_D_GLU425_CA-D_ASN181_CA",
        }
        missing = required_columns.difference(header.columns)
        if missing:
            raise ValueError(f"{main} is missing required columns: {sorted(missing)}")
        rows = sum(1 for _ in main.open()) - 1
        manifest_rows = sum(1 for _ in manifest.open()) - 1
        if rows != 6000 or manifest_rows != 6000:
            raise ValueError(f"{key}: expected 6000 rows, got distance={rows}, manifest={manifest_rows}")
        summary[key] = {
            "distance_rows": rows,
            "manifest_rows": manifest_rows,
            "distance_columns": len(header.columns),
        }
    return summary


def add_upload(plan: list[tuple[Path, str]], source: Path, destination: str) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    plan.append((source, destination))


def public_source_label(path: Path, source_root: Path) -> str:
    for root, label in ((source_root, source_root.name), (ROOT, "repository")):
        try:
            return str(Path(label) / path.relative_to(root))
        except ValueError:
            continue
    return path.name


def build_plan(source_root: Path, stage_root: Path) -> tuple[list[tuple[Path, str]], dict]:
    source_summary = check_source(source_root)
    plan: list[tuple[Path, str]] = []
    release_summary = {"source": source_summary, "cohorts": {}}

    for key, metadata in COHORTS.items():
        main, extra = distance_paths(source_root, key)
        manifest = manifest_path(source_root, key)
        prefix = f"data/Kv2.1/distances/26-11-12_Kv2.1_{key}_masked_v3_50pct_AF2"
        add_upload(plan, main, f"{prefix}_distances.csv")
        add_upload(plan, extra, f"{prefix}_distances_extra.csv")

        staged_main = stage_root / f"26-11-12_Kv2.1_{key}_masked_v3_50pct_AF2_distances_all_ok_rmsd_3A.csv"
        stats = selected_distance_table(main, manifest, staged_main)
        staged_extra = stage_root / f"26-11-12_Kv2.1_{key}_masked_v3_50pct_AF2_distances_extra_all_ok_rmsd_3A.csv"
        extra_frame = pd.read_csv(extra, low_memory=False)
        selected_keys = set(pd.read_csv(staged_main, usecols=["pdb_file"])["pdb_file"].astype(str))
        extra_frame = extra_frame[extra_frame.pdb_file.astype(str).isin(selected_keys)]
        if len(extra_frame) != stats["all_ok_rmsd_3A_rows"]:
            raise ValueError(f"{key}: filtered extra table does not join one-to-one")
        extra_frame.to_csv(staged_extra, index=False)
        stats["all_ok_rmsd_3A_extra_rows"] = int(len(extra_frame))
        add_upload(plan, staged_main, f"{prefix}_distances_all_ok_rmsd_3A.csv")
        add_upload(plan, staged_extra, f"{prefix}_distances_extra_all_ok_rmsd_3A.csv")

        convergence_dir = source_root / "convergence" / f"Kv21_rmsd_convergence_{key}_masked"
        destination_dir = f"data/Kv2.1/convergence/Kv21_rmsd_convergence_{key}_masked_v3_50pct"
        convergence_files = []
        for path in sorted(convergence_dir.iterdir()):
            if path.suffix not in {".csv", ".json"}:
                continue
            add_upload(plan, path, f"{destination_dir}/{path.name}")
            convergence_files.append(path.name)

        release_summary["cohorts"][key] = {
            "label": metadata["label"],
            **stats,
            "convergence_files": convergence_files,
        }

    reconciliation = source_root / "comparisons" / "reconciliation"
    for path in sorted(reconciliation.iterdir()):
        if path.suffix not in {".csv", ".gz", ".json"}:
            continue
        add_upload(
            plan,
            path,
            f"data/Kv2.1/extra/conformation_analysis/{RELEASE_TAG}/{path.name}",
        )

    rmsf_tables = source_root / "figure_data" / "rmsf" / "tables"
    for path in sorted(rmsf_tables.glob("*.csv")):
        add_upload(plan, path, f"data/Kv2.1/rmsf/analysis/{RELEASE_TAG}/{path.name}")
    for path in sorted((source_root / "figures" / "rmsf").glob("*.png")):
        add_upload(plan, path, f"data/Kv2.1/rmsf/figures/{RELEASE_TAG}/{path.name}")
    for path in sorted((reconciliation / "figures").glob("*.png")):
        add_upload(plan, path, f"data/Kv2.1/extra/figures/{RELEASE_TAG}/{path.name}")
    for path in sorted((reconciliation / "figures").glob("*.pdf")):
        add_upload(plan, path, f"data/Kv2.1/extra/figures/{RELEASE_TAG}/{path.name}")

    # These outputs are generated by the repository notebooks after the
    # current tables have been staged.  Include them when present so the
    # dataset carries the exact supplementary artifacts used for publication.
    for figure_id in ("s9", "s10"):
        figure_dir = ROOT / "docs" / "figures" / f"supplementary_figure_{figure_id}"
        table_dir = ROOT / "docs" / "tables" / f"supplementary_figure_{figure_id}"
        for path in sorted(figure_dir.glob("Figure_*.png")) + sorted(figure_dir.glob("Figure_*.pdf")):
            add_upload(plan, path, f"data/Kv2.1/extra/figures/{RELEASE_TAG}/{path.name}")
        for path in sorted(table_dir.glob("*.csv")):
            add_upload(plan, path, f"data/Kv2.1/extra/figure_tables/{figure_id}/{path.name}")

    manifest = {
        "release": RELEASE_TAG,
        "channel": "Kv2.1",
        "protocol": "masked_v3_50pct",
        "source_generated": "2026-11-12",
        "selection": "all_ok_rmsd_3A; successive available-recycle RMSD <= 3.0 A and aligned coverage >= 0.90",
        "structural_qc_status": "not silently claimed; convergence-qualified release",
        "source_root": source_root.name,
        "files": [
            {
                "source": public_source_label(source, source_root),
                "destination": destination,
                "bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
            for source, destination in plan
        ],
        "cohorts": release_summary["cohorts"],
    }
    manifest_path_out = stage_root / "kv21_50pct_derived_release_manifest.json"
    manifest_path_out.write_text(json.dumps(manifest, indent=2) + "\n")
    add_upload(plan, manifest_path_out, f"metadata/{manifest_path_out.name}")
    return plan, manifest


def upload(plan: Iterable[tuple[Path, str]], repo_id: str) -> None:
    from huggingface_hub import HfApi

    api = HfApi()
    for source, destination in plan:
        print(f"UPLOAD {source} -> {repo_id}:{destination}")
        api.upload_file(
            path_or_fileobj=str(source),
            path_in_repo=destination,
            repo_id=repo_id,
            repo_type="dataset",
            commit_message=f"Add Kv2.1 masked 50% derived analyses ({RELEASE_TAG})",
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--repo-id", default=DEFAULT_REPO)
    parser.add_argument("--stage-dir", type=Path)
    parser.add_argument("--upload", action="store_true")
    args = parser.parse_args()

    if not args.source_root.is_dir():
        raise SystemExit(f"Source root does not exist: {args.source_root}")
    temporary = None
    if args.stage_dir is None:
        temporary = tempfile.TemporaryDirectory(prefix="kv21_50pct_release_")
        stage_root = Path(temporary.name)
    else:
        stage_root = args.stage_dir
        stage_root.mkdir(parents=True, exist_ok=True)

    plan, manifest = build_plan(args.source_root, stage_root)
    print(json.dumps({"cohorts": manifest["cohorts"], "planned_files": len(plan)}, indent=2))
    for source, destination in plan:
        print(f"PLAN {source} -> {destination}")
    if args.upload:
        upload(plan, args.repo_id)
        print(f"Uploaded {len(plan)} files to {args.repo_id}")
    elif temporary is not None:
        print("Validation complete; generated staging files were temporary. Use --stage-dir to retain them.")


if __name__ == "__main__":
    main()
