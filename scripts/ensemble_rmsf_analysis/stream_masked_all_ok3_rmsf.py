#!/usr/bin/env python3
"""Final-QC (All-OK-3) RMSF for the 50% masked CaV1.2 and NaV1.5 cohorts.

Masked models are read one at a time from the published ``tar.zst`` coordinate
shards (no extraction to disk). Each model is mapped to the raw residue axis
and fitted to the anchor reference core with the channel package's own
functions (``load_structure``, ``chain_to_raw_ca``, ``kabsch_fit``), using the
50% Hive reference arrays, so the frame is the one used for the published
all-model profiles. Only snapshots that pass the All-OK-3 convergence rule
(``build_all_ok3_manifests.select_all_ok3`` on the Hive convergence manifest)
enter the moments.

Vanilla cohorts are unchanged between releases; their All-OK-3 rows come from
the existing final-QC vanilla profile (same frame: median mean-coordinate
shift versus the 50% all-model run is 0.05 Å, explained by r0 inclusion).
Residue annotations (mask classes, identities) come from the 50% Hive
all-model profile.

Output: ``$VGIC_DATA_ROOT/<channel>/dataRMSF/profiles/<channel>_all_ok_3_per_residue_profiles.csv``
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tarfile
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from ensemble_rmsf_analysis.build_all_ok3_manifests import select_all_ok3  # noqa: E402
from ensemble_rmsf_analysis.recompute_all_ok3_local import add_reference_metrics  # noqa: E402

HF_REPO = "adrishgz/vgic-mutant-structural-ensembles"

CHANNELS = {
    "cav12": {
        "package": "cav12_ensemble_rmsf_cluster_package",
        "module": "cav12_ensemble_rmsf",
        "hive": "Cav1.2",
        "min_matched": 850,
        "min_coverage": 0.80,
        "cohorts": [
            ("wt_masked", "wt", "wt_masked_v3_50pct", "coordinates/CaV1.2/WT/masked/50pct", "Cav12_rmsd_convergence_wt_masked_v3_50pct"),
            ("g402s_masked", "g402s", "g402s_masked_v3_50pct", "coordinates/CaV1.2/G402S/masked/50pct", "Cav12_rmsd_convergence_g402s_masked_v3_50pct"),
            ("g406r_masked", "g406r", "g406r_masked_v3_50pct", "coordinates/CaV1.2/G406R/masked/50pct", "Cav12_rmsd_convergence_g406r_masked_v3_50pct"),
        ],
        "vanilla": ["wt_vanilla", "g402s_vanilla", "g406r_vanilla"],
    },
    "nav15": {
        "package": "nav15_ensemble_rmsf_cluster_package",
        "module": "nav15_ensemble_rmsf",
        "hive": "Nav1.5",
        "min_matched": 700,
        "min_coverage": 0.75,
        "cohorts": [
            ("wt_masked", "wt", "wt_masked_v1_50pct", "coordinates/NaV1.5/WT/masked/v1_50pct", "Nav15_rmsd_convergence_wt_masked_v1_50pct"),
            ("wt_masked_v2", "wt", "wt_masked_v2_50pct", "coordinates/NaV1.5/WT/masked/v2_50pct", "Nav15_rmsd_convergence_wt_masked_v2_50pct"),
            ("qqq_masked", "qqq", "qqq_masked_v1_50pct", "coordinates/NaV1.5/QQQ/masked/v1_50pct", "Nav15_rmsd_convergence_qqq_masked_v1_50pct"),
            ("qqq_masked_v2", "qqq", "qqq_masked_v2_50pct", "coordinates/NaV1.5/QQQ/masked/v2_50pct", "Nav15_rmsd_convergence_qqq_masked_v2_50pct"),
        ],
        "vanilla": ["wt_vanilla", "qqq_vanilla"],
    },
}

METRIC_PREFIXES = (
    "number_of_models", "number_with_residue_resolved", "coverage_fraction",
    "ensemble_mean_", "ensemble_rmsf_A", "mean_distance_to_", "rms_deviation_to_",
    "number_comparable_to_", "coverage_comparable_to_",
)


def allowlist(hive: Path, convergence_dir: str) -> set[str]:
    manifest = pd.read_csv(hive / "rmsd_convergence_filtering" / convergence_dir / "all_models_manifest.csv", low_memory=False)
    keep = manifest.loc[select_all_ok3(manifest)]
    column = "pdb_basename" if "pdb_basename" in keep else "pdb_file"
    return set(keep[column].astype(str).map(lambda p: Path(p).name))


def archive_members(shard_dir: str, allowed_sources: set[str], files: list[str]) -> set[str]:
    """Translate source PDB basenames to the renamed members inside the release shards."""
    from huggingface_hub import hf_hub_download

    for manifest in (f for f in files if f.startswith("metadata/coordinate_pdb_release/") and f.endswith("model_manifest.jsonl")):
        cohort = json.loads(Path(hf_hub_download(HF_REPO, manifest.replace("model_manifest.jsonl", "cohort_manifest.json"), repo_type="dataset")).read_text())
        if cohort.get("release_path") != shard_dir:
            continue
        table = pd.read_json(hf_hub_download(HF_REPO, manifest, repo_type="dataset"), lines=True)
        selected = table.loc[table.source_filename.isin(allowed_sources), "archive_member"]
        return set(selected.map(lambda p: Path(p).name))
    raise FileNotFoundError(f"No release manifest for {shard_dir}")


def process_shard(channel: str, shard: str, allowed: set[str], references_path: str):
    cfg = CHANNELS[channel]
    sys.path.insert(0, str(REPO / "scripts" / cfg["package"] / "src"))
    core = __import__(f"{cfg['module']}.core", fromlist=["load_structure"])
    from huggingface_hub import hf_hub_download

    local = hf_hub_download(HF_REPO, shard, repo_type="dataset")
    refs = np.load(references_path, allow_pickle=False)
    ids = refs["reference_ids"].astype(str).tolist()
    anchor = ids.index(str(refs["anchor_id"].item()))
    anchor_coords = refs["coords"][anchor].astype(np.float64)
    anchor_present = refs["present"][anchor].astype(bool)
    core_mask = refs["core_mask"].astype(bool)
    raw = refs["raw_residue_numbers"].astype(int)
    n = len(raw)
    counts = np.zeros(n, np.int64)
    sums = np.zeros((n, 3))
    sum_sq = np.zeros(n)
    used, failed = [], []
    proc = subprocess.Popen(["zstd", "-dc", str(Path(local).resolve())], stdout=subprocess.PIPE)
    with tempfile.TemporaryDirectory() as tmp, tarfile.open(fileobj=proc.stdout, mode="r|") as tar:
        for member in tar:
            name = Path(member.name).name
            if not member.isfile() or name not in allowed:
                continue
            path = Path(tmp) / name
            path.write_bytes(tar.extractfile(member).read())
            try:
                structure = core.load_structure(path)
                chain = structure.chains["A"]
                coords, present, _, _, _ = core.chain_to_raw_ca(
                    structure, "A", chain.sequence, int(raw[0]), int(raw[-1]), {}
                )
                valid = core_mask & present & anchor_present
                coverage = valid.sum() / (core_mask & anchor_present).sum()
                if valid.sum() < cfg["min_matched"] or coverage < cfg["min_coverage"]:
                    raise ValueError(f"core coverage {coverage:.3f}")
                xyz = core.apply_fit(coords, core.kabsch_fit(coords[valid], anchor_coords[valid]))
                ok = present & np.isfinite(xyz).all(axis=1)
                safe = np.where(ok[:, None], xyz, 0.0)
                counts += ok
                sums += safe
                sum_sq += np.einsum("ij,ij->i", safe, safe)
                used.append(name)
            except Exception as exc:  # recorded, never imputed
                failed.append((name, f"{type(exc).__name__}: {exc}"))
            finally:
                path.unlink(missing_ok=True)
    proc.wait()
    return counts, sums, sum_sq, used, failed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel", choices=sorted(CHANNELS), required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--legacy-vanilla", required=True,
                        help="Existing final-QC profile providing the unchanged vanilla rows.")
    args = parser.parse_args()
    cfg = CHANNELS[args.channel]
    root = Path(os.environ["VGIC_DATA_ROOT"])
    hive = root / "hive" / cfg["hive"]
    references_path = hive / "dataRMSF" / "references" / f"{args.channel}_aligned_references.npz"
    refs = np.load(references_path, allow_pickle=False)
    template = pd.read_csv(hive / "dataRMSF" / "profiles" / f"{args.channel}_all_models_per_residue_profiles.csv")
    legacy = pd.read_csv(args.legacy_vanilla)
    if os.environ.get("HF_HUB_OFFLINE"):
        # Deterministic release layout: five shards per cohort plus one
        # model manifest per cohort (pre-downloaded into the HF cache).
        releases = {
            "cav12": ["cav12_wt_masked_50pct", "cav12_g402s_masked_50pct", "cav12_g406r_masked_50pct"],
            "nav15": ["nav15_wt_masked_v1_50pct", "nav15_wt_masked_v2_50pct",
                      "nav15_qqq_masked_v1_50pct", "nav15_qqq_masked_v2_50pct"],
        }[args.channel]
        files = [f"{c[3]}/shard-{i:03d}.tar.zst" for c in cfg["cohorts"] for i in range(5)]
        files += [f"metadata/coordinate_pdb_release/{r}/model_manifest.jsonl" for r in releases]
    else:
        from huggingface_hub import HfApi
        files = HfApi().list_repo_files(HF_REPO, repo_type="dataset")

    jobs, allowed_by = [], {}
    for out_name, _, hive_name, shard_dir, conv in cfg["cohorts"]:
        sources = allowlist(hive, conv)
        allowed_by[out_name] = archive_members(shard_dir, sources, files)
        print(f"{out_name}: {len(sources)} All-OK-3 snapshots; {len(allowed_by[out_name])} found in release", flush=True)
        for shard in sorted(f for f in files if f.startswith(shard_dir + "/") and f.endswith(".tar.zst")):
            jobs.append((out_name, shard))
    totals = {name: [0, 0, 0, [], []] for name in allowed_by}
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(process_shard, args.channel, shard, allowed_by[name], str(references_path)): name
                   for name, shard in jobs}
        for future in as_completed(futures):
            name = futures[future]
            c, s, q, used, failed = future.result()
            t = totals[name]
            t[0] = t[0] + c; t[1] = t[1] + s; t[2] = t[2] + q
            t[3] += used; t[4] += failed
            print(f"{name}: +{len(used)} models ({len(failed)} failed)", flush=True)

    rows, summary = [], {}
    for out_name, condition, hive_name, _, _ in cfg["cohorts"]:
        counts, sums, sum_sq, used, failed = totals[out_name]
        mean = np.full(sums.shape, np.nan)
        np.divide(sums, counts[:, None], out=mean, where=counts[:, None] > 0)
        mean_sq = np.full(sum_sq.shape, np.nan)
        np.divide(sum_sq, counts, out=mean_sq, where=counts > 0)
        rmsf = np.sqrt(np.maximum(mean_sq - np.einsum("ij,ij->i", mean, mean), 0.0))
        frame = template.loc[template.dataset.eq(hive_name)].sort_values("raw_residue_number").reset_index(drop=True)
        frame["dataset"] = out_name
        frame["protocol"] = out_name.split("_", 1)[1]
        frame["subset"] = "all_ok_3"
        frame["number_of_models"] = len(used)
        frame["number_with_residue_resolved"] = counts
        frame["coverage_fraction"] = counts / max(len(used), 1)
        frame[["ensemble_mean_x_A", "ensemble_mean_y_A", "ensemble_mean_z_A"]] = mean
        frame["ensemble_rmsf_A"] = rmsf
        frame = add_reference_metrics(frame, mean, mean_sq, counts, refs, args.channel)
        rows.append(frame)
        summary[out_name] = {"allowlisted": len(allowed_by[out_name]), "aligned": len(used),
                             "failed": len(failed), "failures": failed[:20]}
    for name in cfg["vanilla"]:
        frame = template.loc[template.dataset.eq(name)].sort_values("raw_residue_number").reset_index(drop=True)
        old = legacy.loc[legacy.dataset.eq(name)].set_index("raw_residue_number")
        for column in [c for c in frame.columns if c.startswith(METRIC_PREFIXES) and c in old]:
            frame[column] = frame.raw_residue_number.map(old[column]).to_numpy()
        frame["subset"] = "all_ok_3"
        rows.append(frame)
        summary[name] = {"source": "legacy final-QC vanilla profile", "models": int(old.number_of_models.iloc[0])}

    out_dir = root / args.channel / "dataRMSF" / "profiles"
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(out_dir / f"{args.channel}_all_ok_3_per_residue_profiles.csv", index=False)
    (out_dir / f"{args.channel}_all_ok_3_stream_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: {kk: v for kk, v in s.items() if kk != "failures"} for k, s in summary.items()}, indent=2))


if __name__ == "__main__":
    main()
