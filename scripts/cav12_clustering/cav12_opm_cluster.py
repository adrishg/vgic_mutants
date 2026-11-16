#!/usr/bin/env python3
"""OPM-framed exploratory clustering for CaV1.2 structural ensembles.

The OPM structure defines the common membrane coordinate frame.  Alignment is
fit only on the OPM-mapped S1--S3 segments, which are used as a stable scaffold;
the mobile S4--S6/pore segments are retained as clustering features.  Recycle
snapshots should not be supplied as independent observations: use one final or
converged model per trajectory.

Raw structures and generated tables are intentionally external inputs/outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import deque
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from scipy.spatial import ConvexHull, cKDTree, distance


# Raw AlphaFold numbering.  These are the curated CaV1.2 topology ranges used
# by the existing RMSF package; OPM supplies the membrane-frame test below.
STABLE_S1_S3 = [
    (125, 143), (159, 179), (189, 209),
    (525, 543), (555, 575), (587, 606),
    (901, 919), (932, 952), (988, 1006),
    (1240, 1261), (1270, 1291), (1302, 1321),
]
MOBILE_S4_S6 = [
    (233, 251), (269, 290), (381, 412),
    (616, 634), (654, 673), (726, 745),
    (1014, 1032), (1052, 1071), (1160, 1181),
    (1373, 1391), (1410, 1430), (1500, 1524),
]
FILTER_RESIDUES = [363, 706, 1135, 1464]
GATE_RESIDUES = [400, 749, 1186, 1519]
OPPOSITE_GATE_PAIRS = [(400, 1186), (749, 1519)]


def ranges_to_numbers(ranges: list[tuple[int, int]]) -> list[int]:
    return [n for start, end in ranges for n in range(start, end + 1)]


def aa(residue) -> str:
    return protein_letters_3to1.get(residue.resname.upper(), "X")


def polymer_residues(chain):
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def load_ca(path: Path, chain_id: str = "A") -> tuple[dict[int, np.ndarray], str]:
    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    model = structure[0]
    chain = model[chain_id] if chain_id in model else next(c for c in model if polymer_residues(c))
    residues = polymer_residues(chain)
    coords = {int(r.id[1]): np.asarray(r["CA"].coord, dtype=float) for r in residues}
    sequence = "".join(aa(r) for r in residues)
    return coords, sequence


def sequence_mapping(model_sequence: str, opm_sequence: str, model_numbers: list[int], opm_numbers: list[int]) -> dict[int, int]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(opm_sequence, model_sequence)[0]
    mapping: dict[int, int] = {}
    for (opm_start, opm_end), (model_start, model_end) in zip(*alignment.aligned):
        for oi, mi in zip(range(opm_start, opm_end), range(model_start, model_end)):
            if opm_sequence[oi] != "X" and model_sequence[mi] != "X":
                mapping[model_numbers[mi]] = opm_numbers[oi]
    return mapping


def kabsch(moving: np.ndarray, fixed: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    moving_center = moving.mean(axis=0)
    fixed_center = fixed.mean(axis=0)
    x = moving - moving_center
    y = fixed - fixed_center
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    transformed = x @ rotation + fixed_center
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - fixed) ** 2, axis=1))))
    return rotation, fixed_center - moving_center @ rotation, rmsd


def transform(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coords @ rotation + translation


def gate_cross_section_area(coords: dict[int, np.ndarray]) -> float:
    """Return the convex-hull area of the four gate Cα positions in x-y."""
    points = np.asarray([[coords[n][0], coords[n][1]] for n in GATE_RESIDUES], dtype=float)
    if len(points) < 3:
        return math.nan
    return float(ConvexHull(points).volume)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "model_path" not in reader.fieldnames:
            raise ValueError("Manifest must be TSV with a model_path column")
        return list(reader)


def robust_scale(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = np.nanmedian(x, axis=0)
    mad = np.nanmedian(np.abs(x - center), axis=0) * 1.4826
    mad[mad < 1e-8] = 1.0
    return (x - center) / mad, center, mad


def pca(x: np.ndarray, variance_target: float = 0.90, max_components: int = 30):
    scaled, center, scale = robust_scale(x)
    u, singular, vt = np.linalg.svd(scaled, full_matrices=False)
    variance = singular**2
    explained = variance / variance.sum() if variance.sum() else variance
    cumulative = np.cumsum(explained)
    n = int(np.searchsorted(cumulative, variance_target) + 1)
    n = max(2, min(n, max_components, u.shape[1]))
    scores = u[:, :n] * singular[:n]
    return scores, center, scale, explained[:n], n


def dbscan(x: np.ndarray, eps: float, min_samples: int) -> np.ndarray:
    neighbors = cKDTree(x).query_ball_point(x, r=eps)
    labels = np.full(len(x), -1, dtype=int)
    visited = np.zeros(len(x), dtype=bool)
    cluster_id = 0
    for start in range(len(x)):
        if visited[start]:
            continue
        visited[start] = True
        seed_neighbors = neighbors[start]
        if len(seed_neighbors) < min_samples:
            continue
        labels[start] = cluster_id
        queue = deque(seed_neighbors)
        queued = set(seed_neighbors)
        while queue:
            point = queue.popleft()
            if not visited[point]:
                visited[point] = True
                point_neighbors = neighbors[point]
                if len(point_neighbors) >= min_samples:
                    for neighbor in point_neighbors:
                        if neighbor not in queued:
                            queue.append(neighbor)
                            queued.add(neighbor)
            if labels[point] == -1:
                labels[point] = cluster_id
        cluster_id += 1
    return labels


def silhouette(x: np.ndarray, labels: np.ndarray) -> float:
    keep = labels >= 0
    if keep.sum() < 3 or len(set(labels[keep])) < 2:
        return math.nan
    data = x[keep]
    labs = labels[keep]
    matrix = distance.cdist(data, data)
    scores = []
    for i, label in enumerate(labs):
        own = labs == label
        own[i] = False
        a = matrix[i, own].mean() if own.any() else 0.0
        b = min(matrix[i, labs == other].mean() for other in set(labs) if other != label)
        scores.append((b - a) / max(a, b) if max(a, b) else 0.0)
    return float(np.mean(scores))


def cluster_panel(name: str, x: np.ndarray, model_rows: list[dict[str, str]], out: Path, min_samples: int):
    scores, center, scale, explained, n_components = pca(x)
    tree = cKDTree(scores)
    k = min(min_samples + 1, len(scores))
    kth = tree.query(scores, k=k)[0][:, -1]
    quantiles = {q: float(np.quantile(kth, q)) for q in (0.80, 0.85, 0.90, 0.95)}
    candidate_rows = []
    candidate_labels = {}
    for q, eps in quantiles.items():
        labels = dbscan(scores, eps, min_samples)
        candidate_labels[q] = labels
        nonnoise = labels[labels >= 0]
        sizes = sorted(np.bincount(nonnoise).tolist(), reverse=True) if len(nonnoise) else []
        candidate_rows.append({
            "panel": name,
            "eps_quantile": q,
            "eps": eps,
            "n_clusters": len(set(nonnoise)),
            "n_noise": int((labels < 0).sum()),
            "noise_fraction": float((labels < 0).mean()),
            "min_cluster_size": min(sizes) if sizes else 0,
            "max_cluster_size": max(sizes) if sizes else 0,
            "silhouette_nonnoise": silhouette(scores, labels),
        })
    primary_q = 0.90 if candidate_rows[2]["n_clusters"] >= 2 else max(
        candidate_rows, key=lambda r: (-1 if math.isnan(r["silhouette_nonnoise"]) else r["silhouette_nonnoise"])
    )["eps_quantile"]
    labels = candidate_labels[primary_q]

    with (out / f"{name}_pca_scores.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["model_path", "cluster", "PC1", "PC2", *[f"PC{i+1}" for i in range(2, n_components)]])
        for row, label, score in zip(model_rows, labels, scores):
            writer.writerow([row["model_path"], int(label), *score.tolist()])
    with (out / f"{name}_cluster_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(candidate_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)
    with (out / f"{name}_cluster_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["panel", "primary_eps_quantile", "primary_eps", "cluster", "n_models", "fraction", "medoid_model"])
        for cluster in sorted(c for c in set(labels) if c >= 0):
            indices = np.flatnonzero(labels == cluster)
            centroid = scores[indices].mean(axis=0)
            medoid = indices[np.argmin(np.sum((scores[indices] - centroid) ** 2, axis=1))]
            writer.writerow([name, primary_q, quantiles[primary_q], cluster, len(indices), len(indices) / len(labels), model_rows[medoid]["model_path"]])
        writer.writerow([name, primary_q, quantiles[primary_q], -1, int((labels < 0).sum()), float((labels < 0).mean()), ""])
    np.savez_compressed(out / f"{name}_pca_model.npz", center=center, scale=scale, explained=explained, primary_eps=quantiles[primary_q])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 5.6))
        unique = sorted(set(labels))
        for cluster in unique:
            keep = labels == cluster
            label = "noise" if cluster < 0 else f"cluster {cluster}"
            ax.scatter(scores[keep, 0], scores[keep, 1], s=16, alpha=0.75, label=label)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title(f"CaV1.2 OPM-framed clustering: {name}")
        if len(unique) <= 12:
            ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(out / f"{name}_pca_clusters.png", dpi=220)
        plt.close(fig)
    except Exception as exc:
        (out / f"{name}_plot_warning.txt").write_text(f"Plot skipped: {type(exc).__name__}: {exc}\n")
    return labels, scores, candidate_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--opm-pdb", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--opm-chain", default="K")
    parser.add_argument("--model-chain", default="A")
    parser.add_argument("--opm-half-thickness", type=float, default=14.5)
    parser.add_argument("--anchor-margin", type=float, default=10.0)
    parser.add_argument("--min-samples", type=int, default=10)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    opm_coords, opm_sequence = load_ca(args.opm_pdb, args.opm_chain)
    opm_numbers = sorted(opm_coords)
    rows = read_manifest(args.manifest)
    if not rows:
        raise ValueError("Manifest is empty")
    first_coords, first_sequence = load_ca(Path(rows[0]["model_path"]), args.model_chain)
    model_numbers = sorted(first_coords)
    mapping = sequence_mapping(first_sequence, opm_sequence, model_numbers, opm_numbers)

    stable_numbers = ranges_to_numbers(STABLE_S1_S3)
    stable_anchors = [n for n in stable_numbers if n in mapping and mapping[n] in opm_coords and abs(opm_coords[mapping[n]][2]) <= args.opm_half_thickness + args.anchor_margin]
    mobile_numbers = ranges_to_numbers(MOBILE_S4_S6)
    feature_numbers = [n for n in mobile_numbers if n in mapping]
    if len(stable_anchors) < 80:
        raise ValueError(f"Only {len(stable_anchors)} OPM-framed stable anchors available")

    aligned_mobile = []
    aligned_geometry = []
    gate_metric_rows = []
    alignment_rows = []
    key_numbers = FILTER_RESIDUES + GATE_RESIDUES
    for row in rows:
        path = Path(row["model_path"])
        coords, sequence = load_ca(path, args.model_chain)
        ref = np.array([opm_coords[mapping[n]] for n in stable_anchors if n in coords], dtype=float)
        moving = np.array([coords[n] for n in stable_anchors if n in coords], dtype=float)
        if len(ref) < 80:
            raise ValueError(f"{path}: insufficient stable anchors ({len(ref)})")
        rotation, translation, fit_rmsd = kabsch(moving, ref)
        mobile = np.array([transform(coords[n], rotation, translation) for n in feature_numbers], dtype=float)
        key = {n: transform(coords[n], rotation, translation) for n in key_numbers if n in coords}
        geometry = []
        for group in (FILTER_RESIDUES, GATE_RESIDUES):
            for i, left in enumerate(group):
                for right in group[i + 1:]:
                    geometry.append(float(np.linalg.norm(key[left] - key[right])))
        for n in FILTER_RESIDUES:
            geometry.append(float(np.linalg.norm(key[n] - key[GATE_RESIDUES[FILTER_RESIDUES.index(n)]])))
        opposite_distances = {
            f"gate_opposite_pair_{left}_{right}_ca_distance_A": float(np.linalg.norm(key[left] - key[right]))
            for left, right in OPPOSITE_GATE_PAIRS
        }
        gate_metric_rows.append({
            "model_path": str(path),
            "dataset": row.get("dataset", row.get("cohort", "")),
            **opposite_distances,
            "gate_opposite_min_ca_distance_A": min(opposite_distances.values()),
            "gate_opposite_mean_ca_distance_A": float(np.mean(list(opposite_distances.values()))),
            "gate_cross_section_hull_area_A2": gate_cross_section_area(key),
        })
        aligned_mobile.append(mobile.reshape(-1))
        aligned_geometry.append(geometry)
        alignment_rows.append({
            "model_path": str(path),
            "dataset": row.get("dataset", row.get("cohort", "")),
            "stable_anchor_count": len(ref),
            "stable_anchor_rmsd_A": fit_rmsd,
        })

    mobile = np.asarray(aligned_mobile, dtype=float)
    geometry = np.asarray(aligned_geometry, dtype=float)
    with (args.output_dir / "opm_core_definition.json").open("w") as handle:
        json.dump({
            "opm_pdb": str(args.opm_pdb),
            "opm_chain": args.opm_chain,
            "opm_half_thickness_A": args.opm_half_thickness,
            "anchor_margin_A": args.anchor_margin,
            "stable_s1_s3_ranges_raw": STABLE_S1_S3,
            "stable_anchor_raw_residues": stable_anchors,
            "stable_anchor_count": len(stable_anchors),
            "mobile_s4_s6_ranges_raw": MOBILE_S4_S6,
            "mobile_feature_residues": feature_numbers,
            "filter_residues": FILTER_RESIDUES,
            "gate_residues": GATE_RESIDUES,
            "note": "Primary alignment excludes S4-S6 and mutation-centered pore motion; these segments are retained for clustering.",
        }, handle, indent=2)
    with (args.output_dir / "opm_alignment_metadata.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(alignment_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(alignment_rows)
    with (args.output_dir / "intracellular_gate_metrics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gate_metric_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(gate_metric_rows)
    labels_mobile, _, _ = cluster_panel("mobile_s4_s6_coordinates", mobile, alignment_rows, args.output_dir, args.min_samples)
    labels_geometry, _, _ = cluster_panel("pore_geometry", geometry, alignment_rows, args.output_dir, args.min_samples)
    with (args.output_dir / "run_summary.json").open("w") as handle:
        json.dump({
            "n_models": len(rows),
            "stable_anchor_count": len(stable_anchors),
            "mobile_feature_count": len(feature_numbers),
            "mobile_primary_clusters": sorted(set(labels_mobile[labels_mobile >= 0].tolist())),
            "geometry_primary_clusters": sorted(set(labels_geometry[labels_geometry >= 0].tolist())),
            "intracellular_gate_metrics": {
                "minimum_opposite_gate_distance_quantiles_A": {
                    str(q): float(np.quantile([r["gate_opposite_min_ca_distance_A"] for r in gate_metric_rows], q))
                    for q in (0.05, 0.50, 0.95)
                },
                "cross_section_hull_area_quantiles_A2": {
                    str(q): float(np.quantile([r["gate_cross_section_hull_area_A2"] for r in gate_metric_rows], q))
                    for q in (0.05, 0.50, 0.95)
                },
                "metrics_tsv": "intracellular_gate_metrics.tsv",
            },
        }, handle, indent=2)


if __name__ == "__main__":
    main()
