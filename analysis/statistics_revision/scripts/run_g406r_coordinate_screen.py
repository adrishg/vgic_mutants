#!/usr/bin/env python3
"""Calculate a reproducible R406-centered screen from model coordinates.

The primary requested rule uses every non-hydrogen atom in ARG406 and the
partner residue. A side-chain-heavy rule is emitted as a prespecified
sensitivity because the all-heavy rule includes the covalent ARG406-GLU407
peptide bond. Results from this script are new coordinate-derived analyses;
they are not presented as recovery of the historical ``shortest_ARG406-*``
implementation.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.statistics_revision.scripts.run_paired_seed_v2 import (  # noqa: E402
    R406_SCREEN_PARTNERS,
    paired_scalar_rows,
    trajectory_summary,
)
from analysis.statistics_revision.scripts.run_seed_block_revision import (  # noqa: E402
    FINAL_DISTANCE_PATHS,
    add_metadata,
    read_csv,
)


HF_REPO = "adrishgz/vgic-mutant-structural-ensembles"
HF_PREFIX = "data/Cav1.2/G406R"
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})
RULES = {
    "all_non_hydrogen_atoms": lambda name: True,
    "side_chain_non_hydrogen_atoms": lambda name: name not in BACKBONE_ATOMS,
}


def coordinate_folder(protocol: str) -> str:
    return "af2_masked_models" if protocol == "masked" else "af2_vanilla_models"


def _element(line: str, atom_name: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    return atom_name.lstrip("0123456789")[:1].upper()


def selected_atoms(path: Path) -> dict[int, list[tuple[str, np.ndarray]]]:
    """Read non-hydrogen atoms for ARG406 and the 19 partner residues."""
    residue_numbers = {406, *(int(partner[3:]) for partner in R406_SCREEN_PARTNERS)}
    expected_names = {406: "ARG", **{
        int(partner[3:]): partner[:3] for partner in R406_SCREEN_PARTNERS
    }}
    atoms: dict[int, list[tuple[str, np.ndarray]]] = {
        number: [] for number in residue_numbers
    }
    observed_names: dict[int, set[str]] = {number: set() for number in residue_numbers}
    observed_chains: dict[int, set[str]] = {number: set() for number in residue_numbers}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                residue_number = int(line[22:26])
            except ValueError:
                continue
            if residue_number not in residue_numbers:
                continue
            altloc = line[16]
            if altloc not in {" ", "A"}:
                continue
            atom_name = line[12:16].strip()
            if _element(line, atom_name) in {"H", "D"}:
                continue
            residue_name = line[17:20].strip().upper()
            observed_names[residue_number].add(residue_name)
            observed_chains[residue_number].add(line[21].strip() or "_")
            atoms[residue_number].append((
                atom_name,
                np.array([
                    float(line[30:38]), float(line[38:46]), float(line[46:54])
                ]),
            ))
    for residue_number, expected_name in expected_names.items():
        if observed_names[residue_number] != {expected_name}:
            raise ValueError(
                f"{path.name}: residue {residue_number} expected {expected_name}; "
                f"observed {sorted(observed_names[residue_number])}"
            )
        if len(observed_chains[residue_number]) != 1:
            raise ValueError(
                f"{path.name}: residue {residue_number} appears in chains "
                f"{sorted(observed_chains[residue_number])}"
            )
    return atoms


def minimum_pair(
    focal_atoms: list[tuple[str, np.ndarray]],
    partner_atoms: list[tuple[str, np.ndarray]],
    rule: str,
) -> tuple[float, str, str]:
    selector = RULES[rule]
    focal = [(name, xyz) for name, xyz in focal_atoms if selector(name)]
    partner = [(name, xyz) for name, xyz in partner_atoms if selector(name)]
    if not focal or not partner:
        raise ValueError(f"No atoms remain under {rule}")
    focal_xyz = np.stack([xyz for _, xyz in focal])
    partner_xyz = np.stack([xyz for _, xyz in partner])
    distances = np.linalg.norm(
        focal_xyz[:, np.newaxis, :] - partner_xyz[np.newaxis, :, :], axis=2
    )
    focal_index, partner_index = np.unravel_index(np.argmin(distances), distances.shape)
    return (
        float(distances[focal_index, partner_index]),
        focal[focal_index][0],
        partner[partner_index][0],
    )


def analyze_snapshot(
    protocol: str,
    metadata: dict[str, object],
    coordinate_path: Path,
) -> list[dict[str, object]]:
    atoms = selected_atoms(coordinate_path)
    rows = []
    for rule in RULES:
        rule_rows = []
        for partner in R406_SCREEN_PARTNERS:
            partner_number = int(partner[3:])
            distance, focal_atom, partner_atom = minimum_pair(
                atoms[406], atoms[partner_number], rule
            )
            rule_rows.append({
                "protocol": protocol,
                "pdb_file": metadata["pdb_file"],
                "seed": metadata["seed"],
                "model_number": metadata["model_number"],
                "recycle_number": metadata["recycle_number"],
                "atom_selection_rule": rule,
                "focal_residue": "ARG406",
                "partner_residue": partner,
                "minimum_distance_A": distance,
                "r406_atom_name": focal_atom,
                "partner_atom_name": partner_atom,
                "pair_pass_at_2A": distance >= 2.0,
                "pair_contact_at_4A": distance <= 4.0,
                "coordinate_dataset_path": (
                    f"{HF_PREFIX}/{coordinate_folder(protocol)}/"
                    f"{metadata['pdb_file']}"
                ),
            })
        screen_pass = all(row["pair_pass_at_2A"] for row in rule_rows)
        for row in rule_rows:
            row["screen_pass_at_2A"] = screen_pass
        rows.extend(rule_rows)
    return rows


def resolve_coordinate(
    protocol: str,
    pdb_file: str,
    coordinate_root: Path | None,
    *,
    download: bool,
) -> Path:
    folder = coordinate_folder(protocol)
    relative = Path(HF_PREFIX) / folder / pdb_file
    if coordinate_root is not None:
        candidates = [
            coordinate_root / relative,
            coordinate_root / folder / pdb_file,
            coordinate_root / pdb_file,
        ]
        for candidate in candidates:
            if candidate.is_file():
                return candidate
    if not download:
        raise FileNotFoundError(
            f"Coordinate unavailable for {pdb_file}; pass --download or --coordinate-root"
        )
    from huggingface_hub import hf_hub_download

    last_error: Exception | None = None
    for attempt in range(1, 6):
        try:
            return Path(hf_hub_download(
                HF_REPO, repo_type="dataset", filename=str(relative)
            ))
        except Exception as error:  # network retries are bounded and visible
            last_error = error
            if attempt == 5:
                break
            wait_seconds = 15 * attempt
            print(
                f"retrying {pdb_file} after {wait_seconds}s "
                f"({type(error).__name__}, attempt {attempt}/5)",
                flush=True,
            )
            time.sleep(wait_seconds)
    assert last_error is not None
    raise last_error


def summarize(records: pd.DataFrame) -> pd.DataFrame:
    snapshot = records.drop_duplicates(
        ["protocol", "pdb_file", "atom_selection_rule"]
    )[[
        "protocol", "pdb_file", "seed", "model_number", "atom_selection_rule",
        "screen_pass_at_2A",
    ]]
    rows = []
    for (rule, protocol), group in snapshot.groupby(
        ["atom_selection_rule", "protocol"], sort=False
    ):
        passed_names = set(group.loc[group.screen_pass_at_2A, "pdb_file"])
        trajectory_any = group.groupby(
            ["seed", "model_number"], as_index=False
        ).screen_pass_at_2A.max()
        base = {
            "atom_selection_rule": rule,
            "protocol": protocol,
            "total_snapshots": len(group),
            "passing_snapshots": int(group.screen_pass_at_2A.sum()),
            "unconditional_pass_fraction": float(group.screen_pass_at_2A.mean()),
            "total_trajectories": len(trajectory_any),
            "trajectories_with_any_pass": int(trajectory_any.screen_pass_at_2A.sum()),
            "trajectory_any_pass_fraction": float(trajectory_any.screen_pass_at_2A.mean()),
            "screen_definition": (
                "all 19 ARG406-partner minimum distances >= 2.0 A"
            ),
            "historical_recovery_status": (
                "new coordinate-derived analysis; not recovery of historical columns"
            ),
        }
        rows.append({**base, "partner_residue": "screen_only"})
        for partner in ("ASP1528", "ASP1533"):
            pair = records[
                records.atom_selection_rule.eq(rule)
                & records.protocol.eq(protocol)
                & records.partner_residue.eq(partner)
                & records.pdb_file.isin(passed_names)
            ]
            numerator = int(pair.pair_contact_at_4A.sum())
            denominator = len(pair)
            rows.append({
                **base,
                "partner_residue": partner,
                "conditional_contact_numerator": numerator,
                "conditional_contact_denominator": denominator,
                "conditional_contact_fraction": (
                    numerator / denominator if denominator else np.nan
                ),
                "contact_definition": "minimum distance <= 4.0 A among screen passes",
            })
    return pd.DataFrame(rows)


def hierarchical_statistics(
    records: pd.DataFrame, *, n_bootstrap: int = 2000
) -> pd.DataFrame:
    """Whole-seed contrasts for each declared coordinate atom rule."""
    rows = []
    for rule_index, rule in enumerate(RULES):
        rule_records = records[records.atom_selection_rule.eq(rule)]
        snapshots = rule_records.drop_duplicates(["protocol", "pdb_file"])[[
            "protocol", "pdb_file", "seed", "model_number", "recycle_number",
            "screen_pass_at_2A",
        ]].copy()
        conditions = {
            protocol: trajectory_summary(
                snapshots[snapshots.protocol.eq(protocol)],
                "screen_pass_at_2A", "mean",
            )
            for protocol in ("vanilla", "masked")
        }
        result = paired_scalar_rows(
            conditions["vanilla"], conditions["masked"], "screen_pass_at_2A",
            contrast="G406R masked - vanilla",
            outcome="unconditional_snapshot_screen_pass_fraction",
            bootstrap=n_bootstrap,
            seed=20269000 + rule_index * 100,
        )
        for row in result:
            row["atom_selection_rule"] = rule
        rows.extend(result)

        any_pass = {
            protocol: snapshots[snapshots.protocol.eq(protocol)].groupby(
                ["seed", "model_number"], as_index=False
            ).screen_pass_at_2A.max()
            for protocol in ("vanilla", "masked")
        }
        result = paired_scalar_rows(
            any_pass["vanilla"], any_pass["masked"], "screen_pass_at_2A",
            contrast="G406R masked - vanilla",
            outcome="trajectory_any_screen_pass_fraction",
            bootstrap=n_bootstrap,
            seed=20269020 + rule_index * 100,
        )
        for row in result:
            row["atom_selection_rule"] = rule
        rows.extend(result)

        for partner_index, partner in enumerate(("ASP1528", "ASP1533")):
            pair = rule_records[rule_records.partner_residue.eq(partner)].copy()
            pair["conditional_contact"] = np.where(
                pair.screen_pass_at_2A,
                pair.pair_contact_at_4A.astype(float),
                np.nan,
            )
            conditions = {
                protocol: trajectory_summary(
                    pair[pair.protocol.eq(protocol)], "conditional_contact", "mean"
                )
                for protocol in ("vanilla", "masked")
            }
            if any(condition.empty for condition in conditions.values()):
                rows.append({
                    "contrast": "G406R masked - vanilla",
                    "outcome": f"conditional_{partner}_contact_fraction",
                    "estimand_id": "not_estimable",
                    "estimate_A": np.nan,
                    "estimate_B": np.nan,
                    "estimate_B_minus_A": np.nan,
                    "CI_low_B_minus_A": np.nan,
                    "CI_high_B_minus_A": np.nan,
                    "bootstrap_replicates": n_bootstrap,
                    "atom_selection_rule": rule,
                    "contact_rule": (
                        "minimum distance <= 4.0 A among screen passes"
                    ),
                    "estimability_note": (
                        "not estimable because at least one protocol has no screen passes"
                    ),
                })
                continue
            result = paired_scalar_rows(
                conditions["vanilla"], conditions["masked"], "conditional_contact",
                contrast="G406R masked - vanilla",
                outcome=f"conditional_{partner}_contact_fraction",
                bootstrap=n_bootstrap,
                seed=20269040 + rule_index * 100 + partner_index * 10,
            )
            for row in result:
                row["atom_selection_rule"] = rule
                row["contact_rule"] = "minimum distance <= 4.0 A among screen passes"
            rows.extend(result)
    table = pd.DataFrame(rows)
    table["resampling_unit"] = "whole recorded nominal seed label"
    table["actual_random_seed_pairing_status"] = (
        "matching numeric labels are nominal design keys; common RNG not verified"
    )
    return table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinate-root", type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=Path("/tmp/vgci_g406r_coordinate_screen.csv"),
        help="Temporary resumable record table kept outside the repository",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "analysis/statistics_revision/paired_seed_v2",
    )
    args = parser.parse_args()

    jobs = []
    for protocol in ("vanilla", "masked"):
        cohort = add_metadata(read_csv(
            ROOT / FINAL_DISTANCE_PATHS[f"cav12|G406R|{protocol}"]
        ))
        for metadata in cohort[
            ["pdb_file", "seed", "model_number", "recycle_number"]
        ].to_dict("records"):
            jobs.append((protocol, metadata))

    def run_job(job: tuple[str, dict[str, object]]) -> list[dict[str, object]]:
        protocol, metadata = job
        coordinate = resolve_coordinate(
            protocol, str(metadata["pdb_file"]), args.coordinate_root,
            download=args.download,
        )
        return analyze_snapshot(protocol, metadata, coordinate)

    completed: set[tuple[str, str]] = set()
    if args.checkpoint.is_file():
        checkpoint = pd.read_csv(args.checkpoint, usecols=["protocol", "pdb_file"])
        counts = checkpoint.groupby(["protocol", "pdb_file"]).size()
        completed = set(counts[counts.eq(len(RULES) * len(R406_SCREEN_PARTNERS))].index)
        print(f"resuming after {len(completed)} complete snapshots", flush=True)
    pending = [
        job for job in jobs
        if (job[0], str(job[1]["pdb_file"])) not in completed
    ]
    processed = len(completed)
    for start in range(0, len(pending), args.chunk_size):
        chunk = pending[start:start + args.chunk_size]
        records: list[dict[str, object]] = []
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            for rows in executor.map(run_job, chunk):
                records.extend(rows)
        pd.DataFrame(records).to_csv(
            args.checkpoint,
            index=False,
            mode="a" if args.checkpoint.exists() else "w",
            header=not args.checkpoint.exists(),
        )
        processed += len(chunk)
        print(f"processed {processed}/{len(jobs)} snapshots", flush=True)

    frame = pd.read_csv(args.checkpoint)
    key_columns = ["protocol", "pdb_file", "atom_selection_rule", "partner_residue"]
    if frame.duplicated(key_columns).any():
        raise ValueError("Checkpoint contains duplicate snapshot-rule-partner rows")
    expected_rows = len(jobs) * len(RULES) * len(R406_SCREEN_PARTNERS)
    if len(frame) != expected_rows:
        raise ValueError(f"Expected {expected_rows} coordinate rows; observed {len(frame)}")
    source = args.output_dir / "source_data"
    source.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        source / "g406r_coordinate_heavy_atom_pair_distances.csv.gz",
        index=False,
    )
    summarize(frame).to_csv(
        args.output_dir / "g406r_coordinate_screen_summary.csv", index=False
    )
    hierarchical_statistics(frame).to_csv(
        args.output_dir / "g406r_coordinate_screen_statistics.csv", index=False
    )


if __name__ == "__main__":
    main()
