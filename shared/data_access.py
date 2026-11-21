"""Resolution of data paths against the working tree, a local cache, or Hugging Face.

Large derived tables are published in the companion dataset rather than tracked
in Git. Small provenance files stay in the repository and are returned as they
are. ``resolve`` always hands back a real local path, so it serves ``np.load``
and ``open`` as readily as ``pandas.read_csv``.

Set ``VGIC_DATA_ROOT`` to a complete existing copy to suppress fetching, and
``VGIC_OFFLINE`` to fail rather than reach the network.
"""

from __future__ import annotations

import os
from pathlib import Path

HF_REPO = "adrishgz/vgic-mutant-structural-ensembles"
HF_REVISION = os.environ.get("VGIC_HF_REVISION", "main")

REPO_ROOT = Path(__file__).resolve().parents[1]
CACHE_ROOT = Path(os.environ.get("VGIC_CACHE", Path.home() / ".cache" / "vgic-mutants"))

_EXTRA_ROOT = os.environ.get("VGIC_DATA_ROOT") or os.environ.get("VGCI_RECONCILIATION_DATA_ROOT")
if _EXTRA_ROOT:
    EXTRA_ROOTS = [Path(_EXTRA_ROOT).expanduser()]
else:
    EXTRA_ROOTS = []

_CHANNELS = {"cav12": "Cav1.2", "kv21": "Kv2.1", "nav15": "Nav1.5"}

# External reconciliation data retain the channel names used by the companion
# release and may be nested below an archive prefix. Repository notebooks use
# the short channel directory names so they remain portable.
_EXTERNAL_CHANNEL_NAMES = {
    "cav12": "Cav1.2",
    "kv21": "Kv2.1",
    "nav15": "Nav1.5",
}

# The dataset already uses ``distances`` for the pre-QC tables, so that name is kept.
_BUCKETS = {
    "dataDistances": "distances",
    "dataRMSD": "rmsd",
    "dataRMSF": "rmsf",
    "dataExtra": "extra",
    "experimental": "experimental",
    "rmsd_threshold_sensitivity": "sensitivity",
    "rmsd_filtered_distances": "archive/rmsd_filtered_distances",
    "rmsd_convergence_filtering": "convergence",
}

_KEPT_TOP_LEVEL = ("docs", "notebooks", "scripts", "tests", "shared", "paperFigures")
_KEPT_SECOND = ("rmsd_convergence_filtering", "experimental")
_KEPT_PREFIXES = ("analysis/statistics_revision/",)


class DataUnavailable(RuntimeError):
    """A required file is neither present locally nor retrievable."""


def is_kept_in_repo(rel_path: str | Path) -> bool:
    """True when this path is version-controlled rather than published."""
    parts = Path(rel_path).parts
    if not parts:
        return False
    if parts[0] in _KEPT_TOP_LEVEL:
        return True
    if str(Path(rel_path)).startswith(_KEPT_PREFIXES):
        return True
    return len(parts) > 1 and parts[1] in _KEPT_SECOND


# The published dataset does not spell every channel the way this repository
# does, and the CaV1.2 buckets carry a cohort level the others do not.  Rather
# than pin one spelling, every known layout is tried in turn (D-081).
_CHANNEL_ALIASES = {
    "cav12": ("CaV1.2", "Cav1.2"),
    "kv21": ("Kv2.1",),
    "nav15": ("Nav1.5", "NaV1.5"),
}

_BUCKET_SUBDIRS = ("", "50pct")


def hf_candidates(rel_path: str | Path) -> list[str]:
    """Every dataset path a repository-relative path might be published under.

    The first entry is the canonical layout; the rest are the spellings and
    cohort subdirectories the dataset actually uses.  Callers try them in
    order, so adding one here can never break a path that already resolves.
    """
    parts = Path(rel_path).parts
    if len(parts) < 2:
        raise ValueError(f"Cannot map {rel_path!r} to a dataset path")

    if parts[0] not in _CHANNELS:
        return [str(Path(*parts))]

    bucket = _BUCKETS.get(parts[1])
    if bucket is None:
        raise ValueError(f"No dataset bucket defined for {parts[0]}/{parts[1]}")

    names = (_CHANNELS[parts[0]], *_CHANNEL_ALIASES.get(parts[0], ()))
    remainder = Path(*parts[2:])
    candidates = []
    for channel in dict.fromkeys(names):
        for subdir in _BUCKET_SUBDIRS:
            base = Path("data") / channel / bucket
            candidate = str((base / subdir / remainder) if subdir else (base / remainder))
            if candidate not in candidates:
                candidates.append(candidate)
    return candidates


def hf_path(rel_path: str | Path) -> str:
    """Map a repository-relative path to its canonical path inside the dataset."""
    return hf_candidates(rel_path)[0]


def _is_lfs_pointer(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(43).startswith(b"version https://git-lfs.github.com/spec/v1")
    except OSError:
        return False


def _usable(path: Path) -> bool:
    return path.is_file() and not _is_lfs_pointer(path)


def _candidate_paths(root: Path, rel: Path):
    """Yield supported layouts for one repository-relative data path."""
    yield root / rel
    parts = rel.parts
    if not parts:
        return
    channel_name = _EXTERNAL_CHANNEL_NAMES.get(parts[0].lower())
    if channel_name is None:
        return
    remainder = Path(*parts[1:])
    yield root / channel_name / remainder
    yield root / "hive" / channel_name / remainder


def resolve(rel_path: str | Path) -> Path:
    """Return a readable local path, fetching it once if it is not already here."""
    rel = Path(rel_path)
    if rel.is_absolute():
        try:
            rel = rel.relative_to(REPO_ROOT)
        except ValueError:
            # Temporary/user-supplied files outside the repository are not
            # companion-dataset keys. Read them directly or report absence.
            if _usable(rel):
                return rel
            raise FileNotFoundError(rel)

    for root in (REPO_ROOT, *EXTRA_ROOTS):
        for candidate in _candidate_paths(root, rel):
            if _usable(candidate):
                return candidate

    cached = CACHE_ROOT / rel
    if _usable(cached):
        return cached

    if os.environ.get("VGIC_OFFLINE"):
        raise DataUnavailable(
            f"{rel} is not available locally and VGIC_OFFLINE is set. "
            f"Clear VGIC_OFFLINE, or point VGIC_DATA_ROOT at a full copy."
        )
    return _fetch(rel, cached)


def _fetch(rel: Path, destination: Path) -> Path:
    from huggingface_hub import hf_hub_download

    candidates = hf_candidates(rel)
    downloaded = None
    last_error: Exception | None = None
    for candidate in candidates:
        try:
            downloaded = hf_hub_download(
                repo_id=HF_REPO,
                repo_type="dataset",
                revision=HF_REVISION,
                filename=candidate,
            )
            break
        except Exception as error:  # noqa: BLE001 - reported with every path tried
            last_error = error
    if downloaded is None:
        tried = "\n  ".join(candidates)
        raise DataUnavailable(
            f"{rel} is not in the working tree and could not be fetched from "
            f"{HF_REPO}. Tried:\n  {tried}\nLast error: {last_error}"
        ) from last_error

    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        try:
            destination.symlink_to(Path(downloaded).resolve())
        except OSError:
            import shutil

            shutil.copy2(downloaded, destination)
    return destination


def read_csv(rel_path: str | Path, **kwargs):
    import pandas as pd

    return pd.read_csv(resolve(rel_path), **kwargs)


def load_npy(rel_path: str | Path, **kwargs):
    """NumPy cannot read an array over HTTP, so the file is materialised first."""
    import numpy as np

    return np.load(resolve(rel_path), **kwargs)
