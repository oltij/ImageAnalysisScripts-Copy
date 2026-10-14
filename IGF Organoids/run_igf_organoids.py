#!/usr/bin/env python3
"""Launch the validated IGF organoid profile through the shared pipeline."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_RUNNER = REPOSITORY_ROOT / "run_stitched_pipeline.py"
IGF_ALIGNMENT_ARGS = (
    "--skip-phase-init",
    "--match-gate-um",
    "25",
    "--iterations",
    "5",
    "--max-points",
    "1500",
)


def validate_igf_alignment_profile(config_path: Path) -> None:
    """Refuse an IGF run that does not use the validated registration profile."""
    raw = json.loads(config_path.read_text(encoding="utf-8"))
    alignment = raw.get("alignment")
    if not isinstance(alignment, dict) or alignment.get("enabled") is not True:
        raise ValueError("IGF configuration must set alignment.enabled to true")
    actual = alignment.get("extra_args")
    if actual != list(IGF_ALIGNMENT_ARGS):
        raise ValueError(
            "IGF alignment.extra_args must be exactly "
            f"{list(IGF_ALIGNMENT_ARGS)}; affine or alternate registration "
            "settings are not part of the validated IGF workflow"
        )


def build_command(
    config_path: Path, *, dry_run: bool = False, resume: bool = False
) -> list[str]:
    """Build the same-interpreter command for the repository pipeline runner."""
    command = [
        sys.executable,
        str(PIPELINE_RUNNER),
        "--config",
        str(config_path.expanduser().resolve()),
    ]
    if dry_run:
        command.append("--dry-run")
    if resume:
        command.append("--resume")
    return command


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)

    config_path = args.config.expanduser().resolve()
    if not PIPELINE_RUNNER.is_file():
        raise FileNotFoundError(f"Shared pipeline runner not found: {PIPELINE_RUNNER}")
    validate_igf_alignment_profile(config_path)
    command = build_command(
        config_path, dry_run=args.dry_run, resume=args.resume
    )
    os.execv(sys.executable, command)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
