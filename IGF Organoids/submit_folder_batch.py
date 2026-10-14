#!/usr/bin/env python3
"""One-command per-field Slurm submission for run_folder_batch.py.

Every selected Field_N becomes exactly one Slurm array task. Each task starts
with its own four IMS tiles and runs IMS export -> BaSiC -> BigStitcher ->
full downstream image analysis, using the existing scientific controllers.

No Slurm job is submitted during --print-only. Avoid duplicate submissions.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import run_folder_batch as batch


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parent
SBATCH_RUNNER = SCRIPT_DIR / 'run_folder_batch.sbatch'
BATCH_RUNNER = SCRIPT_DIR / 'run_folder_batch.py'


def submit(manifest_path: Path, max_parallel: int, *, print_only: bool = False,
           allow_duplicate: bool = False, sbatch: str = 'sbatch') -> str | None:
    manifest_path = batch.abs_path(manifest_path)
    manifest = batch.load_manifest(manifest_path)
    batch.validate_manifest_for_submission(manifest, check_sources=True)
    workspace = Path(manifest['workspace'])
    if manifest_path != workspace / 'manifest.json':
        batch.fail('Manifest must be located at <workspace>/manifest.json')
    jobs = manifest['jobs']
    n = len(jobs)
    if n < 1:
        batch.fail('Manifest has zero fields')
    if max_parallel < 1:
        batch.fail('--max-parallel must be 1 or greater')
    if [j['index'] for j in jobs] != list(range(n)):
        batch.fail('Manifest must have sequential field indices 0..N-1')
    runner = Path(manifest.get('slurm_runner', SBATCH_RUNNER))
    if not runner.is_file():
        batch.fail(f'Cannot find Slurm runner: {runner}')
    python_runner = Path(manifest.get('batch_runner', BATCH_RUNNER))
    if not python_runner.is_file():
        batch.fail(f'Cannot find Python batch runner: {python_runner}')
    recorded = workspace / 'slurm_submission.json'
    if recorded.exists() and not (print_only or allow_duplicate):
        previous = json.loads(recorded.read_text())
        batch.fail(
            f'Workspace already submitted as Slurm job {previous.get("slurm_job_id")}. '
            'Check squeue/sacct; only use --allow-duplicate after the previous tasks have ended.'
        )
    throttle = min(n, max_parallel)
    array = f'0-{n - 1}%{throttle}'
    cmd = [sbatch, '--parsable', f'--array={array}',
           f'--chdir={workspace}', str(runner), str(manifest_path),
           str(python_runner)]
    print(f'Selected {n} Field_N folders => {n} independent Slurm array tasks.')
    print(f'Max simultaneous field jobs: {throttle}')
    print('Per field: 4 IMS tiles -> BaSiC -> BigStitcher -> projections -> full downstream analysis')
    print('Submit command: ' + ' '.join(cmd), flush=True)
    if print_only:
        print('PRINT ONLY: no Slurm job submitted.')
        return None
    if not shutil.which(sbatch):
        batch.fail(f'Could not find Slurm command {sbatch!r}')
    result = subprocess.run(cmd, cwd=workspace, text=True, capture_output=True, check=True)
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)
    value = result.stdout.strip()
    match = re.fullmatch(r'(\d+)(?:;\S+)?', value)
    if not match:
        batch.fail(f'Slurm returned unexpected output {value!r} after successful submission. '
                   'Check squeue/sacct before trying again.')
    job_id = match.group(1)
    batch.atomic_json(recorded, {
        'slurm_job_id': job_id,
        'array': array,
        'workspace': str(workspace),
        'manifest': str(manifest_path),
        'field_count': n,
        'max_parallel': throttle,
        'submitted_unix_time': int(time.time()),
        'command': cmd,
    })
    print(f'\nSUBMITTED Slurm array {job_id}: {n} separate field jobs, at most {throttle} running concurrently.')
    print(f'Check queue: squeue -j {job_id}')
    print(f'Logs: {workspace}/logs/ims_{job_id}_<array_index>.out')
    print(f'Inspect fields: python3 {python_runner} inspect --manifest {manifest_path}')
    return job_id


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Prepare AND submit one full IMS-to-results Slurm job per Field_N, in one command.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument('--input', help='Session, group, or individual Field_N directory; discovers fields recursively')
    sources.add_argument('--manifest', help='Already prepared manifest.json; submit without preparing again')
    parser.add_argument('--spec', help='Generic or reporter profile JSON, required with --input')
    parser.add_argument('--workspace', help='NEW workspace path, required with --input')
    parser.add_argument('--repository', default=str(REPOSITORY_ROOT),
                        help='Repository root containing the existing controllers')
    parser.add_argument('--group', help='Optional group glob, e.g. Exp17_157-2_Set2*')
    parser.add_argument('--session', action='append', help='Direct child session name to include; repeatable')
    parser.add_argument('--field', type=int, action='append', help='Select only Field_N (repeatable)')
    parser.add_argument('--max-depth', type=int, default=5, help='Max number of directory levels searched')
    parser.add_argument('--grid', type=int, nargs=2, default=[2,2], metavar=('X','Y'),
                        help='Stitching grid, must match analysis profile')
    incomplete_mode = parser.add_mutually_exclusive_group()
    incomplete_mode.add_argument('--strict-incomplete', action='store_true',
                        help='Abort if any field has fewer than four raw IMS tiles (default: automatically skip and record)')
    incomplete_mode.add_argument('--skip-incomplete', action='store_true',
                        help='Legacy alias; skipping incomplete fields is already the default')
    parser.add_argument('--require-contiguous', action='store_true',
                        help='Require four consecutive original tile numbers; not necessarily F00-F03')
    parser.add_argument('--max-parallel', type=int, default=2,
                        help='Maximum concurrently running field jobs (Slurm %% throttle)')
    parser.add_argument('--run-plan', help='JSON run plan used to select sessions and validate expected counts')
    parser.add_argument('--allow-dirty-repository', action='store_true',
                        help='Allow run-plan preparation from a dirty checkout (not recommended)')
    parser.add_argument('--print-only', action='store_true',
                        help='Prepare the workspace but print submission command without submitting')
    parser.add_argument('--allow-duplicate', action='store_true',
                        help='With --manifest only: permit resubmission after checking prior Slurm job ended')
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        if args.input:
            if not args.spec or not args.workspace:
                parser.error('--input requires both --spec and --workspace')
            if args.allow_duplicate:
                parser.error('--allow-duplicate is only supported with --manifest')
            # Use precisely the existing discovery and preparation implementation.
            prepare_args = argparse.Namespace(
                input=args.input, session=args.session, group=args.group, field=args.field,
                max_depth=args.max_depth, grid=args.grid,
                require_contiguous=args.require_contiguous,
                strict_incomplete=args.strict_incomplete, skip_incomplete=args.skip_incomplete, spec=args.spec,
                workspace=args.workspace, repository=args.repository,
                run_plan=args.run_plan,
                allow_dirty_repository=args.allow_dirty_repository,
            )
            batch.do_prepare(prepare_args)
            manifest_path = batch.abs_path(args.workspace) / 'manifest.json'
        else:
            if (args.spec or args.workspace or args.session or args.group or args.field
                    or args.skip_incomplete or args.strict_incomplete or args.run_plan
                    or args.allow_dirty_repository):
                parser.error('--spec/--workspace/--session/--group/--field/--run-plan only apply to --input')
            manifest_path = batch.abs_path(args.manifest)
        submit(manifest_path, args.max_parallel,
               print_only=args.print_only, allow_duplicate=args.allow_duplicate)
        return 0
    except (ValueError, FileNotFoundError, KeyError, subprocess.CalledProcessError) as e:
        parser.exit(2, f'ERROR: {e}\n')
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
