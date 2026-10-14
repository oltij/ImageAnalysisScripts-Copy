#!/usr/bin/env python3
"""Submit the validated 94-field IGF manifest as two dependent arrays."""
from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import run_folder_batch as batch


FIRST_ARRAY = '0-46%47'
SECOND_ARRAY = '47-93%47'
EXPECTED_FIELDS = 94
RECORD_NAME = 'slurm_two_batch_submission.json'


def parse_job_id(result: subprocess.CompletedProcess[str]) -> str:
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)
    value = result.stdout.strip()
    match = re.fullmatch(r'(\d+)(?:;\S+)?', value)
    if not match:
        batch.fail(
            f'Slurm returned unexpected output {value!r} after successful submission. '
            'Check squeue/sacct before retrying.'
        )
    return match.group(1)


def load_context(manifest_value: str | Path) -> dict[str, object]:
    manifest_path = batch.abs_path(manifest_value)
    manifest = batch.load_manifest(manifest_path)
    summary = batch.validate_manifest_for_submission(manifest, check_sources=True)
    workspace = Path(manifest['workspace']).resolve(strict=False)
    if manifest_path != workspace / 'manifest.json':
        batch.fail('Manifest must be located at <workspace>/manifest.json')
    if len(manifest['jobs']) != EXPECTED_FIELDS:
        batch.fail(
            f'Production two-batch submission requires exactly {EXPECTED_FIELDS} fields; '
            f'found {len(manifest["jobs"])}'
        )
    plan = manifest.get('run_plan')
    if not plan or not plan.get('validation', {}).get('passed'):
        batch.fail('Production two-batch submission requires a validated run-plan manifest')
    runner = Path(manifest['slurm_runner']).resolve(strict=False)
    python_runner = Path(manifest['batch_runner']).resolve(strict=False)
    if not runner.is_file() or not runner.is_absolute():
        batch.fail(f'Invalid absolute Slurm runner path: {runner}')
    if not python_runner.is_file() or not python_runner.is_absolute():
        batch.fail(f'Invalid absolute Python runner path: {python_runner}')
    return {
        'manifest_path': manifest_path,
        'manifest': manifest,
        'workspace': workspace,
        'runner': runner,
        'python_runner': python_runner,
        'summary': summary,
    }


def slurm_command(
    context: dict[str, object], array: str, *, sbatch: str,
    dependency: str | None = None,
) -> list[str]:
    command = [
        sbatch,
        '--parsable',
        f'--array={array}',
        f'--chdir={context["workspace"]}',
    ]
    if dependency is not None:
        command.append(f'--dependency=afterok:{dependency}')
    command.extend([
        str(context['runner']),
        str(context['manifest_path']),
        str(context['python_runner']),
    ])
    return command


def submit_two_batches(
    manifest_path: str | Path, *, print_only: bool = False,
    allow_duplicate: bool = False, resume_second_batch: bool = False,
    sbatch: str = 'sbatch',
) -> tuple[str, str] | None:
    context = load_context(manifest_path)
    workspace = context['workspace']
    record_path = workspace / RECORD_NAME
    previous = json.loads(record_path.read_text()) if record_path.exists() else None

    if resume_second_batch:
        if not previous or not previous.get('batch1_job_id') or previous.get('batch2_job_id'):
            batch.fail('--resume-second-batch requires a partial record with only Batch 1 submitted')
    elif previous and not (print_only or allow_duplicate):
        batch.fail(
            f'Two-batch submission already recorded in {record_path}. '
            'Use --resume-second-batch only for a recorded partial submission, or '
            '--allow-duplicate only after confirming all prior jobs ended.'
        )

    first_preview = slurm_command(context, FIRST_ARRAY, sbatch=sbatch)
    second_preview = slurm_command(
        context, SECOND_ARRAY, sbatch=sbatch, dependency='<BATCH1_JOB_ID>',
    )
    print(f'Validated manifest: {context["summary"]}')
    print('Batch 1 command: ' + shlex.join(first_preview))
    print('Batch 2 command: ' + shlex.join(second_preview))
    print('Batch 2 is held by afterok until every Batch 1 array task succeeds.')
    if print_only:
        print('PRINT ONLY: no Slurm job submitted and no submission record written.')
        return None
    if not shutil.which(sbatch):
        batch.fail(f'Could not find Slurm command {sbatch!r}')

    history = []
    if previous:
        history = list(previous.get('previous_submissions', []))
        if allow_duplicate:
            history.append({key: value for key, value in previous.items()
                            if key != 'previous_submissions'})

    if resume_second_batch:
        batch1_id = str(previous['batch1_job_id'])
        first_command = previous['batch1_command']
        submitted_at = previous.get('submitted_unix_time', int(time.time()))
    else:
        first_command = slurm_command(context, FIRST_ARRAY, sbatch=sbatch)
        first_result = subprocess.run(
            first_command, cwd=workspace, text=True,
            capture_output=True, check=True,
        )
        batch1_id = parse_job_id(first_result)
        submitted_at = int(time.time())
        batch.atomic_json(record_path, {
            'status': 'batch1_submitted',
            'manifest': str(context['manifest_path']),
            'workspace': str(workspace),
            'batch1_job_id': batch1_id,
            'batch2_job_id': None,
            'batch1_array': FIRST_ARRAY,
            'batch2_array': SECOND_ARRAY,
            'batch1_command': first_command,
            'batch2_command': None,
            'submitted_unix_time': submitted_at,
            'previous_submissions': history,
        })

    second_command = slurm_command(
        context, SECOND_ARRAY, sbatch=sbatch, dependency=batch1_id,
    )
    try:
        second_result = subprocess.run(
            second_command, cwd=workspace, text=True,
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError:
        print(
            f'Batch 1 job {batch1_id} was submitted, but Batch 2 submission failed. '
            f'The partial record is {record_path}; repair the cause and use '
            '--resume-second-batch to avoid duplicating Batch 1.',
            file=sys.stderr,
        )
        raise
    batch2_id = parse_job_id(second_result)
    batch.atomic_json(record_path, {
        'status': 'both_submitted',
        'manifest': str(context['manifest_path']),
        'workspace': str(workspace),
        'batch1_job_id': batch1_id,
        'batch2_job_id': batch2_id,
        'batch1_array': FIRST_ARRAY,
        'batch2_array': SECOND_ARRAY,
        'batch2_dependency': f'afterok:{batch1_id}',
        'batch1_command': first_command,
        'batch2_command': second_command,
        'submitted_unix_time': submitted_at,
        'completed_submission_unix_time': int(time.time()),
        'previous_submissions': history,
    })
    print(f'SUBMITTED Batch 1 job {batch1_id} and dependent Batch 2 job {batch2_id}.')
    print(f'Record: {record_path}')
    print(f'Monitor: squeue -j {batch1_id},{batch2_id}')
    return batch1_id, batch2_id


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, help='Validated 94-field manifest.json')
    parser.add_argument('--print-only', action='store_true', help='Print both commands without submitting')
    parser.add_argument('--resume-second-batch', action='store_true', help='Submit only Batch 2 from a recorded partial submission')
    parser.add_argument('--allow-duplicate', action='store_true', help='Submit a new pair only after confirming prior jobs ended')
    parser.add_argument('--sbatch', default='sbatch', help=argparse.SUPPRESS)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    if args.resume_second_batch and args.allow_duplicate:
        parser.error('--resume-second-batch and --allow-duplicate are mutually exclusive')
    try:
        submit_two_batches(
            args.manifest,
            print_only=args.print_only,
            allow_duplicate=args.allow_duplicate,
            resume_second_batch=args.resume_second_batch,
            sbatch=args.sbatch,
        )
    except (ValueError, FileNotFoundError, KeyError, subprocess.CalledProcessError) as exc:
        parser.exit(2, f'ERROR: {exc}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
