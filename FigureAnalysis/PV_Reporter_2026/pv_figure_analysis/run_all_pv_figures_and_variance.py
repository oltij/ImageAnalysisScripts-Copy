#!/usr/bin/env python3
"""One-command, read-only PV publication-figure + variance/assumption audit.

Runs BOTH publication V1 and V2 (A, A-supp, B, C, ABC, BC-only), then the
full A/B/C descriptive/variance figures. Does not re-segment or re-colocalize.
This controller and underlying scripts do not modify the input metrics CSV.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def command(*args):
    line = [sys.executable, '-u', *(str(arg) for arg in args)]
    print('\n>>> ' + ' '.join(line), flush=True)
    subprocess.run(line, check=True, cwd=ROOT)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-dir', type=Path,
                        default=Path.home()/'neuron_ims_batch50_maskfix',
                        help='Folder with prepared.json; resolves to scratch/_PV_selected_figures')
    parser.add_argument('--include-inference', action='store_true',
                        help='Add exploratory Panel A ANOVA, Welch, Kruskal-Wallis and HC3 results')
    parser.add_argument('--stats-only', action='store_true',
                        help='Generate only variance/statistical figures; leave publication plots untouched')
    parser.add_argument('--publication-only', action='store_true',
                        help='Generate only the publication figures, without variance reports')
    parser.add_argument('--skip-verification', action='store_true')
    args = parser.parse_args(argv)
    if args.stats_only and args.publication_only:
        parser.error('Choose either --stats-only or --publication-only, not both.')
    batch = args.batch_dir.expanduser().resolve()
    prep_path = batch/'prepared.json'
    if not prep_path.is_file():
        parser.error(f'Cannot find {prep_path}. Did you install the maskfix batch?')
    cfg = json.loads(prep_path.read_text(encoding='utf-8'))
    stage_root = cfg.get('stage_root')
    if not stage_root:
        parser.error(f'{prep_path} lacks stage_root.')
    figdir = (Path(stage_root)/'_PV_selected_figures').resolve()
    csv = figdir/'per_acquisition_metrics.csv'
    if not csv.is_file():
        parser.error(f'Cannot find {csv}. Generate the 35-organoid measurements first.')
    # Both underlying scripts revalidate all 35 acquisition-to-block mappings.
    if not args.stats_only:
        command(ROOT/'pv_figures.py', '--batch-dir', batch,
                '--strict', '--redraw-only', '--no-normality')
        if not args.skip_verification:
            command(ROOT/'verify_figure_exports.py', '--figures', figdir)
    if not args.publication_only:
        argv = [ROOT/'pv_variance_diagnostics.py', '--batch-dir', batch, '--strict']
        if args.include_inference:
            argv.append('--include-inference')
        command(*argv)
        if not args.skip_verification:
            command(ROOT/'verify_variance_outputs.py', '--results', figdir/'variance_diagnostics')
    print('\nFINISHED: all outputs derived from existing measurements (no image reprocessing).')
    print('Publication figures:', figdir)
    print('Variance / statistics figures:', figdir/'variance_diagnostics')
    print('CAUTION: check that each acquisition is an independent biological organoid before interpreting p-values.')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f'ERROR: Step failed with exit code {exc.returncode}', file=sys.stderr)
        raise SystemExit(exc.returncode)
