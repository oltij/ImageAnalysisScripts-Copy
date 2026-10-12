#!/usr/bin/env python3
"""Validate the expected statistical tables, plots, and editable SVG text."""
import argparse
from pathlib import Path
from xml.etree import ElementTree as ET
import pandas as pd

STEMS = ['all_metrics_distributions','sample_sd_by_cell_line','sample_sd_by_experiment',
         'additive_model_residual_qq','residuals_vs_fitted',
         'sample_variance_by_cell_line','sample_variance_by_experiment',
         'sample_variance_cell_line_by_experiment_heatmaps',
         'absolute_deviations_by_cell_line',
         'A_PV_area_variance_detail','A_supp_PV_density_variance_detail',
         'B_BiVe3_specificity_variance_detail','B_BiVe3_coverage_variance_detail',
         'C_BiVe4_specificity_variance_detail','C_BiVe4_coverage_variance_detail']
TABLES = ['descriptive_variances.csv','variance_homogeneity_tests.csv','shapiro_normality.csv',
          'additive_model_residual_diagnostics.csv','per_acquisition_residuals.csv',
          'validated_sample_assignments.csv','variance_matrix_cell_line_experiment.csv']
EXPECTED = {'A_PV_area','A_supp_PV_density','B_BiVe3_specificity','B_BiVe3_coverage',
            'C_BiVe4_specificity','C_BiVe4_coverage'}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results',required=True,type=Path)
    a = p.parse_args()
    d = a.results.expanduser()
    for t in TABLES:
        f=d/t
        if not f.is_file() or f.stat().st_size==0:
            raise SystemExit(f'Missing or empty: {f}')
    desc=pd.read_csv(d/TABLES[0])
    if not EXPECTED.issubset(set(desc.metric)) or not {'sample_variance_ddof1','sample_sd','n_valid'}.issubset(desc.columns):
        raise SystemExit('Descriptive CSV missing metrics/variance columns')
    validated=pd.read_csv(d/'validated_sample_assignments.csv')
    if len(validated)!=35 or validated.acquisition.nunique()!=35:
        raise SystemExit('Expected the exact 35 validated acquisitions')
    matrix = pd.read_csv(d/'variance_matrix_cell_line_experiment.csv')
    if len(matrix) != 54 or set(matrix.metric) != EXPECTED:
        raise SystemExit('Variance matrix should have 6 metrics × 3 cell lines × 3 batches (54 rows)')
    # No variance is permitted for n<2.
    if matrix.loc[matrix.n_valid.lt(2), 'sample_variance_ddof1'].notna().any():
        raise SystemExit('Undefined sample variance (n<2) must be NaN, not 0')
    for stem in STEMS:
        for ext in ('svg','pdf','png'):
            f=d/f'{stem}.{ext}'
            if not f.is_file() or f.stat().st_size<=0:
                raise SystemExit(f'Missing/empty: {f}')
        text_nodes=ET.parse(d/f'{stem}.svg').getroot().findall('.//{http://www.w3.org/2000/svg}text')
        if len(text_nodes)<3:
            raise SystemExit(f'Expected editable text in SVG: {stem}')
        print(f'PASS {stem}: SVG/PDF/PNG, {len(text_nodes)} editable text nodes')
    print(f'PASS all {len(TABLES)} CSV outputs, 6 metrics, 35 acquisitions, {len(STEMS)*3} figure files')

if __name__=='__main__': main()
