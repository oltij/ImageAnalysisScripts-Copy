#!/usr/bin/env python3
"""Read-only descriptive variance and statistical-assumption audit for PV figure A/B/C.

Uses the same 35-acquisition per_acquisition_metrics.csv as pv_figures.py.
No image reprocessing, no changes to the original measurements or figures.

Dependencies: Python 3.10+, numpy, pandas, scipy, matplotlib.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SELECTION = {
    '167-4': (1, 3, 5), '177-4': (14, 15), '178-4': (17, 21),
    '169-4': (22, 23), '170-4': (24, 25, 26), '171-4': (27, 28, 29),
    '173-4': (30, 31, 32), '174-4': (33, 34, 35),
    '165-4': (36, 37), '166-4': (38, 39, 40),
    '175-4': (41, 42), '176-4': (44, 45),
    '168-4': (46, 47), '172-4': (48, 49, 50),
}
BLOCKS = {
    165: ('24', 'H1 CL23', 'BiVe3'), 166: ('24', 'H1 CL23', 'BiVe4'),
    167: ('24', 'H1 CL32', 'BiVe3'), 168: ('24', 'H1 CL32', 'BiVe4'),
    169: ('RA02', 'H1 CL23', 'BiVe3'), 170: ('RA02', 'H1 CL23', 'BiVe4'),
    171: ('RA02', 'H1 CL32', 'BiVe3'), 172: ('RA02', 'H1 CL32', 'BiVe4'),
    173: ('RA02', 'H1', 'BiVe3'), 174: ('RA02', 'H1', 'BiVe4'),
    175: ('17', 'H1 CL32', 'BiVe3'), 176: ('17', 'H1 CL32', 'BiVe4'),
    177: ('17', 'H1', 'BiVe3'), 178: ('17', 'H1', 'BiVe4'),
}
GROUPS = ['H1', 'H1 CL32', 'H1 CL23']
EXPERIMENTS = ['17', '24', 'RA02']
ENHANCERS = ['BiVe3', 'BiVe4']
COLORS = {'17': '#F7A1A3', '24': '#DF4C50', 'RA02': '#801C31'}
SHAPES = {'H1': 'o', 'H1 CL32': 's', 'H1 CL23': '^'}
# stem, display name, original CSV column, enhancer filter, units
METRICS = [
    ('A_PV_area', 'A: PV-positive area', 'pv_positive_area_percent', None, '%'),
    ('A_supp_PV_density', 'A supplement: PV+ cell density', 'pv_nuclear_positive_density_mm2', None, 'cells/mm²'),
    ('B_BiVe3_specificity', 'B: BiVe3 specificity', 'specificity_percent', 'BiVe3', '%'),
    ('B_BiVe3_coverage', 'B: BiVe3 coverage', 'coverage_percent', 'BiVe3', '%'),
    ('C_BiVe4_specificity', 'C: BiVe4 specificity', 'specificity_percent', 'BiVe4', '%'),
    ('C_BiVe4_coverage', 'C: BiVe4 coverage', 'coverage_percent', 'BiVe4', '%'),
]
REQUIRED = {'experiment', 'acquisition', 'cell_line', 'enhancer',
            'pv_positive_area_percent', 'pv_nuclear_positive_density_mm2',
            'specificity_percent', 'coverage_percent'}


def resolve_csv(args):
    if args.csv:
        return args.csv.expanduser().resolve()
    batch = args.batch_dir.expanduser().resolve()
    prepared = batch / 'prepared.json'
    if not prepared.is_file():
        raise FileNotFoundError(f'Missing {prepared}. Supply --csv /path/to/per_acquisition_metrics.csv instead.')
    cfg = json.loads(prepared.read_text(encoding='utf-8'))
    stage_root = cfg.get('stage_root')
    records = cfg.get('records', [])
    if not stage_root or len(records) != 50:
        raise ValueError('prepared.json must contain stage_root and exactly 50 records.')
    by_index = {int(r['index']): r for r in records}
    if len(by_index) != 50 or set(by_index) != set(range(1, 51)):
        raise ValueError('prepared.json indices must be unique and cover 1–50.')
    for sample, ids in SELECTION.items():
        for idx in ids:
            if by_index[idx]['sample'] != sample:
                raise ValueError(f'Acquisition #{idx} is not {sample} according to prepared.json.')
    return (Path(stage_root) / '_PV_selected_figures' / 'per_acquisition_metrics.csv').resolve()


def read_metrics(path, strict):
    if not path.is_file():
        raise FileNotFoundError(f'Cannot find measurements: {path}')
    df = pd.read_csv(path, dtype={'experiment': str, 'cell_line': str, 'enhancer': str})
    missing = REQUIRED.difference(df.columns)
    if missing:
        raise ValueError(f'Input CSV missing required columns: {sorted(missing)}')
    if df.empty:
        raise ValueError('Input CSV contains no rows.')
    df['acquisition'] = pd.to_numeric(df['acquisition'], errors='raise').astype(int)
    if df.acquisition.duplicated().any():
        raise ValueError('Duplicate acquisition IDs; must not count repeated rows as independent organoids.')
    expected = {n for nums in SELECTION.values() for n in nums}
    present = set(df.acquisition)
    if not present.issubset(expected):
        raise ValueError(f'Unexpected acquisition IDs: {sorted(present - expected)}')
    if strict and present != expected:
        raise ValueError(f'Expected all 35 selected acquisitions: missing {sorted(expected - present)}')
    labels = {n: s for s, nums in SELECTION.items() for n in nums}
    derived_exp = []
    for row in df.itertuples():
        label = str(row.experiment).strip()
        if labels[row.acquisition] != label:
            raise ValueError(f'Acquisition {row.acquisition} should be block {labels[row.acquisition]}, not {label}')
        m = re.fullmatch(r'(\d+)-4', label)
        if not m or int(m.group(1)) not in BLOCKS:
            raise ValueError(f'Invalid biological block {label!r}')
        exp, line, enhancer = BLOCKS[int(m.group(1))]
        if ' '.join(str(row.cell_line).split()) != line or str(row.enhancer).strip() != enhancer:
            raise ValueError(f'Acquisition {row.acquisition}: line/enhancer disagrees with biological block {label}')
        derived_exp.append(exp)
    df['experiment_id'] = derived_exp
    for _, _, col, enhancer, unit in METRICS:
        values = pd.to_numeric(df.loc[df.enhancer.eq(enhancer) if enhancer else slice(None), col], errors='coerce')
        if np.isinf(values).any():
            raise ValueError(f'{col} contains infinite values.')
        if (values < 0).any() or (unit == '%' and (values > 100).any()):
            raise ValueError(f'{col} contains invalid measurements for unit {unit}.')
    return df


def metric_rows(df, metric):
    _, _, col, enhancer, _ = metric
    sub = df.loc[df.enhancer.eq(enhancer) if enhancer else slice(None)].copy()
    sub['measurement'] = pd.to_numeric(sub[col], errors='coerce')
    return sub


def finite(sub):
    v = sub['measurement'].to_numpy(dtype=float)
    return v[np.isfinite(v)]


def descriptives(v, group_count):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    n = len(v)
    mean = float(np.mean(v)) if n else np.nan
    sd = float(np.std(v, ddof=1)) if n >= 2 else np.nan
    return {
        'n_total': group_count, 'n_valid': n, 'n_missing': group_count - n,
        'mean': mean,
        'median': float(np.median(v)) if n else np.nan,
        'sample_variance_ddof1': float(np.var(v, ddof=1)) if n >= 2 else np.nan,
        'sample_sd': sd,
        'sem': sd / np.sqrt(n) if n >= 2 else np.nan,
        'cv_percent': 100 * sd / abs(mean) if n >= 2 and mean != 0 else np.nan,
        'q1': float(np.quantile(v, .25)) if n else np.nan,
        'q3': float(np.quantile(v, .75)) if n else np.nan,
        'iqr': float(np.quantile(v, .75) - np.quantile(v, .25)) if n else np.nan,
        'min': float(np.min(v)) if n else np.nan,
        'max': float(np.max(v)) if n else np.nan,
        'n_zero': int(np.count_nonzero(v == 0)),
        'n_100': int(np.count_nonzero(v == 100)),
    }


def descriptive_rows(df):
    out = []
    for metric in METRICS:
        stem, title, _, enhancer, unit = metric
        m = metric_rows(df, metric)
        dimensions = [(), ('cell_line',), ('experiment_id',),
                      ('cell_line', 'experiment_id')]
        if enhancer is None:
            dimensions += [('enhancer',), ('cell_line', 'enhancer')]
        for dimension in dimensions:
            if dimension:
                iterator = m.groupby(list(dimension), dropna=False, sort=False)
            else:
                iterator = [('all', m)]
            for key, grp in iterator:
                if not isinstance(key, tuple):
                    key = (key,)
                label = ', '.join(f'{field}={val}' for field, val in zip(dimension, key)) if dimension else 'all'
                out.append(dict(metric=stem, metric_label=title, units=unit,
                                grouping=' x '.join(dimension) if dimension else 'overall',
                                group=label, cell_line=grp['cell_line'].iloc[0] if 'cell_line' in dimension else '',
                                experiment_id=grp['experiment_id'].iloc[0] if 'experiment_id' in dimension else '',
                                enhancer=grp['enhancer'].iloc[0] if 'enhancer' in dimension else (enhancer or 'both'),
                                **descriptives(finite(grp), len(grp))))
    return pd.DataFrame(out)


def shapiro_details(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    n = len(v)
    if n < 3:
        return dict(W=np.nan, p=np.nan, note='not testable (n<3)')
    if n > 5000:
        return dict(W=np.nan, p=np.nan, note='n>5000; Q-Q plot preferred')
    if np.ptp(v) == 0:
        return dict(W=np.nan, p=np.nan, note='not testable (constant values)')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        w, p = stats.shapiro(v)
    return dict(W=float(w), p=float(p),
                note='very low power (n<8)' if n < 8 else 'visual checks also needed')


def normality_rows(df):
    rows = []
    for metric in METRICS:
        m = metric_rows(df, metric)
        for axis in ['overall', 'cell_line', 'experiment_id']:
            iterator = [('all', m)] if axis == 'overall' else m.groupby(axis)
            for grp_name, subset in iterator:
                vals = finite(subset)
                d = shapiro_details(vals)
                rows.append(dict(metric=metric[0], grouping=axis, group=grp_name,
                                 n=len(vals), shapiro_W=d['W'], shapiro_p=d['p'], note=d['note']))
    return pd.DataFrame(rows)


def variance_homogeneity(df):
    rows = []
    for metric in METRICS:
        m = metric_rows(df, metric)
        for axis, order in [('cell_line', GROUPS), ('experiment_id', EXPERIMENTS)]:
            grps = [(name, finite(m[m[axis] == name])) for name in order]
            grps = [(name, v) for name, v in grps if len(v) > 0]
            raw = dict(metric=metric[0], grouping=axis,
                       group_names=' | '.join(name for name, _ in grps),
                       group_ns=' | '.join(str(len(v)) for _, v in grps))
            tested = len(grps) >= 2 and all(len(v) >= 2 for _, v in grps) and not all(np.ptp(v) == 0 for _, v in grps)
            if not tested:
                rows.append(dict(**raw, levene_brown_forsythe_F=np.nan,
                                 levene_brown_forsythe_p=np.nan, fligner_chi2=np.nan,
                                 fligner_p=np.nan, max_sample_var_div_min=np.nan,
                                 note='Insufficient n or all groups constant.'))
                continue
            arrays = [v for _, v in grps]
            sample_vars = [np.var(v, ddof=1) for v in arrays]
            minvar = min(sample_vars)
            ratio = max(sample_vars) / minvar if minvar > 0 else np.inf if max(sample_vars) else np.nan
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    lv = stats.levene(*arrays, center='median')
                    fl = stats.fligner(*arrays)
                rows.append(dict(**raw, levene_brown_forsythe_F=float(lv.statistic),
                                 levene_brown_forsythe_p=float(lv.pvalue),
                                 fligner_chi2=float(fl.statistic), fligner_p=float(fl.pvalue),
                                 max_sample_var_div_min=float(ratio),
                                 note='Exploratory; small n limits power; neither test proves equality.'))
            except ValueError as err:
                rows.append(dict(**raw, levene_brown_forsythe_F=np.nan,
                                 levene_brown_forsythe_p=np.nan, fligner_chi2=np.nan,
                                 fligner_p=np.nan, max_sample_var_div_min=float(ratio),
                                 note=f'Cannot evaluate: {err}'))
    return pd.DataFrame(rows)


def design_matrix(m, include_enhancer=False):
    cols = [np.ones(len(m), dtype=float)]
    labels = ['Intercept']
    for group in GROUPS[1:]:
        cols.append((m.cell_line == group).to_numpy(dtype=float))
        labels.append(f'cell_line[{group}]')
    for exp in EXPERIMENTS[1:]:
        cols.append((m.experiment_id == exp).to_numpy(dtype=float))
        labels.append(f'experiment[{exp}]')
    if include_enhancer:
        cols.append((m.enhancer == 'BiVe4').to_numpy(dtype=float))
        labels.append('enhancer[BiVe4]')
    return np.column_stack(cols), labels


def fit_additive_model(df, metric, include_enhancer=False):
    x = metric_rows(df, metric)
    x = x[np.isfinite(x.measurement)].copy()
    if len(x) < 8:
        return None
    X, labels = design_matrix(x, include_enhancer=include_enhancer)
    n, p = X.shape
    if n - p < 3 or np.linalg.matrix_rank(X) != p:
        return None
    y = x.measurement.to_numpy(dtype=float)
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ y
    fitted = X @ beta
    resid = y - fitted
    rss = float(resid @ resid)
    leverages = np.einsum('ij,jk,ik->i', X, inv, X)
    hc3_meat = X.T @ ((resid / (1 - leverages))[:, None] ** 2 * X)
    hc3 = inv @ hc3_meat @ inv
    usual = inv * rss / (n-p)
    x = x.copy()
    x['fitted'] = fitted
    x['residual'] = resid
    x['standardized_residual'] = resid / math.sqrt(rss / (n-p)) if rss else np.nan
    x['leverage'] = leverages
    return dict(data=x, X=X, coef=beta, labels=labels,
                cov_hc3=hc3, cov_ols=usual, dof=n-p, rss=rss, n=n, p=p)


def residual_rows(df):
    diag, detail, models = [], [], {}
    for metric in METRICS:
        mod = fit_additive_model(df, metric)
        if mod is None:
            diag.append(dict(metric=metric[0], n=0, note='model not estimable'))
            continue
        models[metric[0]] = mod
        x = mod['data']
        s = shapiro_details(x.residual.to_numpy())
        diag.append(dict(metric=metric[0], n=mod['n'], model='cell_line + experiment',
                         n_parameters=mod['p'], residual_df=mod['dof'],
                         residual_variance_RSS_div_df=mod['rss'] / mod['dof'],
                         residual_SD=math.sqrt(mod['rss'] / mod['dof']),
                         residual_shapiro_W=s['W'], residual_shapiro_p=s['p'],
                         max_leverage=float(x.leverage.max()), note=s['note']))
        for row in x.itertuples():
            detail.append(dict(metric=metric[0], acquisition=int(row.acquisition),
                               cell_line=row.cell_line, experiment_id=row.experiment_id,
                               enhancer=row.enhancer, observed=float(row.measurement),
                               fitted=float(row.fitted), residual=float(row.residual),
                               standardized_residual=float(row.standardized_residual),
                               leverage=float(row.leverage)))
    return pd.DataFrame(diag), pd.DataFrame(detail), models


def grouped_metric_plot(df, out):
    fig, axes = plt.subplots(2, 3, figsize=(12.1, 6.9), constrained_layout=True)
    for ax, metric in zip(axes.flat, METRICS):
        m = metric_rows(df, metric)
        values = [finite(m[m.cell_line == g]) for g in GROUPS]
        positions = np.arange(3)
        if any(len(v) for v in values):
            plotted = [v if len(v) else np.asarray([np.nan]) for v in values]
            ax.boxplot(plotted, positions=positions, widths=.43, whis=1.5,
                       showfliers=False, patch_artist=True,
                       medianprops={'color':'black'},
                       boxprops={'facecolor':'white','edgecolor':'#444444'})
        for j, line in enumerate(GROUPS):
            sub = m[(m.cell_line == line) & np.isfinite(m.measurement)]
            for row in sub.itertuples():
                # Stable, deterministic horizontal placement based on global image number.
                jitter = (((int(row.acquisition) * 4177) % 97) / 96 - .5) * .24
                ax.scatter(j+jitter, row.measurement, marker=SHAPES[line],
                           s=30, c=COLORS[row.experiment_id], linewidths=0,
                           edgecolors='none', alpha=.96, zorder=3)
        ax.set_xticks(positions, ['H1','CL32','CL23'])
        ax.set_title(metric[1], fontsize=10, loc='left', weight='bold')
        ax.set_ylabel(metric[-1])
        ax.spines[['right', 'top']].set_visible(False)
    handles = [plt.Line2D([],[], marker='o', ms=6, linestyle='None',
                          markerfacecolor=COLORS[k], markeredgecolor='none', label=k)
               for k in EXPERIMENTS]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, -.025), ncol=3,
               title='Experiment ID', frameon=False)
    export(fig, out, 'all_metrics_distributions')


def sd_plot(desc, out, grouping):
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 6.7), constrained_layout=True)
    for ax, metric in zip(axes.flat, METRICS):
        tab = desc[(desc.metric == metric[0]) & (desc.grouping == grouping)]
        labels = GROUPS if grouping == 'cell_line' else EXPERIMENTS
        for j, label in enumerate(labels):
            row = tab[tab.cell_line.eq(label)] if grouping == 'cell_line' else tab[tab.experiment_id.eq(label)]
            if row.empty:
                continue
            d = row.iloc[0]
            if np.isfinite(d.sample_sd):
                ax.bar(j, d.sample_sd, width=.62, color='#999999' if grouping == 'cell_line' else COLORS[label])
                ax.text(j, d.sample_sd, f'n={int(d.n_valid)}\nvar={d.sample_variance_ddof1:.3g}',
                        ha='center', va='bottom', fontsize=7)
        ax.set_title(metric[1], fontsize=10, loc='left', weight='bold')
        ax.set_xticks(range(3), ['H1','CL32','CL23'] if grouping == 'cell_line' else EXPERIMENTS)
        ax.set_ylabel(f'Sample SD ({metric[-1]})')
        ax.spines[['top','right']].set_visible(False)
        ax.margins(y=.22)
    export(fig, out, 'sample_sd_by_cell_line' if grouping == 'cell_line' else 'sample_sd_by_experiment')


def qq_plot(models, out):
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.3), constrained_layout=True)
    for ax, metric in zip(axes.flat, METRICS):
        mod = models.get(metric[0])
        if mod is not None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                stats.probplot(mod['data'].residual.to_numpy(), dist='norm', plot=ax)
            ax.set_title(metric[1] + ' (adjusted residuals)', loc='left', fontsize=10, weight='bold')
        else:
            ax.text(.5,.5,'Model not estimable',ha='center',va='center',transform=ax.transAxes)
    export(fig, out, 'additive_model_residual_qq')


def residual_scatter(models, out):
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.2), constrained_layout=True)
    for ax, metric in zip(axes.flat, METRICS):
        mod = models.get(metric[0])
        if mod is None:
            ax.text(.5,.5,'Model not estimable',ha='center',va='center',transform=ax.transAxes)
            continue
        sub = mod['data']
        for exp in EXPERIMENTS:
            x = sub[sub.experiment_id == exp]
            for row in x.itertuples():
                ax.scatter(row.fitted, row.residual, marker=SHAPES[row.cell_line],
                           c=COLORS[exp], s=30, edgecolors='none', linewidths=0)
        ax.axhline(0, color='black', linewidth=.9)
        ax.set_xlabel('Fitted value')
        ax.set_ylabel('Residual')
        ax.set_title(metric[1], fontsize=10, loc='left', weight='bold')
        ax.spines[['top','right']].set_visible(False)
    export(fig, out, 'residuals_vs_fitted')


def export(fig, out, name):
    for ext, dpi in [('svg', 300), ('pdf', 300), ('png', 300)]:
        fig.savefig(out / f'{name}.{ext}', dpi=dpi if ext != 'png' else 600,
                    bbox_inches='tight')
    plt.close(fig)


# ============================================================================
# Direct variance figures (sample s², ddof=1). Added without changing any
# published measurement, previous plot or inferential calculation.
# ============================================================================


def variance_axis_units(metric):
    return 'percentage points squared' if metric[-1] == '%' else '(cells/mm²) squared'


def variance_matrix(df):
    """Explicitly include empty and singleton line×batch cells as undefined."""
    records = []
    for metric in METRICS:
        m = metric_rows(df, metric)
        for line in GROUPS:
            for experiment in EXPERIMENTS:
                s = m[(m.cell_line == line) & (m.experiment_id == experiment)]
                d = descriptives(finite(s), len(s))
                records.append(dict(metric=metric[0], metric_label=metric[1],
                                    units=metric[-1], cell_line=line,
                                    experiment_id=experiment, **d))
    return pd.DataFrame(records)


def variance_bar_ax(ax, desc, metric, grouping, *, compact=False):
    """Visualize *variance*, not SD; n<2 remains undefined, never zero."""
    stem, title, _, _, units = metric
    labels = GROUPS if grouping == 'cell_line' else EXPERIMENTS
    short = ['H1', 'CL32', 'CL23'] if grouping == 'cell_line' else labels
    rows = desc[(desc.metric == stem) & (desc.grouping == grouping)]
    variances = []
    counts = []
    for label in labels:
        match = rows[rows[grouping].eq(label)]
        row = match.iloc[0] if len(match) else None
        variances.append(float(row.sample_variance_ddof1) if row is not None else np.nan)
        counts.append(int(row.n_valid) if row is not None else 0)
    heights = [value if np.isfinite(value) else 0.0 for value in variances]
    colors = ['#828282', '#A6A6A6', '#C9C9C9'] if grouping == 'cell_line' else [COLORS[e] for e in labels]
    for i, (value, height, count, color) in enumerate(zip(variances, heights, counts, colors)):
        if np.isfinite(value):
            ax.bar(i, height, width=.66, color=color, edgecolor='none', zorder=2)
        else:
            ax.bar(i, 0, width=.66, facecolor='none', edgecolor='#bbbbbb', hatch='///')
    ymax = max(heights, default=0.0)
    ax.set_ylim(0, ymax * 1.50 if ymax > 0 else 1)
    for i, (value, count) in enumerate(zip(variances, counts)):
        annotation = f'n={count}\ns²={value:.3g}' if np.isfinite(value) else f'n={count}\nN/A'
        ax.text(i, (value if np.isfinite(value) else 0.0) + ax.get_ylim()[1]*.025,
                annotation, ha='center', va='bottom', fontsize=6.7 if compact else 7.1)
    ax.set_xticks(range(3), short)
    ax.set_ylabel('Sample variance\n(' + variance_axis_units(metric) + ')', fontsize=8)
    ax.set_title(title, fontsize=9 if compact else 10, loc='left', weight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', labelsize=8)
    ax.yaxis.grid(True, color='#e1e1e1', linewidth=.5, zorder=0)
    ax.set_axisbelow(True)


def variance_heatmap_ax(ax, matrix, metric, *, compact=False):
    """Variance across experiment × line; undefined sample variance shown grey."""
    stem = metric[0]
    m = matrix[matrix.metric.eq(stem)]
    var = np.full((len(GROUPS), len(EXPERIMENTS)), np.nan, float)
    ns = np.zeros(var.shape, int)
    for i, line in enumerate(GROUPS):
        for j, exp in enumerate(EXPERIMENTS):
            row = m[(m.cell_line == line) & (m.experiment_id == exp)]
            if len(row):
                var[i, j] = float(row.iloc[0].sample_variance_ddof1)
                ns[i, j] = int(row.iloc[0].n_valid)
    maxval = np.nanmax(var) if np.any(np.isfinite(var)) else 0.0
    vmin = 0.0
    vmax = maxval if maxval > 0 else 1.0
    cmap = plt.get_cmap('Greys').copy()
    cmap.set_bad('#eeeeee')
    ax.imshow(np.ma.masked_invalid(var), cmap=cmap, vmin=vmin, vmax=vmax,
              aspect='equal', interpolation='nearest')
    ax.set_xticks(range(3), EXPERIMENTS, fontsize=8)
    ax.set_yticks(range(3), ['H1', 'CL32', 'CL23'], fontsize=8)
    ax.set_xlabel('Experiment', fontsize=8)
    ax.set_ylabel('Cell line', fontsize=8)
    ax.set_title(metric[1], fontsize=9 if compact else 10, loc='left', weight='bold')
    ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 3, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.8)
    ax.tick_params(which='minor', length=0)
    for i in range(3):
        for j in range(3):
            value = var[i,j]
            if np.isfinite(value):
                content = f'n={ns[i,j]}\ns²={value:.3g}'
                color = 'white' if (maxval > 0 and value / maxval >= .57) else 'black'
            else:
                content = f'n={ns[i,j]}\nN/A'
                color = '#666666'
            ax.text(j, i, content, ha='center', va='center', fontsize=7.0 if compact else 7.4, color=color)


def variance_overview_plots(desc, matrix, out):
    """Three explicit sample-variance overviews, covering all A/B/C outcomes."""
    for grouping, name in [('cell_line', 'sample_variance_by_cell_line'),
                           ('experiment_id', 'sample_variance_by_experiment')]:
        fig, axes = plt.subplots(2, 3, figsize=(13.7, 7.7), constrained_layout=True)
        for ax, metric in zip(axes.flat, METRICS):
            variance_bar_ax(ax, desc, metric, grouping, compact=False)
        fig.suptitle('Sample variance (s²; n≥2 required) — ' +
                     ('by cell line' if grouping == 'cell_line' else 'by experiment'),
                     fontsize=13, fontweight='bold')
        export(fig, out, name)
    fig, axes = plt.subplots(2, 3, figsize=(13.3, 8.4), constrained_layout=True)
    for ax, metric in zip(axes.flat, METRICS):
        variance_heatmap_ax(ax, matrix, metric)
    fig.suptitle('Sample variance by cell line × experiment (grey = n<2; within-panel shading)',
                 fontsize=12, fontweight='bold')
    export(fig, out, 'sample_variance_cell_line_by_experiment_heatmaps')


def deviation_from_median_plot(df, out):
    """Individual absolute deviations underlying Brown–Forsythe dispersion."""
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.4), constrained_layout=True)
    for ax, metric in zip(axes.flat, METRICS):
        m = metric_rows(df, metric)
        for i, group in enumerate(GROUPS):
            sub = m[(m.cell_line == group) & np.isfinite(m.measurement)].copy()
            if sub.empty:
                continue
            median = sub.measurement.median()
            diffs = np.abs(sub.measurement.to_numpy(dtype=float) - median)
            for row, delta in zip(sub.itertuples(), diffs):
                jitter = (((int(row.acquisition) * 4177) % 97) / 96 - .5) * .30
                ax.scatter(i + jitter, delta, s=32, color=COLORS[row.experiment_id],
                           marker=SHAPES[group], linewidths=0, edgecolors='none', zorder=3)
            ax.hlines(float(np.mean(diffs)), i-.17, i+.17, linewidth=2.2,
                      color='black', zorder=4)
        ax.set_xticks(range(3), ['H1', 'CL32', 'CL23'])
        ax.set_ylabel('Absolute deviation from group median\n(' + metric[-1] + ')', fontsize=8)
        ax.set_title(metric[1], fontsize=10, loc='left', weight='bold')
        ax.spines[['top', 'right']].set_visible(False)
    legend = [plt.Line2D([], [], marker='o', linestyle='none', ms=6,
                         markerfacecolor=COLORS[e], markeredgecolor='none', label=e)
              for e in EXPERIMENTS]
    fig.legend(handles=legend, loc='lower center', ncol=3,
               bbox_to_anchor=(.5, -.023), title='Experiment', frameon=False)
    fig.suptitle('Within-line absolute deviations (black line = mean absolute deviation)',
                 fontsize=12, fontweight='bold')
    export(fig, out, 'absolute_deviations_by_cell_line')


def variance_per_metric_plots(desc, matrix, out):
    """One detailed variance illustration for each of the six measurements."""
    for metric in METRICS:
        fig = plt.figure(figsize=(12.5, 3.85))
        gs = fig.add_gridspec(1, 3, left=.065, right=.987, top=.78, bottom=.20,
                              width_ratios=[1, 1, 1.06], wspace=.46)
        variance_bar_ax(fig.add_subplot(gs[0,0]), desc, metric, 'cell_line', compact=True)
        variance_bar_ax(fig.add_subplot(gs[0,1]), desc, metric, 'experiment_id', compact=True)
        variance_heatmap_ax(fig.add_subplot(gs[0,2]), matrix, metric, compact=True)
        fig.suptitle(metric[1] + ' — sample variance (s²)', fontsize=13, fontweight='bold', y=.975)
        fig.text(.5, .055,
                 'Variance is in squared measurement units. n<2: undefined, not zero. '
                 'Heatmap shades are relative within this outcome.',
                 ha='center', va='center', fontsize=8)
        export(fig, out, f'{metric[0]}_variance_detail')


def welch_anova(groups):
    # Welch 1951 one-way F and Satterthwaite denominator df; no equal-variance assumption.
    groups = [np.asarray(x, float) for x in groups]
    k = len(groups)
    n = np.array([len(x) for x in groups], dtype=float)
    if k < 2 or np.any(n < 2):
        return np.nan, np.nan, np.nan, 'insufficient group sizes'
    variances = np.array([np.var(x, ddof=1) for x in groups])
    if np.any(variances <= 0):
        return np.nan, np.nan, np.nan, 'zero within-group variance; Welch not defined'
    w = n / variances
    ws = w.sum()
    mean = np.array([x.mean() for x in groups])
    wm = float(np.dot(w, mean) / ws)
    # Welch correction and denominator df.
    base = np.sum((1 - w/ws)**2/(n-1))
    f = (np.sum(w*(mean-wm)**2) / (k-1)) / (1+2*(k-2)/(k*k-1)*base)
    df2 = (k*k-1)/(3*base)
    return float(f), float(stats.f.sf(f, k-1, df2)), float(df2), 'Welch unequal-variance ANOVA'


def holm(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    prev = 0.
    for rank, i in enumerate(order):
        prev = max(prev, min(1., p[i]*(m-rank)))
        adj[i] = prev
    return adj


def wald_test(mod, cov):
    # Joint cell-line effect for the two coefficients relative to H1.
    beta = mod['coef'][1:3]
    C = cov[1:3, 1:3]
    if np.linalg.matrix_rank(C) < 2:
        return np.nan, np.nan
    F = float(beta.T @ np.linalg.inv(C) @ beta / 2)
    return F, float(stats.f.sf(F, 2, mod['dof']))


def a_inference(df):
    metric = METRICS[0]
    m = metric_rows(df, metric)
    subgroups = [finite(m[m.cell_line == g]) for g in GROUPS]
    rows = []
    if all(len(x) >= 2 for x in subgroups):
        f = stats.f_oneway(*subgroups)
        rows.append(dict(test='One-way ordinary ANOVA', F_or_H=float(f.statistic), p=float(f.pvalue),
                         df_num=2, df_den=sum(map(len, subgroups))-3, note='Unadjusted; assumes equal group variances'))
        wf, wp, wd2, msg = welch_anova(subgroups)
        rows.append(dict(test="Welch's one-way ANOVA", F_or_H=wf, p=wp,
                         df_num=2, df_den=wd2, note=msg+'; unadjusted for experiment'))
        kw = stats.kruskal(*subgroups)
        rows.append(dict(test='Kruskal-Wallis H', F_or_H=float(kw.statistic), p=float(kw.pvalue),
                         df_num=2, df_den=np.nan, note='Unadjusted rank-distribution test; not necessarily a median test'))
    contrasts = []
    for with_enhancer in (False, True):
        mod = fit_additive_model(df, metric, include_enhancer=with_enhancer)
        if mod is None:
            rows.append(dict(test='Batch-adjusted linear model', F_or_H=np.nan,
                             p=np.nan, note='model rank-deficient/not estimable'))
            continue
        modname = 'cell_line + experiment' + (' + enhancer' if with_enhancer else '')
        for cov_label, cov in [('ordinary covariance', mod['cov_ols']), ('HC3 robust covariance', mod['cov_hc3'])]:
            fval, pval = wald_test(mod, cov)
            rows.append(dict(test=f'Joint cell-line effect: {modname} ({cov_label})',
                             F_or_H=fval, p=pval, df_num=2, df_den=mod['dof'],
                             note='Approximate F test, conditional on additive model; independent experimental units required'))
        combos = [ ('H1 CL32 - H1', np.array([0, 1, 0])),
                   ('H1 CL23 - H1', np.array([0, 0, 1])),
                   ('H1 CL32 - H1 CL23', np.array([0, 1, -1])) ]
        local = []
        for label, vec in combos:
            c = np.zeros(mod['p'])
            c[:3] = vec
            difference = float(c @ mod['coef'])
            se = math.sqrt(max(0., float(c @ mod['cov_hc3'] @ c)))
            t = difference / se if se > 0 else np.nan
            pv = float(2 * stats.t.sf(abs(t), mod['dof'])) if np.isfinite(t) else np.nan
            ci = stats.t.ppf(.975, mod['dof']) * se
            local.append(dict(model=modname, comparison=label,
                              adjusted_difference_percentage_points=difference,
                              hc3_se=se, t=t, df=mod['dof'], raw_p=pv,
                              ci95_lower=difference-ci, ci95_upper=difference+ci,
                              note='Exploratory additive-model HC3 contrast; Holm across 3 cell-line comparisons'))
        good = [j for j, d in enumerate(local) if np.isfinite(d['raw_p'])]
        if good:
            corrected = holm([local[j]['raw_p'] for j in good])
            for j, ap in zip(good, corrected):
                local[j]['holm_p'] = float(ap)
        contrasts.extend(local)
    return pd.DataFrame(rows), pd.DataFrame(contrasts)


def write_summary(out, src, df, desc, homog, residual, inference_enabled):
    lines = [
        '# PV figure variance and assumption checks', '',
        f'Source CSV: `{src}`', f'Acquisitions present: **{len(df)}** (maximum selected: 35)',
        'Variance = **sample variance** (`ddof=1`, units squared). SD is the square root of the variance.',
        'For A, PV area is a percentage, so SD is in percentage points and sample variance in percentage-points squared.',
        '', '## Variance by cell line', '',
        '| Metric | H1 (n; SD; variance) | CL32 (n; SD; variance) | CL23 (n; SD; variance) |',
        '|---|---|---|---|'
    ]
    for metric in METRICS:
        cells = []
        for name in GROUPS:
            x = desc[(desc.metric==metric[0]) & (desc.grouping=='cell_line') & (desc.cell_line==name)]
            if x.empty:
                cells.append('missing')
            else:
                row = x.iloc[0]
                cells.append(f"{int(row.n_valid)}; {row.sample_sd:.4g}; {row.sample_variance_ddof1:.4g}")
        lines.append(f"| {metric[1]} | {' | '.join(cells)} |")
    lines += ['', '## Brown–Forsythe/Levene checks', '',
              'Levene test with `center=median` (Brown–Forsythe variant). Small groups provide limited power.',
              'A high p-value **does not establish equal variances**. The experiment-adjusted residuals are also relevant.', '',
              '| Metric | Grouping | Brown–Forsythe p | Fligner p | Max/min variance |',
              '|---|---|---:|---:|---:|']
    for row in homog.itertuples():
        lines.append(f'| {row.metric} | {row.grouping} | {row.levene_brown_forsythe_p:.4g} | {row.fligner_p:.4g} | {row.max_sample_var_div_min:.4g} |')
    lines += ['', '## Additive-model residual checks', '',
              'Model for each metric: value ~ cell_line + experiment ID (independent observations assumed).',
              'This is an exploratory diagnostic, not evidence that assumptions hold.', '',
              '| Metric | Residual df | Residual SD | Residual Shapiro p |',
              '|---|---:|---:|---:|']
    for row in residual.itertuples():
        lines.append(f'| {row.metric} | {getattr(row,"residual_df",float("nan"))} | {getattr(row,"residual_SD",float("nan")):.4g} | {getattr(row,"residual_shapiro_p",float("nan")):.4g} |')
    lines += ['', '## Direct variance figures', '',
              '- `sample_variance_by_cell_line`: sample variance for every A, B and C outcome grouped by cell line.',
              '- `sample_variance_by_experiment`: sample variance for every outcome grouped by batch.',
              '- `sample_variance_cell_line_by_experiment_heatmaps`: sample variances for each line × batch, with sample sizes and undefined values identified.',
              '- `absolute_deviations_by_cell_line`: absolute deviations from within-line medians; mean absolute-deviation bars; relates to Brown–Forsythe.',
              '- Six `*_variance_detail` figures give a full-size, outcome-specific view of all three variance breakdowns.',
              '- `variance_matrix_cell_line_experiment.csv` includes missing cells and n=1 groups explicitly. Variance is *undefined* for n<2; it is not zero.',
              '', '## Interpretation notes', '',
              '- The experiment × cell-line cells are sparse or absent; use `descriptive_variances.csv` to inspect them. Do not use missing cells as zero observations.',
              '- For A, percentages are pooled across BiVe3 and BiVe4 as in the figure; an optional enhancer-adjusted analysis is included when `--include-inference` is set.',
              '- For B/C, specificity and coverage are distinct outcomes; each uses the 17 or 18 organoids with the relevant enhancer.',
              '- Variances for metrics in different units cannot be compared numerically.',
              '- Most cell-line subsets are too small for a decisive normality test. Inspect Q–Q plots and residuals, not just Shapiro p-values.',
              '- Very high or low values near the 0%/100% boundaries can create mean–variance dependence.',
              '- **Important:** An acquisition is not automatically an independent biological replicate. If multiple fields/sections originate from one organoid, these tests pseudoreplicate. Confirm independence before publishing inferential p-values.',
              '- Experiment ID may denote a small number of biological batches. Fixed adjustment does not establish generalizability across batches.',
              '- All outputs are descriptive/exploratory. Choose the primary statistical model and multiplicity correction before formal significance reporting.',
              '', f'A univariate exploratory comparison for A was {"enabled" if inference_enabled else "not run (add --include-inference to enable)"}.', '']
    (out/'README_results.md').write_text('\n'.join(lines), encoding='utf-8')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch-dir', type=Path, default=Path.home()/'neuron_ims_batch50_maskfix',
                   help='Folder containing prepared.json from the maskfix batch')
    p.add_argument('--csv', type=Path, help='Override: input per_acquisition_metrics.csv path')
    p.add_argument('--output-dir', type=Path,
                   help='Default: folder named variance_diagnostics next to the input CSV')
    p.add_argument('--strict', action='store_true', help='Require the exact 35 expected acquisitions')
    p.add_argument('--include-inference', action='store_true',
                   help='Additionally compute exploratory Panel A ordinary, Welch, Kruskal and batch-adjusted HC3 tests')
    p.add_argument('--no-plots', action='store_true', help='Generate CSV/report only')
    args = p.parse_args(argv)
    src = resolve_csv(args)
    df = read_metrics(src, args.strict)
    out = args.output_dir.expanduser().resolve() if args.output_dir else src.parent / 'variance_diagnostics'
    out.mkdir(parents=True, exist_ok=True)
    desc = descriptive_rows(df)
    normal = normality_rows(df)
    homo = variance_homogeneity(df)
    diagnostics, detail, models = residual_rows(df)
    matrix = variance_matrix(df)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        desc.to_csv(out/'descriptive_variances.csv', index=False)
        matrix.to_csv(out/'variance_matrix_cell_line_experiment.csv', index=False)
        normal.to_csv(out/'shapiro_normality.csv', index=False)
        homo.to_csv(out/'variance_homogeneity_tests.csv', index=False)
        diagnostics.to_csv(out/'additive_model_residual_diagnostics.csv', index=False)
        detail.to_csv(out/'per_acquisition_residuals.csv', index=False)
        df[['acquisition','experiment','experiment_id','cell_line','enhancer']].to_csv(out/'validated_sample_assignments.csv',index=False)
        if args.include_inference:
            tests, contrasts = a_inference(df)
            tests.to_csv(out/'A_exploratory_omnibus_tests.csv',index=False)
            contrasts.to_csv(out/'A_exploratory_HC3_pairwise.csv',index=False)
    if not args.no_plots:
        grouped_metric_plot(df, out)
        sd_plot(desc, out, 'cell_line')
        sd_plot(desc, out, 'experiment_id')
        qq_plot(models, out)
        residual_scatter(models, out)
        variance_overview_plots(desc, matrix, out)
        deviation_from_median_plot(df, out)
        variance_per_metric_plots(desc, matrix, out)
    write_summary(out, src, df, desc, homo, diagnostics, args.include_inference)
    print('\nPV figure A/B/C variance and assumption analysis completed.')
    print(f'Source (read only): {src}')
    print(f'Output: {out}')
    print(f'Acquisitions validated: {len(df)}; plots: {"none" if args.no_plots else "15 x SVG/PDF/PNG"}')
    print('Summary: README_results.md; sample variance: descriptive_variances.csv')
    print('Direct variance figures: cell-line and experiment overviews, heatmap, 6 individual metric details.')
    print('Brown-Forsythe and Fligner: variance_homogeneity_tests.csv')
    print('Residual checks: additive_model_residual_diagnostics.csv')
    if args.include_inference:
        print('Exploratory A inference: A_exploratory_omnibus_tests.csv and A_exploratory_HC3_pairwise.csv')
    print('NOTE: Confirm independent organoids before interpreting inferential p-values.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, FileNotFoundError, KeyError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(2)
