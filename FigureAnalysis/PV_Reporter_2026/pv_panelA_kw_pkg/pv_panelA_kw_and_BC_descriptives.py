#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

# ---------- helpers ----------

REDS = ["#f3b0ae", "#d96a6a", "#8f1d1d"]
EXPERIMENT_ORDER = ["17", "24", "RA02"]
CELL_LINE_ORDER = ["H1", "H1 CL32", "H1 CL23"]
SHAPES = {"H1": "o", "H1 CL32": "s", "H1 CL23": "^"}


def load_stage_root(batch_dir: Path) -> Path:
    prepared = batch_dir / "prepared.json"
    if prepared.exists():
        try:
            data = json.loads(prepared.read_text())
            stage_root = Path(data["stage_root"])
            if stage_root.exists():
                return stage_root
        except Exception:
            pass
    if (batch_dir / "_PV_selected_figures").exists():
        return batch_dir
    # maybe batch_dir is actually stage_root sibling parent; use as-is
    return batch_dir


def discover_metrics_csv(stage_root: Path) -> Path:
    candidates = [
        stage_root / "_PV_selected_figures" / "per_acquisition_metrics.csv",
        stage_root / "_PV_selected_figures" / "selected_metrics.csv",
        stage_root / "per_acquisition_metrics.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    figdir = stage_root / "_PV_selected_figures"
    if figdir.exists():
        matches = list(figdir.rglob("*.csv"))
        for m in matches:
            low = m.name.lower()
            if "per_acquisition" in low or "selected" in low or "metric" in low:
                return m
    matches = list(stage_root.rglob("*.csv"))
    for m in matches:
        low = m.name.lower()
        if "per_acquisition" in low or "selected" in low or "metric" in low:
            return m
    raise FileNotFoundError("Could not find per_acquisition_metrics.csv under stage_root/_PV_selected_figures")


def first_present(df: pd.DataFrame, names: Sequence[str]) -> Optional[str]:
    cols = {c.lower(): c for c in df.columns}
    for name in names:
        if name.lower() in cols:
            return cols[name.lower()]
    return None


def require_col(df: pd.DataFrame, names: Sequence[str], label: str) -> str:
    c = first_present(df, names)
    if c is None:
        raise KeyError(f"Could not find column for {label}. Tried: {names}")
    return c


def holm_adjust(pvals: Sequence[float]) -> List[float]:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    adjusted = np.empty(n, dtype=float)
    running = 0.0
    for rank, idx in enumerate(order):
        val = (n - rank) * pvals[idx]
        running = max(running, val)
        adjusted[idx] = min(1.0, running)
    return adjusted.tolist()


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "p = NA"
    if p < 1e-4:
        return "p < 1e-4"
    return f"p = {p:.4f}"


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_all_formats(fig: plt.Figure, out_base: Path, dpi: int = 600) -> None:
    fig.savefig(out_base.with_suffix('.svg'), bbox_inches='tight')
    fig.savefig(out_base.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(out_base.with_suffix('.png'), dpi=dpi, bbox_inches='tight')


def iqr_vals(x: Iterable[float]) -> Tuple[float,float,float]:
    arr = np.asarray(list(x), dtype=float)
    q1 = np.nanpercentile(arr, 25)
    med = np.nanpercentile(arr, 50)
    q3 = np.nanpercentile(arr, 75)
    return q1, med, q3


def desc_stats(series: pd.Series) -> Dict[str, float]:
    x = pd.to_numeric(series, errors='coerce').dropna().astype(float)
    q1, med, q3 = iqr_vals(x)
    return {
        'n': int(x.shape[0]),
        'mean': float(x.mean()),
        'sd': float(x.std(ddof=1)) if len(x) > 1 else float('nan'),
        'variance': float(x.var(ddof=1)) if len(x) > 1 else float('nan'),
        'sem': float(x.std(ddof=1) / math.sqrt(len(x))) if len(x) > 1 else float('nan'),
        'median': float(med),
        'q1': float(q1),
        'q3': float(q3),
        'iqr': float(q3 - q1),
        'min': float(x.min()),
        'max': float(x.max()),
    }


def make_metric_from_counts(df: pd.DataFrame, numerator_col: str, denom_col: str, as_percent: bool = True) -> pd.Series:
    num = pd.to_numeric(df[numerator_col], errors='coerce')
    den = pd.to_numeric(df[denom_col], errors='coerce')
    out = num / den
    out = out.where(den > 0)
    return out * 100.0 if as_percent else out


def maybe_compute_density(df: pd.DataFrame) -> Tuple[pd.Series, str]:
    c = first_present(df, [
        'A_supp_PV_density', 'pv_density', 'pv_positive_density', 'pv_positive_cells_per_mm2',
        'pv_nuclear_density', 'pv_hoechst_colocalized_density', 'pv_positive_nuclear_density_per_mm2'
    ])
    if c:
        return pd.to_numeric(df[c], errors='coerce'), c
    count_col = first_present(df, [
        'pv_nuclear_count', 'pv_positive_nuclear_count', 'pv_hoechst_colocalized_count',
        'pv_positive_cells', 'pv_colocalized_count'
    ])
    area_col = first_present(df, [
        'organoid_area_mm2', 'hoechst_organoid_area_mm2', 'projected_organoid_area_mm2'
    ])
    if count_col and area_col:
        return pd.to_numeric(df[count_col], errors='coerce') / pd.to_numeric(df[area_col], errors='coerce'), f'{count_col}/{area_col}'
    raise KeyError('Could not find or compute density metric for panel A.')


def maybe_compute_nuclear_fraction(df: pd.DataFrame) -> Tuple[pd.Series, str]:
    c = first_present(df, [
        'pv_positive_of_hoechst_percent', 'pv_nuclear_fraction_percent', 'pv_cells_pct_of_hoechst',
        'pv_hoechst_fraction_percent', 'pv_percent_hoechst'
    ])
    if c:
        return pd.to_numeric(df[c], errors='coerce'), c
    num_col = first_present(df, [
        'pv_nuclear_count', 'pv_positive_nuclear_count', 'pv_hoechst_colocalized_count',
        'pv_positive_cells', 'pv_colocalized_count'
    ])
    den_col = first_present(df, [
        'hoechst_count', 'hoechst_nuclear_count', 'nuclei_count', 'hoechst_positive_count'
    ])
    if num_col and den_col:
        return make_metric_from_counts(df, num_col, den_col, as_percent=True), f'100*{num_col}/{den_col}'
    raise KeyError('Could not find or compute nuclear PV fraction metric (PV+/Hoechst+ %).')


def get_BC_metrics(df: pd.DataFrame) -> Dict[str, str]:
    mapping = {}
    for label, names in {
        'B_BiVe3_specificity': ['B_BiVe3_specificity', 'bive3_specificity', 'bive3_specificity_percent'],
        'B_BiVe3_coverage': ['B_BiVe3_coverage', 'bive3_coverage', 'bive3_coverage_percent'],
        'C_BiVe4_specificity': ['C_BiVe4_specificity', 'bive4_specificity', 'bive4_specificity_percent'],
        'C_BiVe4_coverage': ['C_BiVe4_coverage', 'bive4_coverage', 'bive4_coverage_percent'],
    }.items():
        c = first_present(df, names)
        if c:
            mapping[label] = c
    missing = [k for k in ['B_BiVe3_specificity','B_BiVe3_coverage','C_BiVe4_specificity','C_BiVe4_coverage'] if k not in mapping]
    if missing:
        raise KeyError(f'Missing B/C metrics: {missing}. Available columns: {list(df.columns)}')
    return mapping


def pairwise_mwu(df: pd.DataFrame, group_col: str, value_col: str, order: Sequence[str]) -> pd.DataFrame:
    rows = []
    raw = []
    pairs = []
    for i in range(len(order)):
        for j in range(i+1, len(order)):
            g1, g2 = order[i], order[j]
            x = pd.to_numeric(df.loc[df[group_col]==g1, value_col], errors='coerce').dropna().values
            y = pd.to_numeric(df.loc[df[group_col]==g2, value_col], errors='coerce').dropna().values
            if len(x) == 0 or len(y) == 0:
                p = np.nan; u = np.nan
            else:
                u, p = stats.mannwhitneyu(x, y, alternative='two-sided', method='auto')
            pairs.append((g1, g2))
            raw.append(p)
            rows.append({'group1': g1, 'group2': g2, 'n1': len(x), 'n2': len(y), 'U': u, 'p_raw': p})
    adj = holm_adjust([1.0 if pd.isna(p) else p for p in raw])
    for r, p_adj in zip(rows, adj):
        r['p_holm'] = p_adj
    return pd.DataFrame(rows)


def overall_kw(df: pd.DataFrame, group_col: str, value_col: str, order: Sequence[str]) -> Dict[str, float]:
    groups = [pd.to_numeric(df.loc[df[group_col]==g, value_col], errors='coerce').dropna().values for g in order]
    groups = [g for g in groups if len(g) > 0]
    stat, p = stats.kruskal(*groups, nan_policy='omit')
    return {'H': float(stat), 'p': float(p)}


def add_stat_bracket(ax, x1, x2, y, text, h, fontsize=9):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.1, c='black')
    ax.text((x1+x2)/2, y+h, text, ha='center', va='bottom', fontsize=fontsize)


def panel_plot(df: pd.DataFrame, group_col: str, value_col: str, y_label: str, title: str,
               out_base: Path, pairwise: pd.DataFrame, kw: Dict[str,float],
               use_shapes: bool = True) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.1))
    positions = np.arange(1, len(CELL_LINE_ORDER)+1)
    data = [pd.to_numeric(df.loc[df[group_col]==g, value_col], errors='coerce').dropna().values for g in CELL_LINE_ORDER]

    bp = ax.boxplot(data, positions=positions, widths=0.5, patch_artist=True, whis=1.5, showfliers=False)
    for box in bp['boxes']:
        box.set(facecolor='white', edgecolor='black', linewidth=1.2)
    for med in bp['medians']:
        med.set(color='black', linewidth=1.3)
    for w in bp['whiskers'] + bp['caps']:
        w.set(color='black', linewidth=1.0)

    exp_color = {e: REDS[i] for i, e in enumerate(EXPERIMENT_ORDER)}
    rng = np.random.default_rng(0)
    for i, group in enumerate(CELL_LINE_ORDER, start=1):
        sub = df.loc[df[group_col]==group].copy()
        if sub.empty:
            continue
        xs = i + rng.uniform(-0.12, 0.12, size=len(sub))
        for (_, row), x in zip(sub.iterrows(), xs):
            exp = str(row['experiment_id'])
            color = exp_color.get(exp, '#999999')
            marker = SHAPES[group] if use_shapes else 'o'
            y = pd.to_numeric(row[value_col], errors='coerce')
            ax.scatter(x, y, s=42, marker=marker, facecolor=color, edgecolor='none', zorder=3)

        vals = pd.to_numeric(sub[value_col], errors='coerce').dropna().values
        if len(vals):
            mean = float(np.mean(vals))
            sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
            ax.errorbar(i, mean, yerr=sd, fmt='_', markersize=22, color='black', elinewidth=2.0, capsize=4, capthick=1.5, zorder=4)
            ax.scatter([i], [mean], s=0)

    ax.set_xticks(positions)
    ax.set_xticklabels(CELL_LINE_ORDER)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ymax = np.nanmax(pd.to_numeric(df[value_col], errors='coerce').values)
    ymin = np.nanmin(pd.to_numeric(df[value_col], errors='coerce').values)
    yr = ymax - ymin if ymax > ymin else 1.0
    base = ymax + 0.08*yr
    h = 0.03*yr if yr > 0 else 0.5
    for k, row in pairwise.iterrows():
        i = CELL_LINE_ORDER.index(row['group1']) + 1
        j = CELL_LINE_ORDER.index(row['group2']) + 1
        add_stat_bracket(ax, i, j, base + k*0.09*yr, fmt_p(row['p_holm']), h)
    ax.text(0.03, 0.97, f"Kruskal–Wallis H = {kw['H']:.3f}\n{fmt_p(kw['p'])}", transform=ax.transAxes,
            va='top', ha='left', fontsize=9,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.25'))

    # Legends
    exp_handles = [Line2D([0], [0], marker='o', linestyle='', markerfacecolor=REDS[i], markeredgecolor='none', label=EXPERIMENT_ORDER[i], markersize=7)
                   for i in range(len(EXPERIMENT_ORDER))]
    leg1 = ax.legend(handles=exp_handles, title='Experiment', loc='upper right', frameon=False)
    ax.add_artist(leg1)
    if use_shapes:
        shape_handles = [Line2D([0], [0], marker=SHAPES[g], linestyle='', markerfacecolor='black', markeredgecolor='none', label=g, markersize=7)
                         for g in CELL_LINE_ORDER]
        ax.legend(handles=shape_handles, title='Cell line', loc='upper right', bbox_to_anchor=(1.0, 0.72), frameon=False)

    save_all_formats(fig, out_base)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Panel A Kruskal–Wallis figures + B/C descriptive stats')
    ap.add_argument('--batch-dir', required=True, help='Batch directory (e.g. ~/neuron_ims_batch50_maskfix) or stage_root path')
    ap.add_argument('--strict', action='store_true', help='Raise on missing columns')
    args = ap.parse_args()

    batch_dir = Path(args.batch_dir).expanduser().resolve()
    stage_root = load_stage_root(batch_dir)
    metrics_csv = discover_metrics_csv(stage_root)
    figdir = ensure_dir(stage_root / '_PV_selected_figures')
    outdir = ensure_dir(figdir / 'panelA_kw')

    df = pd.read_csv(metrics_csv)
    cell_col = require_col(df, ['cell_line'], 'cell line')
    exp_col = require_col(df, ['experiment_id', 'experiment', 'batch', 'exp_id'], 'experiment id')
    df = df.rename(columns={cell_col: 'cell_line', exp_col: 'experiment_id'})
    df['cell_line'] = df['cell_line'].astype(str)
    df['experiment_id'] = df['experiment_id'].astype(str)

    # Metrics
    density, density_source = maybe_compute_density(df)
    nucfrac, frac_source = maybe_compute_nuclear_fraction(df)
    df['A_density_kw'] = density
    df['A_nuclear_fraction_kw'] = nucfrac

    # panel A stats
    panelA_rows = []
    pairwise_paths = []
    for metric, title, ylabel, base_name in [
        ('A_density_kw', 'PV-positive cell density', 'PV+ cells / mm²', 'A_PV_density_kw_publication'),
        ('A_nuclear_fraction_kw', 'PV-positive Hoechst-colocalized cells', 'PV+ / Hoechst+ cells (%)', 'A_PV_nuclear_fraction_kw_publication'),
    ]:
        sub = df[['cell_line','experiment_id', metric]].copy()
        sub = sub.loc[sub['cell_line'].isin(CELL_LINE_ORDER)].dropna()
        kw = overall_kw(sub, 'cell_line', metric, CELL_LINE_ORDER)
        pair = pairwise_mwu(sub, 'cell_line', metric, CELL_LINE_ORDER)
        pair['metric'] = metric
        pair['overall_H'] = kw['H']
        pair['overall_p'] = kw['p']
        pair.to_csv(outdir / f'{base_name}_pairwise_Holm.csv', index=False)
        pairwise_paths.append(outdir / f'{base_name}_pairwise_Holm.csv')
        panel_plot(sub, 'cell_line', metric, ylabel, title, outdir / base_name, pair, kw, use_shapes=True)

        descs = []
        for g in CELL_LINE_ORDER:
            vals = sub.loc[sub['cell_line']==g, metric]
            st = desc_stats(vals)
            st.update({'metric': metric, 'grouping': 'cell_line', 'group': g})
            descs.append(st)
        panelA_rows.append(pd.DataFrame(descs))
        panelA_rows.append(pd.DataFrame([{ 'metric': metric, 'grouping': 'overall', 'group': 'Kruskal-Wallis', 'H': kw['H'], 'p': kw['p'], 'source_column_or_formula': density_source if metric=='A_density_kw' else frac_source }]))

    panelA_desc = pd.concat(panelA_rows, ignore_index=True, sort=False)
    panelA_desc.to_csv(outdir / 'panelA_kw_descriptives_and_tests.csv', index=False)

    # BC descriptive only
    bcmap = get_BC_metrics(df)
    bc_rows = []
    for label, col in bcmap.items():
        sub = df[['cell_line', 'experiment_id', col]].copy().rename(columns={col: 'value'})
        sub = sub.loc[sub['value'].notna() & sub['cell_line'].isin(CELL_LINE_ORDER)]
        for g in CELL_LINE_ORDER:
            vals = sub.loc[sub['cell_line']==g, 'value']
            if len(vals) == 0:
                continue
            st = desc_stats(vals)
            st.update({'metric': label, 'grouping': 'cell_line', 'group': g})
            bc_rows.append(st)
        # overall per metric
        st = desc_stats(sub['value'])
        st.update({'metric': label, 'grouping': 'overall', 'group': 'all'})
        bc_rows.append(st)
    bc_df = pd.DataFrame(bc_rows)
    bc_df.to_csv(outdir / 'BC_descriptive_mean_SD.csv', index=False)

    # md report
    md = []
    md.append('# Panel A Kruskal–Wallis results and B/C descriptive statistics\n')
    md.append(f'- Input metrics CSV: `{metrics_csv}`')
    md.append(f'- Density metric source: `{density_source}`')
    md.append(f'- Nuclear PV fraction source: `{frac_source}`\n')
    md.append('## Panel A tests\n')
    for metric, name in [('A_density_kw', 'PV-positive cell density'), ('A_nuclear_fraction_kw', 'PV-positive Hoechst-colocalized cells (%)')]:
        sub = df[['cell_line', metric]].dropna()
        kw = overall_kw(sub, 'cell_line', metric, CELL_LINE_ORDER)
        md.append(f'### {name}')
        md.append(f'- Kruskal–Wallis: H = {kw["H"]:.4f}, {fmt_p(kw["p"])}')
        pair = pairwise_mwu(sub, 'cell_line', metric, CELL_LINE_ORDER)
        for _, r in pair.iterrows():
            md.append(f'  - {r.group1} vs {r.group2}: U = {r.U:.3f}, raw p = {r.p_raw:.4g}, Holm-adjusted {fmt_p(r.p_holm)}')
        md.append('')
    md.append('## B/C descriptive statistics')
    md.append('See `BC_descriptive_mean_SD.csv` for mean, SD, variance, SEM, median, IQR, and n by cell line.\n')
    (outdir / 'panelA_kw_results_and_BC_descriptives.md').write_text('\n'.join(md))

    print(f'Stage root: {stage_root}')
    print(f'Input metrics CSV: {metrics_csv}')
    print(f'Output directory: {outdir}')
    print('Wrote:')
    for p in [
        outdir / 'A_PV_density_kw_publication.svg',
        outdir / 'A_PV_nuclear_fraction_kw_publication.svg',
        outdir / 'panelA_kw_descriptives_and_tests.csv',
        outdir / 'BC_descriptive_mean_SD.csv',
        outdir / 'panelA_kw_results_and_BC_descriptives.md',
    ]:
        print(f'  {p}')

if __name__ == '__main__':
    main()
