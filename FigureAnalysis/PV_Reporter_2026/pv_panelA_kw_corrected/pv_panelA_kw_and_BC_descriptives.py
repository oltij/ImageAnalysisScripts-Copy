#!/usr/bin/env python3
"""Publication figures: KW test for nuclear-associated PV density and PV area; B/C descriptive only.

Designed for the exact columns emitted by the 2026-10-09 selected-acquisition batch.
Pure downstream analysis: no image processing, segmentation, or rematching.
"""
from __future__ import annotations

import argparse
import io
import json
import math
import warnings
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

LINES = ['H1', 'H1 CL32', 'H1 CL23']
EXPERIMENTS = ['17', '24', 'RA02']
EXP_COLORS = {'17':'#F7A1A3', '24':'#DF4C50', 'RA02':'#801C31'}
MARKERS = {'H1':'o', 'H1 CL32':'s', 'H1 CL23':'^'}
# The acquisition's experiment field is a numbered block, not the experimental batch.
BLOCK_TO_EXP = {
    '165':'24', '166':'24', '167':'24', '168':'24',
    '169':'RA02', '170':'RA02', '171':'RA02', '172':'RA02', '173':'RA02', '174':'RA02',
    '175':'17', '176':'17', '177':'17', '178':'17',
}
PLOT_METRICS = [
    ('pv_nuclear_positive_density_mm2', 'PV-positive, Hoechst-associated cell density',
     'PV-positive cells / mm²', 'A_density_KruskalWallis'),
    ('pv_nuclear_assoc_area_pct', 'PV-positive area from Hoechst-associated PV objects',
     'Hoechst-associated PV area / organoid area (%)', 'A_nuclear_associated_area_KruskalWallis'),
]


def batch_root(path: Path) -> Path:
    prepared = path / 'prepared.json'
    if prepared.is_file():
        obj = json.loads(prepared.read_text())
        if 'stage_root' not in obj:
            raise KeyError(f'{prepared} lacks stage_root')
        root = Path(obj['stage_root']).expanduser()
        if not root.is_absolute():
            root = path / root
        return root.resolve()
    return path.resolve()


def describe(s: pd.Series) -> dict:
    a = pd.to_numeric(s, errors='coerce').dropna().astype(float).to_numpy()
    if len(a) == 0:
        return {'n':0, 'mean':np.nan, 'sd':np.nan, 'median':np.nan, 'q1':np.nan,
                'q3':np.nan, 'sem':np.nan, 'min':np.nan, 'max':np.nan}
    return {'n':int(len(a)), 'mean':float(np.mean(a)),
            'sd':float(np.std(a,ddof=1)) if len(a)>1 else np.nan,
            'median':float(np.median(a)), 'q1':float(np.quantile(a,.25)),
            'q3':float(np.quantile(a,.75)),
            'sem':float(np.std(a,ddof=1)/np.sqrt(len(a))) if len(a)>1 else np.nan,
            'min':float(np.min(a)), 'max':float(np.max(a))}


def fmtp(p:float) -> str:
    if p < .0001: return 'p < 0.0001'
    return f'p = {p:.4f}'


def hc(pvals:list[float]) -> list[float]:
    order = np.argsort(pvals)
    adjusted = np.zeros(len(pvals))
    running = 0.0
    for i, j in enumerate(order):
        running = max(running, (len(pvals)-i)*pvals[j])
        adjusted[j] = min(1.0, running)
    return adjusted.tolist()


def h5_pv_assoc_area(row:pd.Series) -> dict:
    """Take unique XY pixels for all PV objects matched to Hoechst.

    The H5 'filtered_object_a_pixels' table contains whole, segmented PV ROI objects
    *selected* by nucleus colocalization, not just the PV/nucleus overlapping pixels.
    Percent area = union(PV-associated-object pixels) / Hoechst tissue mask area * 100.
    """
    h5path = Path(str(row['source_pv_nuclear_h5']))
    if not h5path.is_file():
        raise FileNotFoundError(f'Nuclear PV HDF5 missing: {h5path} ({row["acquisition_label"]})')
    with h5py.File(h5path,'r') as h5:
        key = 'tables/filtered_object_a_pixels'
        if key not in h5:
            raise KeyError(f'{h5path}: expected {key}; keys: {list(h5.keys())}')
        pixels = pd.read_csv(io.BytesIO(h5[key][()].tobytes()))
        # Cross-check that this file is indeed PV (object A) vs Hoechst (object B).
        meta_key = 'metadata/json'
        if meta_key in h5:
            raw = h5[meta_key][()]
            if isinstance(raw, bytes): obj = json.loads(raw.decode('utf-8'))
            else: obj = json.loads(bytes(raw).decode('utf-8'))
            a = str(obj.get('objects',{}).get('object_a_name','')).lower()
            b = str(obj.get('objects',{}).get('object_b_name','')).lower()
            if a and a != 'pv': raise ValueError(f'{h5path}: object A should be PV, got {a}')
            if b and b not in ('hoechst','dapi'): raise ValueError(f'{h5path}: object B should be Hoechst, got {b}')

    if not {'CellID','X','Y'}.issubset(pixels.columns):
        raise KeyError(f'{h5path}: filtered_object_a_pixels columns should contain CellID, X, Y; got {list(pixels.columns)}')
    unique_cells = pixels['CellID'].dropna().nunique()
    expected = int(row['pv_nuclear_positive_cells'])
    if unique_cells != expected:
        raise ValueError(f'{row["acquisition_label"]}: HDF5 matched PV cell count ({unique_cells}) != CSV count ({expected}); refusing output')
    unique_xy = len(pixels.drop_duplicates(subset=['X','Y']))
    total_pixels = int(row['organoid_area_pixels'])
    if total_pixels <= 0: raise ValueError(f'Invalid organoid_area_pixels for {row["acquisition_label"]}')
    if unique_xy > total_pixels:
        raise ValueError(f'{row["acquisition_label"]}: nuclear-associated PV pixels {unique_xy} > organoid pixels {total_pixels}')
    orig_union = int(row['pv_positive_union_pixels'])
    if unique_xy > orig_union:
        raise ValueError(f'{row["acquisition_label"]}: nuclear-associated PV pixels {unique_xy} > all PV-positive union pixels {orig_union}')
    return {'pv_nuclear_assoc_union_pixels':unique_xy,
            'pv_nuclear_assoc_area_pct':100*unique_xy/total_pixels,
            'pv_nuclear_h5_matched_cells':int(unique_cells)}


def stats_for(df:pd.DataFrame, value_col:str):
    arrays=[df.loc[df.cell_line==g, value_col].dropna().to_numpy(dtype=float) for g in LINES]
    if any(len(a)<2 for a in arrays):
        raise ValueError(f'Insufficient data for KW {value_col}: {[len(a) for a in arrays]}')
    h,p = kruskal(*arrays)
    pairs=[]
    for i in range(3):
        for j in range(i+1,3):
            u,raw = mannwhitneyu(arrays[i],arrays[j],alternative='two-sided',method='auto')
            pairs.append({'metric':value_col,'group1':LINES[i],'group2':LINES[j],
                          'n1':len(arrays[i]),'n2':len(arrays[j]),'U':float(u),'p_raw':float(raw)})
    adjusted=hc([p['p_raw'] for p in pairs])
    for r,a in zip(pairs,adjusted): r['p_holm']=a
    return {'metric':value_col,'H':float(h),'df':2,'p':float(p),'n_total':sum(map(len,arrays))},pd.DataFrame(pairs)


def draw_panel(df:pd.DataFrame,metric:str,title:str,ylabel:str,base:Path,kw:dict,pairs:pd.DataFrame):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none', 'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(7.1,6.4))
    arr=[df.loc[df.cell_line==g,metric].dropna().to_numpy(dtype=float) for g in LINES]
    b=ax.boxplot(arr, positions=[1,2,3],widths=.49,whis=1.5,showfliers=False,patch_artist=True,
                 medianprops={'color':'#222222','linewidth':1.6},
                 boxprops={'edgecolor':'#222222','linewidth':1.2},
                 whiskerprops={'color':'#222222','linewidth':1.2},
                 capprops={'color':'#222222','linewidth':1.2})
    for patch in b['boxes']:patch.set_facecolor('white')
    for j,g in enumerate(LINES):
        sub=df[df.cell_line==g].reset_index(drop=True)
        # Fixed, deterministic offsets (avoid moving points on reruns).
        offsets=np.random.default_rng(123+j).uniform(-.16,.16,size=len(sub))
        for k,(_,r) in enumerate(sub.iterrows()):
            ax.scatter(j+1+offsets[k],r[metric],s=43, marker=MARKERS[g],
                       c=EXP_COLORS[r['experiment_id']],edgecolors='none',alpha=.92,zorder=5)
        mean=float(np.mean(arr[j])); sd=float(np.std(arr[j],ddof=1))
        ax.errorbar(j+1,mean,yerr=sd,fmt='_',color='#000000',capsize=5,markersize=20,
                    capthick=2.0,elinewidth=2.0,zorder=7)

    ymax=float(np.max(np.concatenate(arr))); ymin=float(np.min(np.concatenate(arr)))
    span=max(ymax-ymin, 0.15*max(abs(ymax),1))
    # Brackets at staggered heights; avoid legend overlay via bbox_to_anchor.
    heights=[ymax+.10*span, ymax+.21*span, ymax+.32*span]
    for idx,r in pairs.reset_index(drop=True).iterrows():
        i=LINES.index(r.group1)+1; j=LINES.index(r.group2)+1
        y=heights[idx];rise=.025*span
        ax.plot([i,i,j,j],[y,y+rise,y+rise,y],lw=1.15,color='black',clip_on=False)
        ax.text((i+j)/2,y+rise+.012*span,fmtp(r.p_holm),ha='center',va='bottom',fontsize=8.8)
    ax.set_ylim(min(0,ymin-.08*span),ymax+.48*span)
    ax.set_xlim(.45,3.55)
    ax.set_xticks([1,2,3],LINES)
    ax.set_ylabel(ylabel)
    ax.set_title(title,fontsize=12,weight='bold',pad=14)
    ax.spines[['top','right']].set_visible(False)
    ax.text(.02,.97,f'Kruskal–Wallis: H(2) = {kw["H"]:.2f}; {fmtp(kw["p"])}',
            transform=ax.transAxes,va='top',ha='left',fontsize=9.2)
    # Put experiment legend outside plot, preserving consistent colors.
    handles=[Line2D([0],[0],marker='o',linestyle='None',markersize=7,color=EXP_COLORS[e],label=f'Exp. {e}') for e in EXPERIMENTS]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.085),ncol=3,frameon=False,fontsize=9)
    fig.text(.5,.017,'Boxes: median and IQR; whiskers: 1.5×IQR; black marks: mean ± SD.\nPairwise: two-sided Mann–Whitney tests with Holm correction.',ha='center',fontsize=8.5)
    fig.subplots_adjust(bottom=.23)
    for ext in ['svg','pdf','png']:
        fig.savefig(base.with_suffix('.'+ext),dpi=600 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--batch-dir',type=Path,required=True)
    ap.add_argument('--strict',action='store_true',help='Require 35 observations and expected B/C samples')
    a=ap.parse_args()
    root=batch_root(a.batch_dir.expanduser().resolve())
    csv=root/'_PV_selected_figures'/'per_acquisition_metrics.csv'
    if not csv.is_file():raise FileNotFoundError(f'Missing expected metrics CSV: {csv}')
    data=pd.read_csv(csv)
    required=['experiment','acquisition_label','cell_line','enhancer','organoid_area_pixels',
              'pv_positive_union_pixels','pv_nuclear_positive_cells','pv_nuclear_positive_density_mm2',
              'source_pv_nuclear_h5','specificity_percent','coverage_percent']
    missing=[c for c in required if c not in data.columns]
    if missing: raise ValueError(f'Missing expected columns: {missing}')
    if data.acquisition_label.duplicated().any():raise ValueError('Duplicate acquisition_label; expected one record per selected acquisition')
    data['experiment_id']=data['experiment'].astype(str).str.extract(r'^(\d{3})',expand=False).map(BLOCK_TO_EXP)
    if data.experiment_id.isna().any():
        invalid=data.loc[data.experiment_id.isna(),'experiment'].unique().tolist()
        raise ValueError(f'Unknown experiment block IDs: {invalid}')
    if (~data.cell_line.isin(LINES)).any():raise ValueError(f'Unexpected cell_line: {data.cell_line.unique().tolist()}')
    if (~data.enhancer.isin(['BiVe3','BiVe4'])).any():raise ValueError(f'Unexpected enhancer: {data.enhancer.unique().tolist()}')
    if a.strict and len(data)!=35: raise ValueError(f'Expected 35 selected rows, got {len(data)}')
    observed=data.enhancer.value_counts().to_dict()
    if a.strict and (observed.get('BiVe3')!=17 or observed.get('BiVe4')!=18):
        raise ValueError(f'Expected BiVe3 n=17 and BiVe4 n=18; observed {observed}')
    if data[['pv_nuclear_positive_density_mm2','organoid_area_pixels']].isna().any().any():
        raise ValueError('Missing density or organoid area values')
    # Check the existing density calculation for consistency with counts and pixel calibration.
    if not (data['pv_nuclear_positive_density_mm2']>0).all():
        raise ValueError('Unexpected zero/negative nuclear PV cell density; review input')
    # Compute whole-object area only from the PV-then-Hoechst association result HDF5.
    # Do not substitute raw PV positive area or PV+ nuclei / Hoechst nuclei fraction.
    areas=[]
    for k,r in data.iterrows():
        derived=h5_pv_assoc_area(r)
        areas.append(derived)
        print(f'[{len(areas):02d}/{len(data)}] {r["acquisition_label"]}: nuclear PV cells={derived["pv_nuclear_h5_matched_cells"]}, associated area={derived["pv_nuclear_assoc_area_pct"]:.3f}%',flush=True)
    for key in areas[0]: data[key]=[r[key] for r in areas]
    outdir=root/'_PV_selected_figures'/'panelA_kw'
    outdir.mkdir(parents=True,exist_ok=True)
    data.to_csv(outdir/'per_acquisition_metrics_with_nuclear_associated_PV_area.csv',index=False)

    omnibus=[]; pair_all=[]; desc=[]
    for metric,title,ylab,stem in PLOT_METRICS:
        kw,pairs=stats_for(data,metric)
        omnibus.append(kw)
        pair_all.append(pairs)
        pairs.to_csv(outdir/f'{stem}_pairwise_Holm.csv',index=False)
        for group in LINES:
            desc.append({'metric':metric,'grouping':'cell_line','group':group,**describe(data.loc[data.cell_line==group,metric])})
        draw_panel(data,metric,title,ylab,outdir/stem,kw,pairs)
    kw_df=pd.DataFrame(omnibus)
    kw_df.to_csv(outdir/'Panel_A_KruskalWallis_omnibus.csv',index=False)
    pd.concat(pair_all,ignore_index=True).to_csv(outdir/'Panel_A_pairwise_MannWhitney_Holm.csv',index=False)
    pd.DataFrame(desc).to_csv(outdir/'Panel_A_descriptive_mean_SD.csv',index=False)

    bc=[]
    for enhancer,label in [('BiVe3','B'),('BiVe4','C')]:
        subset=data[data.enhancer==enhancer]
        for metric in ['specificity_percent','coverage_percent']:
            for group in LINES:
                vals=subset.loc[subset.cell_line==group,metric]
                if vals.empty:raise ValueError(f'Missing {enhancer}/{metric}/{group}')
                bc.append({'panel':label,'enhancer':enhancer,'metric':metric,
                           'cell_line':group, **describe(vals)})
            bc.append({'panel':label,'enhancer':enhancer,'metric':metric,
                       'cell_line':'ALL',**describe(subset[metric])})
    bc_df=pd.DataFrame(bc)
    bc_df.to_csv(outdir/'BC_descriptive_mean_SD.csv',index=False)
    lines=['# PV figure statistics','',
           'Panel A: Kruskal–Wallis omnibus test (3 cell lines) and 3 two-sided Mann–Whitney tests per endpoint (Holm correction).',
           'BiVe3 and BiVe4 measurements pooled for Panel A. B/C: descriptive only.',
           'Each CSV row is an acquisition; results assume acquisitions are independent organoids.',
           'Kruskal–Wallis does NOT adjust for experimental batch. Batch composition is unbalanced.',
           '', '## Panel A results','']
    for _,r in kw_df.iterrows():
        lines.append(f'### {r.metric}')
        lines.append(f'Kruskal–Wallis H(2) = {r.H:.4f}; {fmtp(r.p)}; n = {int(r.n_total)}')
        for _,q in pd.concat(pair_all,ignore_index=True).query('metric == @r.metric').iterrows():
            lines.append(f'- {q.group1} vs {q.group2}: U = {q.U:.3f}, unadjusted {fmtp(q.p_raw)}, Holm-adjusted {fmtp(q.p_holm)}')
        lines.append('')
    lines+=['## Area definition','',
            'PV-associated area (%) = 100 × unique (X,Y) pixels among full, segmented PV objects matched to Hoechst / Hoechst-defined organoid mask pixels.',
            'The PV objects come from the HDF5 dataset tables/filtered_object_a_pixels in each source_pv_nuclear_h5.',
            'This is NOT raw all-PV area, and it is NOT the fraction of Hoechst nuclei that are PV positive.',
            '', '## B/C mean ± SD', '',bc_df.to_markdown(index=False,floatfmt='.3f'),'']
    (outdir/'panelA_kw_results_and_BC_descriptives.md').write_text('\n'.join(lines))
    print('\n'+'='*72)
    print('COMPLETE. Results:',outdir)
    print(kw_df.to_string(index=False))
    print('\nB/C descriptive statistics: BC_descriptive_mean_SD.csv')

if __name__=='__main__':main()
