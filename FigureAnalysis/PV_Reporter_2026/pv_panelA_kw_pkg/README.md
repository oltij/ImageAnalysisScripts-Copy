# Panel A Kruskal–Wallis + B/C descriptive statistics

This package generates:

1. **Panel A density figure** with:
   - boxplots (1.5×IQR whiskers)
   - individual organoid points
   - experiment colors (17, 24, RA02)
   - cell-line shapes
   - mean ± SD black error bars
   - **Kruskal–Wallis** omnibus p-value
   - **Holm-adjusted pairwise Mann–Whitney** p-values

2. **A second Panel A figure** based on **Hoechst-colocalized PV-positive cells**, plotted as:
   - **PV+ / Hoechst+ cells (%)**
   - same styling and Kruskal–Wallis annotations

3. **B/C descriptive statistics only**:
   - mean
   - SD
   - variance
   - SEM
   - median
   - IQR
   - n

## Files produced

Inside:

`<stage_root>/_PV_selected_figures/panelA_kw/`

The script writes:

- `A_PV_density_kw_publication.svg/.pdf/.png`
- `A_PV_nuclear_fraction_kw_publication.svg/.pdf/.png`
- `A_PV_density_kw_publication_pairwise_Holm.csv`
- `A_PV_nuclear_fraction_kw_publication_pairwise_Holm.csv`
- `panelA_kw_descriptives_and_tests.csv`
- `BC_descriptive_mean_SD.csv`
- `panelA_kw_results_and_BC_descriptives.md`

## Great Lakes install

From your Mac:

```bash
scp ~/Downloads/pv_panelA_kw_and_BC_descriptives.zip \
  oltij@greatlakes-xfer.arc-ts.umich.edu:~/
```

On Great Lakes:

```bash
conda activate ims-mosaic
cd ~
unzip -oq pv_panelA_kw_and_BC_descriptives.zip -d pv_panelA_kw_pkg
cd pv_panelA_kw_pkg

python pv_panelA_kw_and_BC_descriptives.py \
  --batch-dir ~/neuron_ims_batch50_maskfix \
  --strict
```

If your `--batch-dir` contains `prepared.json`, the script will use `stage_root` from that file.
Otherwise it treats the directory as the stage root directly.

## Download results back to your Mac

```bash
mkdir -p "$HOME/Desktop/PV_KW_Final"

rsync -avh --progress \
  'oltij@greatlakes-xfer.arc-ts.umich.edu:/gpfs/accounts/parent_root/parent0/oltij/Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/_PV_selected_figures/panelA_kw/' \
  "$HOME/Desktop/PV_KW_Final/"
```
