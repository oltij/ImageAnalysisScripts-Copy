"""Synthetic smoke tests. Never use simulated images as biological results."""
import importlib.util
import io
import json
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import tifffile
from matplotlib.backends.backend_pdf import PdfPages

SCRIPT = Path(__file__).resolve().parents[1] / 'build_per_organoid_reports.py'
spec = importlib.util.spec_from_file_location('organoid_builder', SCRIPT)
mod = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_synthetic_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        analysis = root / 'analysis'
        shape = (55, 65)
        image = (np.arange(shape[0]*shape[1], dtype=np.uint16).reshape(shape)%309)*10
        mask = np.zeros(shape, np.uint8)
        mask[3:51, 4:61] = 1
        for channel in ('Hoechst','GFP','LHX6','PV'):
            p=analysis/'00_max_projections'/f'{channel}.tif'
            p.parent.mkdir(parents=True, exist_ok=True)
            tifffile.imwrite(p, image)
            if channel != 'Hoechst':
                pp=analysis/'02_registration'/channel/f'aligned_{channel}_max.tif'
                pp.parent.mkdir(parents=True, exist_ok=True)
                tifffile.imwrite(pp, image)
            else:
                p=analysis/'03_aligned_segmentation'/channel/'OrganoidMask.tiff'
                p.parent.mkdir(parents=True, exist_ok=True)
                tifffile.imwrite(p,mask)
        for channel in ('GFP','LHX6','PV'):
            folder = analysis/'06_final_shape_filtered'/channel
            folder.mkdir(parents=True)
            pixels=[]
            metrics=[]
            for oid in range(1,21):
                x=8+(oid%10)*5; y=8+(oid//10)*9
                for dx in range(3):
                    for dy in range(3):
                        pixels.append(dict(ROI=oid,X=x+dx,Y=y+dy))
                metrics.append(dict(ObjectID=oid,Area_um2=oid*.7,**{
                    'Mean intensity': float(oid*10),'Integrated intensity':float(oid*90),
                    'Circularity':.1+.02*oid,'Solidity':.5+.01*oid,'Eccentricity':.01+.03*oid}))
            pd.DataFrame(pixels).to_csv(folder/f'{channel}_ROI_pixels_FILTERED.csv',index=False)
            pd.DataFrame(metrics).to_csv(folder/f'{channel}_kept_ROI_metrics.csv',index=False)
        h5p=analysis/'08_marker_colocalization'/'LHX6_vs_PV'/'analysis.h5'
        h5p.parent.mkdir(parents=True)
        ap=[];bp=[]
        for oid in (1,2,3):
            for dx in range(3):
                for dy in range(3):
                    ap.append(dict(ObjectID=oid,X=10+oid*8+dx,Y=12+dy))
                    bp.append(dict(ObjectID=oid+10,X=11+oid*8+dx,Y=12+dy))
        overlaps=pd.DataFrame([
            dict(A_ObjectID=1,B_ObjectID=11,A_Percent_Inside_B=66.67,B_Percent_Inside_A=66.67,Shared_Area=6),
            dict(A_ObjectID=2,B_ObjectID=12,A_Percent_Inside_B=66.67,B_Percent_Inside_A=66.67,Shared_Area=6),
            dict(A_ObjectID=3,B_ObjectID=13,A_Percent_Inside_B=66.67,B_Percent_Inside_A=66.67,Shared_Area=6),
        ])
        matching=overlaps.iloc[:2].copy()
        meta={'objects':{'object_a_name':'LHX6','object_b_name':'PV'},'matching':{'threshold_percent':50}}
        with h5py.File(h5p,'w') as h:
            h.create_dataset('metadata/json',data=np.bytes_(json.dumps(meta)))
            h.create_dataset('images/object_a_original',data=image)
            h.create_dataset('images/object_b_original',data=image*2)
            for key,frame in {'object_a_pixels':pd.DataFrame(ap),'object_b_pixels':pd.DataFrame(bp),'stable_matching':matching,'overlap_results':overlaps}.items():
                h.create_dataset('tables/'+key,data=np.frombuffer(frame.to_csv(index=False).encode(),dtype='uint8'))
        p=mod.load_pair(h5p)
        assert len(p.matching)==2
        intens=mod.intensity_rows(p,0,'151-2')
        assert len(intens)==2 and (intens['overlap_area_pixels']==6).all()
        assert intens['A_roi_integrated_intensity_native_sum'].gt(0).all()
        assert intens['B_overlap_integrated_intensity_native_sum'].gt(0).all()
        allowed=mod.stable_ids_by_channel([h5p])
        assert allowed['LHX6']=={1,2} and allowed['PV']=={11,12}
        q=mod.quartile_representatives(pd.DataFrame([{'ObjectID':i,'Area_um2':i} for i in range(1,9)]),'Area_um2',['ObjectID'])
        assert [int(row['ObjectID']) for _,row,_,_ in q]==[2,4,6,8]
        # Genuine final stable matched pair example (exclude unmatched 3/13).
        # ROI area/intensities and plots refer to the same A/B identifiers.
        pdf=root/'test_synthetic_report.pdf'
        with PdfPages(pdf) as out:
            mod.whole_page(out,analysis,[])
            rows=mod.roi_quartile_pages(out,analysis,'LHX6',['integrated_intensity'],{1,2},0,'151-3',[])
            vals=mod.matched_pair_quartile_page(out,p,0,'151-3')
        assert len(rows)==2
        assert len(vals)==2 and set(vals['A_ObjectID'])=={1,2}
        assert set(vals['B_ObjectID'])=={11,12}
        import fitz
        doc=fitz.open(pdf)
        txt='\n'.join(page.get_text() for page in doc)
        assert doc.page_count==3
        assert 'FINAL stable-matched' in txt
        assert 'Integrated fluorescence' in txt
        assert 'ROI #1' in txt or 'ROI #2' in txt
        assert 'Q1' in txt and 'Q2' in txt
        doc.close()
        # Integration: invoke CLI main() with a small stub for upstream field
        # selection while keeping the real TIFF/ROI/HDF5 readers and PDF writer.
        stub = root / 'plot_igf_exp17_ventral_with_gfp_plots.py'
        record = {
            'index': 0, 'sample': '151-3', 'block_id': 151, 'source_field':
                'Exp17_151-3_Set3_DAPI_GFP_LHX6_PV_Field_1', 'output_dir': str(analysis),
            'cell_line': 'H1', 'condition': '+IGF', 'anti_gfp_status': 'No anti-GFP (-3)',
            'pv_positive_density_per_mm2': 50.0, 'gfp_density_per_mm2': 100.0,
        }
        stub.write_text(
            'CELL_LINE_ORDER=("H1","CL23","CL32")\n'
            'CONDITIONS=("No IGF","+IGF")\n'
            'ANTI_GFP_ORDER=("Anti-GFP stained (-2)","No anti-GFP (-3)")\n'
            'METRIC_LEGENDS={}\n'
            'PLOT_METRICS=(("pv_positive_density_per_mm2","All PV+","Cells/mm2","pv"),)\n'
            'GFP_PLOT_METRICS=(("gfp_density_per_mm2","All GFP+","ROIs/mm2","gfp"),)\n'
            'def load_sample_info(path): return {}\n'
            'def collect_fields(original, other, selection, info):\n'
            f'    return [{record!r}], [], []\n'
        )
        orig_here = mod.HERE
        mod.HERE = root
        try:
            code = mod.main(['--h1-workspace', str(root), '--remaining-workspace', str(root),
                             '--out', str(root/'reports'), '--only-sample','151-3',
                             '--max-fields','1', '--max-examples','1',
                             '--roi-metrics','area','--compact'])
        finally:
            mod.HERE = orig_here
        assert code == 0
        produced = root/'reports'/'151-3'/'Exp17_151-3_Set3_DAPI_GFP_LHX6_PV_Field_1_organoid_report_quartiles.pdf'
        assert produced.is_file() and produced.stat().st_size>10000
        intensity_csv = root/'reports'/'151-3'/'Exp17_151-3_Set3_DAPI_GFP_LHX6_PV_Field_1_displayed_stable_matched_pair_quartiles.csv'
        assert len(pd.read_csv(intensity_csv)) == 2 and set(pd.read_csv(intensity_csv)['A_ObjectID'])=={1,2}
        assert (root/'reports'/'report_manifest.csv').is_file()
        print(f'Synthetic test passed: {pdf.stat().st_size} bytes')
