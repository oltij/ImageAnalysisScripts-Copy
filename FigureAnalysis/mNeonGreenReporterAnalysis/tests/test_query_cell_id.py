"""Synthetic query tests: scoped Stage-07 CellIDs and independent Stage-08/09 outcomes."""
from __future__ import annotations
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "query_cell_id.py"


def store(h5, key, frame):
    h5.create_dataset(f"tables/{key}", data=np.frombuffer(frame.to_csv(index=False).encode(), dtype="uint8"))


def source(path, marker, source_pixels, mappings, nuclei):
    path.parent.mkdir(parents=True, exist_ok=True)
    pixels = pd.DataFrame([{"ObjectID": oid, "X": x, "Y": y, "ROI": oid} for oid, xy in source_pixels.items() for x,y in xy])
    mp = pd.DataFrame(mappings, columns=["CellID", "A_ObjectID", "B_ObjectID"])
    filt = pixels.merge(mp, left_on="ObjectID", right_on="A_ObjectID")[["CellID", "ROI", "X", "Y"]]
    npix = pd.DataFrame([{"ObjectID": oid, "X":x,"Y":y} for oid,xy in nuclei.items() for x,y in xy])
    with h5py.File(path, "w") as h5:
        h5.create_dataset("metadata/json", data=np.bytes_(json.dumps({"objects":{"object_a_name":marker,"object_b_name":"Hoechst"}})))
        h5.create_dataset("images/object_a_original", data=np.arange(100,dtype="uint16").reshape(10,10))
        h5.create_dataset("images/object_b_original", data=np.arange(100,dtype="uint16").reshape(10,10))
        store(h5,"cell_id_mapping",mp)
        store(h5,"object_a_pixels",pixels)
        store(h5,"object_b_pixels",npix)
        store(h5,"filtered_object_a_pixels",filt)


def pair(path, stage, a_pixels, b_pixels, mapping_a, mapping_b, matches):
    path.parent.mkdir(parents=True,exist_ok=True)
    def pix(d,m):
        return pd.DataFrame([{"ObjectID":oid,"ROI":oid,"CellID":m[oid],"X":x,"Y":y} for oid,xy in d.items() for x,y in xy])
    a=pix(a_pixels,mapping_a);b=pix(b_pixels,mapping_b)
    overlaps = pd.DataFrame([[10,20,2,66.66,100.],[10,21,3,100.,100.]],
                            columns=["A_ObjectID","B_ObjectID","Shared_Area","A_Percent_Inside_B","B_Percent_Inside_A"])
    match = pd.DataFrame(matches, columns=["A_ObjectID","B_ObjectID"])
    out_mapping=match.copy();out_mapping.insert(0,"CellID",np.arange(1,len(match)+1))
    with h5py.File(path,"w") as h5:
        h5.create_dataset("metadata/json",data=np.bytes_(json.dumps({"objects":{"object_a_name":"LHX6","object_b_name":"PV"},
                                                                      "matching":{"threshold_percent":50},
                                                                      **({"analysis_variant":"nucleus_constrained_marker_marker"} if stage==9 else {})})))
        h5.create_dataset("images/object_a_original",data=np.arange(100,dtype="uint16").reshape(10,10))
        h5.create_dataset("images/object_b_original",data=np.arange(100,dtype="uint16").reshape(10,10))
        store(h5,"object_a_pixels",a)
        store(h5,"object_b_pixels",b)
        store(h5,"overlap_results",overlaps)
        store(h5,"stable_matching",match)
        store(h5,"cell_id_mapping",out_mapping)
        if stage==9:
            store(h5,"nuclear_constraint_audit",pd.DataFrame([
                {"A_ObjectID":10,"B_ObjectID":20,"Final_Constrained_Eligibility":True,"Exclusion_Reason":""},
                {"A_ObjectID":10,"B_ObjectID":21,"Final_Constrained_Eligibility":False,"Exclusion_Reason":"different_canonical_nuclear_objectids"},
            ]))


def fixture(tmp_path):
    a={10:[(1,1),(2,1),(3,1)],11:[(5,1),(6,1),(7,1)]}
    b={20:[(1,1),(2,1)],21:[(1,1),(2,1),(3,1)]}
    nuclei={101:[(1,3),(2,3)],102:[(5,3),(6,3)]}
    source(tmp_path/"07_nuclear_colocalization/LHX6_vs_Hoechst/analysis.h5","LHX6",a,[(1,10,101),(2,11,102)],nuclei)
    source(tmp_path/"07_nuclear_colocalization/PV_vs_Hoechst/analysis.h5","PV",b,[(1,20,101),(2,21,102)],nuclei)
    pair(tmp_path/"08_marker_colocalization/LHX6_vs_PV/analysis.h5",8,a,b,{10:1,11:2},{20:1,21:2},[(10,21)])
    pair(tmp_path/"09_marker_colocalization_nucleus_constrained/LHX6_vs_PV/analysis.h5",9,a,b,{10:1,11:2},{20:1,21:2},[(10,20)])
    for marker, roi in (("LHX6",10),("PV",20)):
        folder=tmp_path/"06_final_shape_filtered"/marker
        folder.mkdir(parents=True)
        pd.DataFrame([{"ObjectID":roi,"Area_um2":3.5,"Circularity":.8,"Solidity":.9}]).to_csv(folder/f"{marker}_kept_ROI_metrics.csv",index=False)


def test_requires_marker_when_source_cell_id_ambiguous(tmp_path):
    fixture(tmp_path)
    r=subprocess.run([sys.executable,str(SCRIPT),"--analysis",str(tmp_path),"--cell-id","1"],capture_output=True,text=True)
    assert r.returncode != 0
    assert "specify --marker" in r.stderr


def test_cell_id_query_both_matching_variants(tmp_path):
    fixture(tmp_path)
    r=subprocess.run([sys.executable,str(SCRIPT),"--analysis",str(tmp_path),"--cell-id","1","--marker","LHX6"],capture_output=True,text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    folder=tmp_path/"cell_id_queries/nucleus_101"
    for fname in ("cell_identity_summary.json", "cell_id_report.pdf", "marker_fluorescence.png",
                  "stage08_stage09_overlap_candidates.csv", "stage08_stage09_variant_summary.csv", "marker_measurements.csv"):
        assert (folder/fname).is_file() and (folder/fname).stat().st_size > 0
    rows=pd.read_csv(folder/"stage08_stage09_overlap_candidates.csv")
    assert len(rows)==4
    assert len(rows[(rows.Stage==8)&(rows.Final_Stable_Match)&(rows.Same_Hoechst_Nucleus==False)])==1
    assert len(rows[(rows.Stage==9)&(rows.Final_Stable_Match)&(rows.Same_Hoechst_Nucleus==True)])==1
    m=pd.read_csv(folder/"marker_measurements.csv")
    assert set(m.Marker)=={"LHX6","PV"}
    assert set(m.Stage07_CellID.astype(int))=={1}
    assert int(m.iloc[0].Canonical_Hoechst_Nuclear_ObjectID)==101
    summary=json.loads((folder/"cell_identity_summary.json").read_text())
    assert summary["stage08_final_pairs_involving_query_nucleus"] == 1
    assert summary["stage09_final_pairs_involving_query_nucleus"] == 1


def test_same_numeric_cellid_in_two_markers_can_mean_different_nuclei(tmp_path):
    fixture(tmp_path)
    # Deliberately reassign PV's Stage-07 LOCAL CellID numbering without
    # changing the PV original ObjectID -> nuclear ObjectID associations.
    b={20:[(1,1),(2,1)],21:[(1,1),(2,1),(3,1)]}
    nuclei={101:[(1,3),(2,3)],102:[(5,3),(6,3)]}
    source(tmp_path/"07_nuclear_colocalization/PV_vs_Hoechst/analysis.h5", "PV", b,
           [(2,20,101),(1,21,102)], nuclei)
    a={10:[(1,1),(2,1),(3,1)],11:[(5,1),(6,1),(7,1)]}
    for marker, roi in (("LHX6",11),("PV",21)):
        file = tmp_path / "06_final_shape_filtered" / marker / f"{marker}_kept_ROI_metrics.csv"
        metrics = pd.read_csv(file)
        metrics = pd.concat([metrics, pd.DataFrame([{"ObjectID":roi,"Area_um2":3.5,"Circularity":.8,"Solidity":.9}])], ignore_index=True)
        metrics.to_csv(file,index=False)
    pair(tmp_path/"08_marker_colocalization/LHX6_vs_PV/analysis.h5",8,a,b,{10:1,11:2},{20:2,21:1},[(10,21)])
    pair(tmp_path/"09_marker_colocalization_nucleus_constrained/LHX6_vs_PV/analysis.h5",9,a,b,{10:1,11:2},{20:2,21:1},[(10,20)])
    r=subprocess.run([sys.executable,str(SCRIPT),"--analysis",str(tmp_path),"--cell-id","1","--marker","PV"],capture_output=True,text=True)
    assert r.returncode==0, r.stdout+r.stderr
    result=json.loads((tmp_path/"cell_id_queries/nucleus_102/cell_identity_summary.json").read_text())
    assert result["canonical_Hoechst_nuclear_ObjectID"]==102
    assert result["number_assigned_markers"]==2  # LHX6 C2 + PV C1 are both on nucleus 102.
    r=subprocess.run([sys.executable,str(SCRIPT),"--analysis",str(tmp_path),"--cell-id","1","--marker","LHX6"],capture_output=True,text=True)
    assert r.returncode==0, r.stdout+r.stderr
    result=json.loads((tmp_path/"cell_id_queries/nucleus_101/cell_identity_summary.json").read_text())
    assert result["canonical_Hoechst_nuclear_ObjectID"]==101


def test_stage09_missing_is_explicitly_reported(tmp_path):
    fixture(tmp_path)
    (tmp_path/"09_marker_colocalization_nucleus_constrained/LHX6_vs_PV/analysis.h5").unlink()
    r=subprocess.run([sys.executable,str(SCRIPT),"--analysis",str(tmp_path),"--cell-id","1","--marker","LHX6"],capture_output=True,text=True)
    assert r.returncode == 0, r.stdout+r.stderr
    result=json.loads((tmp_path/"cell_id_queries/nucleus_101/cell_identity_summary.json").read_text())
    assert result["stage08_available"]
    assert result["stage09_available"] is False
    s=pd.read_csv(tmp_path/"cell_id_queries/nucleus_101/stage08_stage09_variant_summary.csv")
    assert s.loc[s.Stage==9,"Status"].iloc[0]=="stage_not_available"


def test_pair_id_queries_are_variant_scoped_and_export_measurements(tmp_path):
    fixture(tmp_path)
    for pair_id, stage, other_present in (
        ("S08_LHX6_PV_P000001", 8, False),
        ("S09_LHX6_PV_P000001", 9, False),
    ):
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--analysis", str(tmp_path), "--pair-id", pair_id],
            capture_output=True, text=True,
        )
        assert r.returncode == 0, r.stdout + r.stderr
        folder = tmp_path / "pair_id_queries" / pair_id
        for name in (
            "pair_query_summary.json", "pair_query_measurements.csv",
            "pair_query_stage08_stage09_cross_reference.csv", "pair_query_report.pdf",
        ):
            assert (folder / name).is_file() and (folder / name).stat().st_size > 0
        summary = json.loads((folder / "pair_query_summary.json").read_text())
        assert summary["queried_stage"] == stage
        assert summary["present_in_stage08"] is (stage == 8 or other_present)
        assert summary["present_in_stage09"] is (stage == 9 or other_present)
        measurements = pd.read_csv(folder / "pair_query_measurements.csv")
        assert measurements.iloc[0].PairID == pair_id
        assert measurements.iloc[0].Overlap_Area_Pixels > 0
        assert measurements.iloc[0].A_Integrated_Fluorescence_Native > 0
        cross = pd.read_csv(folder / "pair_query_stage08_stage09_cross_reference.csv")
        assert {"A_Stage07_Source_CellID", "B_Stage07_Source_CellID",
                "A_Canonical_Nuclear_ObjectID", "B_Canonical_Nuclear_ObjectID"}.issubset(cross.columns)


def test_pair_id_query_cross_references_common_pair(tmp_path):
    fixture(tmp_path)
    a={10:[(1,1),(2,1),(3,1)],11:[(5,1),(6,1),(7,1)]}
    b={20:[(1,1),(2,1)],21:[(1,1),(2,1),(3,1)]}
    pair(tmp_path/"08_marker_colocalization/LHX6_vs_PV/analysis.h5",8,a,b,{10:1,11:2},{20:1,21:2},[(10,20)])
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--analysis", str(tmp_path),
         "--pair-id", "S08_LHX6_PV_P000001"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    folder = tmp_path / "pair_id_queries" / "S08_LHX6_PV_P000001"
    summary = json.loads((folder / "pair_query_summary.json").read_text())
    assert summary["present_in_stage08"] and summary["present_in_stage09"]
    assert summary["stage09_pair_id"] == "S09_LHX6_PV_P000001"
    cross = pd.read_csv(folder / "pair_query_stage08_stage09_cross_reference.csv")
    assert cross.Present.all()
    assert set(cross.Tracking_Category) == {"common_to_both"}


def test_query_rejects_mixed_cell_nuclear_and_pair_selectors(tmp_path):
    fixture(tmp_path)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--analysis", str(tmp_path),
         "--cell-id", "1", "--marker", "LHX6",
         "--pair-id", "S08_LHX6_PV_P000001"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "not allowed with argument" in result.stderr
