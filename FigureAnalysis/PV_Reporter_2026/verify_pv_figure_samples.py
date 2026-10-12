#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

ROOT = Path(
    "/gpfs/accounts/parent_root/parent0/oltij/"
    "Neuron_Rebuttal_IMS_Batch50_MaskFix_20261009/"
    "_PV_selected_figures"
)

SOURCE = ROOT / "per_acquisition_metrics.csv"
KW_SOURCE = (
    ROOT / "panelA_kw" /
    "per_acquisition_metrics_with_nuclear_associated_PV_area.csv"
)

BLOCK_TO_EXPERIMENT = {
    "165": "24", "166": "24",
    "167": "24", "168": "24",
    "169": "RA02", "170": "RA02",
    "171": "RA02", "172": "RA02",
    "173": "RA02", "174": "RA02",
    "175": "17", "176": "17",
    "177": "17", "178": "17",
}

LINES = ["H1", "H1 CL32", "H1 CL23"]
EXPERIMENTS = ["17", "24", "RA02"]
ENHANCERS = ["BiVe3", "BiVe4"]


def load_csv(path):
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        block = row["experiment"].split("-")[0]
        row["_experiment_id"] = BLOCK_TO_EXPERIMENT[block]

    return rows


rows = load_csv(SOURCE)

print("\nSOURCE FILE:")
print(SOURCE)
print("\nTOTAL SELECTED ACQUISITIONS:", len(rows))

ids = [r["acquisition_label"] for r in rows]
assert len(ids) == len(set(ids)), "Duplicate acquisition IDs"

print("\n" + "=" * 75)
print("COMPLETE SAMPLE BREAKDOWN")
print("=" * 75)

for enhancer in ENHANCERS:
    print(f"\n{enhancer}")
    print("-" * 75)

    for line in LINES:
        for exp in EXPERIMENTS:
            selected = [
                r for r in rows
                if r["enhancer"] == enhancer
                and r["cell_line"] == line
                and r["_experiment_id"] == exp
            ]

            labels = [r["acquisition_label"] for r in selected]

            print(
                f"{line:10s} | Experiment {exp:4s} | "
                f"n = {len(selected):2d} | "
                f"{', '.join(labels) if labels else 'NONE'}"
            )

    total = sum(r["enhancer"] == enhancer for r in rows)
    print(f"\nTOTAL {enhancer}: {total}")

print("\n" + "=" * 75)
print("SPECIFIC COUNTS FROM YOUR FIGURE LEGENDS")
print("=" * 75)

checks = [
    ("BiVe4", "H1 CL32", "24", 2),
    ("BiVe4", "H1 CL23", "RA02", 3),
    ("BiVe4", "H1 CL32", "RA02", 3),
    ("BiVe3", "H1 CL23", "RA02", 2),
]

for enhancer, line, exp, expected in checks:
    selected = [
        r for r in rows
        if r["enhancer"] == enhancer
        and r["cell_line"] == line
        and r["_experiment_id"] == exp
    ]

    actual = len(selected)
    status = "PASS" if actual == expected else "MISMATCH"

    print(
        f"{status}: {enhancer}, {line}, Exp {exp}: "
        f"expected {expected}, actual {actual}"
    )

    for row in selected:
        print(
            f"    {row['acquisition_label']} "
            f"| {row['source_ims_filename']}"
        )

print("\n" + "=" * 75)
print("CHECK AGAINST KRUSKAL-WALLIS FIGURE INPUT")
print("=" * 75)

if KW_SOURCE.exists():
    kw_rows = load_csv(KW_SOURCE)

    source_ids = set(ids)
    kw_ids = {r["acquisition_label"] for r in kw_rows}

    missing = sorted(source_ids - kw_ids)
    extra = sorted(kw_ids - source_ids)

    print("Source rows:", len(rows))
    print("KW analysis rows:", len(kw_rows))
    print("Missing in KW:", missing)
    print("Extra in KW:", extra)

    source_metadata = {
        r["acquisition_label"]: (
            r["experiment"],
            r["cell_line"],
            r["enhancer"],
        )
        for r in rows
    }

    mismatched = [
        r["acquisition_label"]
        for r in kw_rows
        if source_metadata.get(r["acquisition_label"]) != (
            r["experiment"],
            r["cell_line"],
            r["enhancer"],
        )
    ]

    print("Metadata mismatches:", mismatched)

    assert not missing and not extra and not mismatched
    print("\nPASS: All source acquisitions match the KW analysis input.")

else:
    print("KW-derived input file not found:", KW_SOURCE)

print("\nAUDIT COMPLETE")
