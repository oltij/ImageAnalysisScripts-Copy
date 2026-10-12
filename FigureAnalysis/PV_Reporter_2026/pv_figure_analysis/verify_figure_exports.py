#!/usr/bin/env python3
"""Verify twelve publication figure sets and editable SVG text in V1 and V2."""
import argparse
from pathlib import Path
from xml.etree import ElementTree as ET

BASE_STEMS = [
    'ABC_PV_area_specificity_coverage_publication',
    'BC_only_specificity_coverage_publication',
    'A_PV_area_publication',
    'A_supp_PV_density_publication',
    'B_BiVe3_specificity_coverage_publication',
    'C_BiVe4_specificity_coverage_publication',
]
STEMS = [f'{name}_{version}' for version in ('v1', 'v2') for name in BASE_STEMS]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--figures', type=Path, required=True)
    args = p.parse_args()
    for stem in STEMS:
        for ext in ('svg','pdf','png'):
            f = args.figures / f'{stem}.{ext}'
            if not f.is_file() or f.stat().st_size <= 0:
                raise SystemExit(f'Missing or empty: {f}')
        svg = args.figures / f'{stem}.svg'
        txt = ET.parse(svg).getroot().findall('.//{http://www.w3.org/2000/svg}text')
        if len(txt) < 3:
            raise SystemExit(f'No editable text found: {svg}')
        print(f'OK  {stem}: SVG/PDF/PNG + {len(txt)} live text objects')
    print('All twelve figure sets passed (36 files).')

if __name__ == '__main__': main()
