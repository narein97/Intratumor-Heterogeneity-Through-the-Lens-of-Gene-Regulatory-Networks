#!/usr/bin/env python3
"""Writes copies of clone-labelled patient h5ads with all malignant clones pooled into one tumour population 'T' (normal cells stay 'N').

Usage:
    python make_pooled_tumor_h5ad.py --in-dir stage4_out/Breast/Stage4 --out-dir stage4_out/Breast/Stage4_pooled
"""
import argparse
import glob
import os

import anndata as ad
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", required=True, help="directory with <sample>_clones.h5ad files")
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(a.in_dir, "*_clones.h5ad"))):
        ad_ = ad.read_h5ad(f)
        for k in [k for k in ad_.uns if k.startswith("celloracle") or k in ("clone_categories", "bootstrap", "shuffle")]:
            del ad_.uns[k]
        ad_.obs["clone"] = pd.Categorical(ad_.obs["clone"].astype(str).where(ad_.obs["clone"].astype(str) == "N", "T"))
        out = os.path.join(a.out_dir, os.path.basename(f))
        ad_.write_h5ad(out)
        print(f"[ok] {f} -> {out}  ({(ad_.obs['clone'] == 'T').sum()} T, {(ad_.obs['clone'] == 'N').sum()} N)")


if __name__ == "__main__":
    main()
