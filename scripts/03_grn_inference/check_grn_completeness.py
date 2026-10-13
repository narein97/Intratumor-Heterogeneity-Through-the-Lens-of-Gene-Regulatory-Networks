#!/usr/bin/env python3
"""Checks that every patient h5ad has the original, 100 bootstrap and 100 shuffle GRNs.

Usage:
    python check_grn_completeness.py --root clones_out
"""

import argparse
import glob
import anndata as ad

EXPECTED_BOOTSTRAP = 100
EXPECTED_SHUFFLE = 100
EXPECTED_ORIGINAL = 1


def check_one(h5ad_path):
    """Returns (n_bootstrap, n_shuffle, n_original, error_or_None)."""
    try:
        adata = ad.read_h5ad(h5ad_path, backed="r")
    except Exception as e:
        return None, None, None, f"{type(e).__name__}: {e}"

    n_bootstrap = len(adata.uns["celloracle_bootstrap"]) if "celloracle_bootstrap" in adata.uns else 0
    n_shuffle = len(adata.uns["celloracle_shuffle"]) if "celloracle_shuffle" in adata.uns else 0
    n_original = 1 if ("celloracle_original" in adata.uns and adata.uns["celloracle_original"] is not None
                       and len(adata.uns["celloracle_original"]) > 0) else 0

    return n_bootstrap, n_shuffle, n_original, None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--exclude-cancer-type", nargs="*", default=["Ovarian"])
    args = p.parse_args()

    files = sorted(glob.glob(f"{args.root}/**/*.h5ad", recursive=True))
    files = [f for f in files if "_with_support" not in f]
    if args.exclude_cancer_type:
        files = [f for f in files if not any(f"/{ct}/" in f for ct in args.exclude_cancer_type)]

    print(f"Checking {len(files)} h5ad(s) under {args.root}...\n")

    n_complete = 0
    n_incomplete = 0

    for f in files:
        n_boot, n_shuf, n_orig, error = check_one(f)

        if error is not None:
            n_incomplete += 1
            print(f"{f}")
            print(f"    COULD NOT OPEN: {error}")
            continue

        is_complete = (n_boot == EXPECTED_BOOTSTRAP and
                      n_shuf == EXPECTED_SHUFFLE and
                      n_orig == EXPECTED_ORIGINAL)

        if is_complete:
            n_complete += 1
        else:
            n_incomplete += 1
            print(f"{f}")
            print(f"    bootstrap: {n_boot}/{EXPECTED_BOOTSTRAP}   "
                  f"shuffle: {n_shuf}/{EXPECTED_SHUFFLE}   "
                  f"original: {n_orig}/{EXPECTED_ORIGINAL}")

    print(f"\n{'='*60}")
    print(f"  Complete   : {n_complete}")
    print(f"  Incomplete : {n_incomplete}")
    print(f"  Total      : {len(files)}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
