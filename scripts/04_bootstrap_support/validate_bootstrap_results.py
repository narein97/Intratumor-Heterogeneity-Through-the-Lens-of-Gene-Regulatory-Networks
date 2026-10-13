#!/usr/bin/env python3
"""Quality-control pass over the patient h5ads and their *_with_support companions.

Usage:
    python validate_bootstrap_results.py --root clones_out
"""

import argparse
import glob
import os
import numpy as np
import pandas as pd
import anndata as ad


REQUIRED_SOURCE_KEYS = ["celloracle_original", "celloracle_bootstrap", "celloracle_shuffle",
                        "celloracle_columns", "celloracle_genes", "clone_categories",
                        "cnp_genes"]
REQUIRED_OBS_COLUMNS = ["clone", "cell_type", "patient", "sample"]


def support_output_path(source_path, suffix="_with_support"):
    base, ext = os.path.splitext(source_path)
    return f"{base}{suffix}{ext}"


def validate_source(h5ad_path):
    issues = []
    try:
        adata = ad.read_h5ad(h5ad_path, backed="r")
    except Exception as e:
        issues.append(("ERROR", f"could not open source file: {type(e).__name__}: {e}"))
        return issues, None

    missing_uns = [k for k in REQUIRED_SOURCE_KEYS if k not in adata.uns]
    if missing_uns:
        issues.append(("ERROR", f"source missing .uns key(s): {missing_uns}"))

    missing_obs = [c for c in REQUIRED_OBS_COLUMNS if c not in adata.obs.columns]
    if missing_obs:
        issues.append(("ERROR", f"source missing obs column(s): {missing_obs}"))

    if "cnp" not in adata.obsm:
        issues.append(("ERROR", "source missing obsm['cnp']"))
    elif "cnp_genes" in adata.uns:
        cnp = adata.obsm["cnp"]
        genes = adata.uns["cnp_genes"]
        if cnp.shape[1] != len(genes):
            issues.append(("ERROR", f"obsm['cnp'] has {cnp.shape[1]} columns but "
                                   f"uns['cnp_genes'] has {len(genes)} entries -- mismatch"))

    if "clone" in adata.obs.columns:
        n_nan_clone = adata.obs["clone"].isna().sum()
        if n_nan_clone > 0:
            issues.append(("ERROR", f"{n_nan_clone} cell(s) have NaN in obs['clone']"))

    return issues, adata


def validate_support_output(source_path, p_thresh=0.05, expected_n_boot=100, expected_n_shuffle=100):
    """Checks the SEPARATE companion file, not the source."""
    issues = []
    out_path = support_output_path(source_path)

    if not os.path.exists(out_path):
        issues.append(("INFO", f"no support output yet at {os.path.basename(out_path)} "
                              f"(compute_bootstrap_support.py hasn't been run for this patient)"))
        return issues

    try:
        adata = ad.read_h5ad(out_path, backed="r")
    except Exception as e:
        issues.append(("ERROR", f"could not open support output file: {type(e).__name__}: {e}"))
        return issues

    n_boot = len(adata.uns.get("celloracle_bootstrap", {}))
    n_shuffle = len(adata.uns.get("celloracle_shuffle", {}))
    if n_boot != expected_n_boot:
        issues.append(("WARN", f"only {n_boot}/{expected_n_boot} bootstrap replicates -- "
                               f"support values here are out of {n_boot}, not {expected_n_boot}"))
    if n_shuffle != expected_n_shuffle:
        issues.append(("WARN", f"only {n_shuffle}/{expected_n_shuffle} shuffle replicates"))

    if "celloracle_original_filtered" not in adata.uns:
        issues.append(("ERROR", "support output file exists but has no celloracle_original_filtered"))
        return issues

    filtered = adata.uns["celloracle_original_filtered"]
    if not isinstance(filtered, pd.DataFrame):
        filtered = pd.DataFrame(filtered)

    for col in ["clone", "source", "target", "p", "bootstrap_support"]:
        if col not in filtered.columns:
            issues.append(("ERROR", f"celloracle_original_filtered missing column '{col}'"))
    if any("missing column" in msg for _, msg in issues):
        return issues

    bad_p = filtered[filtered["p"] >= p_thresh]
    if len(bad_p) > 0:
        issues.append(("ERROR", f"{len(bad_p)} row(s) have p >= {p_thresh} -- filter failed"))

    support = filtered["bootstrap_support"]
    if (support < 0).any():
        issues.append(("ERROR", f"{(support < 0).sum()} row(s) have negative bootstrap_support"))
    if (support > n_boot).any():
        issues.append(("ERROR", f"{(support > n_boot).sum()} row(s) exceed the actual "
                               f"bootstrap replicate count ({n_boot})"))

    real_clones = set(np.asarray(adata.uns.get("clone_categories", [])).astype(str))
    bad_clones = set(filtered["clone"].astype(str)) - real_clones
    if bad_clones:
        issues.append(("ERROR", f"clone value(s) not in clone_categories: {bad_clones}"))

    if len(filtered) == 0:
        issues.append(("WARN", "celloracle_original_filtered is empty"))

    return issues


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--h5ad", default=None)
    p.add_argument("--p-thresh", type=float, default=0.05)
    p.add_argument("--expected-n-boot", type=int, default=100)
    p.add_argument("--expected-n-shuffle", type=int, default=100)
    p.add_argument("--exclude-cancer-type", nargs="*", default=["Ovarian"])
    p.add_argument("--skip-support-check", action="store_true",
                   help="only validate source files, don't look for companion output")
    args = p.parse_args()

    if args.h5ad:
        files = [args.h5ad]
    else:
        files = sorted(glob.glob(f"{args.root}/**/*_clones.h5ad", recursive=True))
        # never match a support-output file as if it were a source file
        files = [f for f in files if "_with_support" not in f]
        if args.exclude_cancer_type:
            files = [f for f in files if not any(f"/{ct}/" in f for ct in args.exclude_cancer_type)]

    print(f"Checking {len(files)} source patient h5ad(s)...\n")

    n_source_clean = n_source_error = 0
    n_support_missing = n_support_clean = n_support_warn = n_support_error = 0

    for f in files:
        src_issues, _ = validate_source(f)
        src_errors = [i for i in src_issues if i[0] == "ERROR"]
        if src_errors:
            n_source_error += 1
            print(f"[SOURCE ERROR] {f}")
            for _, msg in src_errors:
                print(f"    {msg}")
        else:
            n_source_clean += 1

        if not args.skip_support_check:
            sup_issues = validate_support_output(f, args.p_thresh, args.expected_n_boot, args.expected_n_shuffle)
            sup_errors = [i for i in sup_issues if i[0] == "ERROR"]
            sup_warns = [i for i in sup_issues if i[0] == "WARN"]
            sup_info = [i for i in sup_issues if i[0] == "INFO"]

            if sup_errors:
                n_support_error += 1
                print(f"[SUPPORT ERROR] {support_output_path(f)}")
                for _, msg in sup_errors:
                    print(f"    {msg}")
            elif sup_warns:
                n_support_warn += 1
                print(f"[SUPPORT WARN]  {support_output_path(f)}")
                for _, msg in sup_warns:
                    print(f"    {msg}")
            elif sup_info:
                n_support_missing += 1
            else:
                n_support_clean += 1

    print(f"\n{'='*60}")
    print("SOURCE FILES")
    print(f"  Clean : {n_source_clean}   Errors : {n_source_error}   Total : {len(files)}")
    if not args.skip_support_check:
        print("\nSUPPORT OUTPUT FILES")
        print(f"  Clean : {n_support_clean}   Warnings : {n_support_warn}   "
              f"Errors : {n_support_error}   Not yet generated : {n_support_missing}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
