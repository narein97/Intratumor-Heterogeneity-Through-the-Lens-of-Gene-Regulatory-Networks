#!/usr/bin/env python3
"""Adds 100 bootstrap and 100 label-shuffle definitions to each patient h5ad.

Usage:
    python add_bootstrap_shuffle.py --root clones_out --inplace
"""

import os
import glob
import argparse
import logging

import numpy as np
import anndata as ad

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("boot_shuffle")


def clone_categories(labels):
    """Order categories: malignant integer clones ascending, then 'N' last."""
    uniq = list(dict.fromkeys(labels))
    mal = sorted((c for c in uniq if c != "N"),
                 key=lambda c: (int(c) if str(c).isdigit() else 1e9, c))
    return mal + (["N"] if "N" in uniq else [])


def build_bootstrap(labels, cats, n_boot, frac, rng):
    """Selected cell indices per replicate: frac sampled w/o replacement per clone."""
    per_clone_idx = {c: np.where(labels == c)[0] for c in cats}
    # fixed replicate size => a single 2-D array
    sizes = {c: max(1, int(frac * len(idx))) for c, idx in per_clone_idx.items()}
    k = int(sum(sizes.values()))
    out = np.empty((n_boot, k), dtype=np.int32)
    for i in range(n_boot):
        sel = [rng.choice(per_clone_idx[c], size=sizes[c], replace=False)
               for c in cats]
        out[i] = np.concatenate(sel)
    return out, k, sizes


def build_shuffle(labels, cats, n_shuffle, rng):
    """Permuted clone codes per replicate: non-N labels permuted, N fixed."""
    code = {c: i for i, c in enumerate(cats)}
    codes = np.array([code[c] for c in labels], dtype=np.int16)
    nonN = np.where(labels != "N")[0]
    out = np.empty((n_shuffle, len(labels)), dtype=np.int16)
    for i in range(n_shuffle):
        row = codes.copy()
        row[nonN] = rng.permutation(codes[nonN])
        out[i] = row
    return out


def augment(path, out_path, clone_col, n_boot, n_shuffle, frac, seed):
    adata = ad.read_h5ad(path)
    if clone_col not in adata.obs:
        log.warning("[skip] %s -- no obs['%s']", path, clone_col)
        return False
    labels = adata.obs[clone_col].astype(str).to_numpy()
    cats = clone_categories(labels)
    if len(cats) < 2:
        log.warning("[skip] %s -- < 2 clone categories", path)
        return False

    rng = np.random.default_rng(seed)
    boot_idx, k, sizes = build_bootstrap(labels, cats, n_boot, frac, rng)
    shuf = build_shuffle(labels, cats, n_shuffle, rng)

    adata.uns["clone_categories"] = np.array(cats, dtype=object)
    adata.uns["bootstrap"] = {"indices": boot_idx,
                              "frac": float(frac), "n": int(n_boot)}
    adata.uns["shuffle"] = {"labels": shuf, "n": int(n_shuffle)}

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    adata.write_h5ad(out_path)
    log.info("[ok]   %s  cats=%s boot=(%d,%d) shuffle=(%d,%d) -> %s",
             os.path.basename(path), cats, n_boot, k, n_shuffle,
             adata.n_obs, out_path)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", help="single patient .h5ad")
    p.add_argument("--output", help="output .h5ad (default: overwrite --input)")
    p.add_argument("--root", help="walk this dir for *_clones.h5ad instead")
    p.add_argument("--inplace", action="store_true",
                   help="[--root mode] overwrite each file in place")
    p.add_argument("--out-root", help="[--root mode] mirror tree under here")
    p.add_argument("--clone-col", default="clone")
    p.add_argument("--n-boot", type=int, default=100)
    p.add_argument("--n-shuffle", type=int, default=100)
    p.add_argument("--boot-frac", type=float, default=0.9)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    if args.input:
        out = args.output or args.input
        augment(args.input, out, args.clone_col,
                args.n_boot, args.n_shuffle, args.boot_frac, args.seed)
    elif args.root:
        files = sorted(glob.glob(os.path.join(args.root, "**", "*_clones.h5ad"),
                                 recursive=True))
        log.info("found %d h5ad file(s) under %s", len(files), args.root)
        for f in files:
            if args.out_root:
                rel = os.path.relpath(f, args.root)
                out = os.path.join(args.out_root, rel)
            else:  # default in-place
                out = f
            augment(f, out, args.clone_col,
                    args.n_boot, args.n_shuffle, args.boot_frac, args.seed)
    else:
        p.error("provide --input FILE or --root DIR")


if __name__ == "__main__":
    main()
