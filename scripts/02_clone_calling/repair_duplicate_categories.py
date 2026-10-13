"""Repairs h5ad files that fail to load with 'Categorical categories must be unique'.

Usage:
    python repair_duplicate_categories.py broken.h5ad repaired.h5ad
"""

import sys
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from anndata._io.specs import read_elem


def decode(val):
    return val.decode() if isinstance(val, bytes) else str(val)


def read_obs_with_repair(obs_group):
    """Read every obs column individually; for any categorical column that
    fails with the duplicate-categories error, deduplicate its categories
    and remap codes to point at the first occurrence of each duplicate's
    decoded value, rather than letting the whole file fail to load."""
    index_key = obs_group.attrs.get("_index", "_index")
    obs_index = read_elem(obs_group[index_key])

    cols = {}
    for key in obs_group.keys():
        if key == index_key:
            continue
        try:
            cols[key] = read_elem(obs_group[key])
        except ValueError as e:
            if "Categorical categories must be unique" not in str(e):
                raise
            print(f"  [repair] obs['{key}'] has duplicate categories -- fixing")
            raw_categories = read_elem(obs_group[key]["categories"])
            codes = read_elem(obs_group[key]["codes"])

            decoded = [decode(c) for c in raw_categories]
            seen, new_categories, remap = {}, [], {}
            for old_code, val in enumerate(decoded):
                if val not in seen:
                    seen[val] = len(new_categories)
                    new_categories.append(val)
                remap[old_code] = seen[val]

            new_codes = np.array([remap.get(c, -1) for c in codes])
            cols[key] = pd.Categorical.from_codes(new_codes, categories=new_categories)

            n_dup = len(raw_categories) - len(new_categories)
            print(f"    {len(raw_categories)} categories -> {len(new_categories)} "
                  f"({n_dup} duplicate(s) removed)")

    return pd.DataFrame(cols, index=obs_index)


def load_with_repair(path):
    with h5py.File(path, "r") as f:
        X = read_elem(f["X"])
        var = read_elem(f["var"])
        obs = read_obs_with_repair(f["obs"])
        obsm = {k: read_elem(f["obsm"][k]) for k in f["obsm"]} if "obsm" in f else {}
        uns = {k: read_elem(f["uns"][k]) for k in f["uns"]} if "uns" in f else {}
        layers = {k: read_elem(f["layers"][k]) for k in f["layers"]} if "layers" in f else {}
    return ad.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns, layers=layers)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    in_path = sys.argv[1]
    diagnose_only = "--diagnose-only" in sys.argv

    print(f"Reading {in_path} (working around the crash, not triggering it)...")
    adata = load_with_repair(in_path)
    print(f"Loaded successfully: {adata.n_obs} cells x {adata.n_vars} genes")

    if diagnose_only:
        sys.exit(0)

    if len(sys.argv) < 3:
        print("Provide an output path as the second argument to write the repaired file.")
        sys.exit(1)

    out_path = sys.argv[2]
    adata.write_h5ad(out_path)
    print(f"-> {out_path}")
