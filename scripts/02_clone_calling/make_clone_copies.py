"""Writes <patient>_clones.h5ad with obs['clone'] set to the clone-sweep clones (obs['clone_leiden']) for the EOC patients.

Usage:
    cd clones_out && python make_clone_copies.py
"""
import glob
import anndata as ad

for path in glob.glob("*.h5ad"):
    if path.endswith("_clones.h5ad"):
        continue  # don't re-process outputs from a previous run

    patient = path[:-len(".h5ad")]
    a = ad.read_h5ad(path)

    if "clone_leiden" not in a.obs:
        print(f"[skip] {path} -- no obs['clone_leiden']")
        continue

    a.obs["clone"] = a.obs["clone_leiden"]
    out_path = f"{patient}_clones.h5ad"
    a.write_h5ad(out_path)
    print(f"[ok] {path} -> {out_path}")
