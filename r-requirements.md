# R requirements (EOC inferCNV step only)

- R >= 4.0
- `infercnv` **1.4.0** (the version used by Zhang et al. 2022). Newer inferCNV releases no longer write
  `infercnv.observations.txt`; `scripts/01_eoc_infercnv/extract_cnv_from_infercnv_obj.R` recovers the matrix
  from the serialized `.infercnv_obj` in that case (this extraction path was written without access to a
  real newer object, so check its output before relying on it).
- Argument parsing is base R only; no other R packages are needed.
