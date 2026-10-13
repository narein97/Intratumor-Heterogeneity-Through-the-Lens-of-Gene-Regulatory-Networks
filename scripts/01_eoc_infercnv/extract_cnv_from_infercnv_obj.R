#!/usr/bin/env Rscript
# Writes infercnv.observations.txt from the final inferCNV object (for inferCNV versions that do not write it).
# Usage: Rscript extract_cnv_from_infercnv_obj.R infercnv_out/EOC349

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("Usage: Rscript extract_cnv_from_infercnv_obj.R <infercnv_output_dir>\n")
    quit(status = 1)
}
infercnv_dir <- args[1]
log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# Find the final-stage object by highest step number, rather than hardcoding
# the exact filename -- the parameter-encoding suffix (NF_NA.SD_1.5...) could
# vary slightly and shouldn't need a code change to match.
candidates <- list.files(infercnv_dir, pattern = "^[0-9]+.*denoise.*\\.infercnv_obj$",
                         full.names = TRUE)
if (length(candidates) == 0) {
    log("[ERROR] No '*denoise*.infercnv_obj' file found in ", infercnv_dir)
    log("        Files present: ", paste(list.files(infercnv_dir, pattern = "infercnv_obj$"), collapse = ", "))
    quit(status = 1)
}
step_numbers <- as.integer(sub("^([0-9]+).*", "\\1", basename(candidates)))
obj_path <- candidates[which.max(step_numbers)]
log("Using: ", obj_path)

infercnv_obj <- readRDS(obj_path)

if (!"expr.data" %in% slotNames(infercnv_obj)) {
    log("[ERROR] Object has no 'expr.data' slot -- the @expr.data guess was wrong.")
    log("        Actual slot names: ", paste(slotNames(infercnv_obj), collapse = ", "))
    log("        Please share this output so the correct slot can be identified.")
    quit(status = 1)
}

expr_matrix <- infercnv_obj@expr.data
if (is.null(expr_matrix) || nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    log("[ERROR] @expr.data exists but is empty (", nrow(expr_matrix), " x ", ncol(expr_matrix), ")")
    quit(status = 1)
}

log("Extracted matrix: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")

out_path <- file.path(infercnv_dir, "infercnv.observations.txt")
# Same format the old plain-text output used, and what load_cnv_observations()
# in build_h5ad_per_patient.py already parses: quoted names, space-delimited,
# genes as rows, cells as columns.
write.table(expr_matrix, out_path, sep = " ", quote = TRUE, row.names = TRUE, col.names = TRUE)
log("-> ", out_path)
