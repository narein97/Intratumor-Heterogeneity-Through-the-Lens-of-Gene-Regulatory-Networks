#!/usr/bin/env Rscript
# Runs inferCNV on one EOC patient with the parameters of Zhang et al. 2022 (inferCNV 1.4.0).
# Usage: Rscript run_infercnv_zhang2022.R --counts split_patients/EOC372_counts.tsv --annotations split_patients/EOC372_annotations.tsv --gene_order_file gencode_v25_gene_order.txt --stromal_labels Stromal --malignant_labels EOC --out_dir infercnv_out/EOC372 --seed 1

suppressPackageStartupMessages({
    library(infercnv)
})

# =============================================================================
# 0. Argument parsing (base R only)
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL, required = FALSE) {
    idx <- which(args == flag)
    if (length(idx) == 0) {
        if (required) {
            cat(sprintf("[ERROR] Missing required argument: %s\n", flag), file = stderr())
            quit(status = 1)
        }
        return(default)
    }
    args[idx + 1]
}

counts_file       <- get_arg("--counts", required = TRUE)
annotations_file  <- get_arg("--annotations", required = TRUE)
gene_order_file   <- get_arg("--gene_order_file", required = TRUE)
stromal_labels    <- get_arg("--stromal_labels", required = TRUE)
malignant_labels  <- get_arg("--malignant_labels", required = TRUE)
out_dir           <- get_arg("--out_dir", required = TRUE)
seed              <- as.integer(get_arg("--seed", default = "1"))
n_reference_cells <- as.integer(get_arg("--n_reference_cells", default = "150"))  # "up to 150" per paper
min_subclone_size <- as.integer(get_arg("--min_subclone_size", default = "5"))    # paper: "less than five cells"

stromal_labels <- strsplit(stromal_labels, ",")[[1]]
malignant_labels <- strsplit(malignant_labels, ",")[[1]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

log("═══════════════════════════════════════════════════")
log("inferCNV per Zhang et al. 2022 Methods")
log("═══════════════════════════════════════════════════")
log("counts             : ", counts_file)
log("annotations        : ", annotations_file)
log("gene_order_file    : ", gene_order_file)
log("stromal labels     : ", paste(stromal_labels, collapse = ", "))
log("malignant labels   : ", paste(malignant_labels, collapse = ", "),
    " (anything NOT in stromal_labels or malignant_labels -- e.g. Immune -- is EXCLUDED entirely)")
log("out_dir            : ", out_dir)
log("random seed        : ", seed, " (NOT specified in the paper -- see script header)")
log("reference cells     up to ", n_reference_cells, " stromal cells (paper: 'up to 150')")
log("min subclone size  : ", min_subclone_size, " cells (paper: 'less than five cells' filtered out)")
log("═══════════════════════════════════════════════════")

log("Loading annotations...")
annot <- read.delim(annotations_file, header = FALSE, sep = "\t",
                    col.names = c("cell", "cell_type"), stringsAsFactors = FALSE)

stromal_cells <- annot$cell[annot$cell_type %in% stromal_labels]
if (length(stromal_cells) == 0) {
    log("[ERROR] No cells matched stromal_labels (", paste(stromal_labels, collapse=", "),
        ") in ", annotations_file, " -- check --stromal_labels matches your own annotation strings")
    quit(status = 1)
}

n_sample <- min(n_reference_cells, length(stromal_cells))
reference_cells <- sample(stromal_cells, n_sample)
log("Stromal cells available: ", length(stromal_cells),
    "  |  sampled as reference: ", n_sample)

annot_filtered <- annot[annot$cell %in% reference_cells | annot$cell_type %in% malignant_labels, ]

n_excluded <- sum(!(annot$cell_type %in% stromal_labels) & !(annot$cell_type %in% malignant_labels))
if (n_excluded > 0) {
    excluded_types <- unique(annot$cell_type[!(annot$cell_type %in% stromal_labels) &
                                             !(annot$cell_type %in% malignant_labels)])
    log("Excluding ", n_excluded, " cell(s) that are neither stromal nor malignant ",
        "(cell_type: ", paste(excluded_types, collapse = ", "), ") -- not fed to inferCNV at all")
}

annot_out <- file.path(out_dir, "infercnv_annotations_filtered.tsv")
write.table(annot_filtered, annot_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
log("Filtered annotation file (reference-subsampled) -> ", annot_out)

reference_group_names <- unique(annot_filtered$cell_type[annot_filtered$cell %in% reference_cells])
log("Reference group name(s) passed to inferCNV: ", paste(reference_group_names, collapse = ", "))

log("Creating inferCNV object...")
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_file,
    annotations_file  = annot_out,
    delim             = "\t",
    gene_order_file   = gene_order_file,
    ref_group_names   = reference_group_names
)

log("Running inferCNV (this can take a long time depending on cell count)...")
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff                          = 0.1,
    out_dir                         = out_dir,
    denoise                         = TRUE,
    HMM                             = TRUE,
    hclust_method                   = "ward.D2",
    analysis_mode                   = "subclusters",
    tumor_subcluster_partition_method = "random_trees",
    tumor_subcluster_pval           = 0.05,
    num_threads                     = 10
)
log("inferCNV run complete. Raw output -> ", out_dir)

log("Extracting subclone assignments from infercnv.observation_groupings.txt...")

groupings_path <- file.path(out_dir, "infercnv.observation_groupings.txt")
if (!file.exists(groupings_path)) {
    log("[ERROR] ", groupings_path, " not found -- cannot extract subclones")
    quit(status = 1)
}

groupings <- read.table(groupings_path, header = TRUE, row.names = 1,
                        sep = " ", quote = "\"", stringsAsFactors = FALSE)
if (!"Dendrogram.Group" %in% colnames(groupings)) {
    log("[ERROR] Expected column 'Dendrogram.Group' not found in ", groupings_path,
        " -- found instead: ", paste(colnames(groupings), collapse = ", "),
        " (inferCNV's output column naming may differ across versions -- see script header)")
    quit(status = 1)
}

subclone_df <- data.frame(
    cell = rownames(groupings),
    observation_group = "observation",   # all cells in this file are, by
                                          # definition, observation-group
                                          # cells; reference cells never
                                          # appear here
    subclone = paste0("subcluster_", groupings$Dendrogram.Group),
    stringsAsFactors = FALSE
)

n_before <- length(unique(subclone_df$subclone))
log("Subclones found (before size filter): ", n_before, " -- ",
    paste(names(table(subclone_df$subclone)), table(subclone_df$subclone), sep="=", collapse=", "))

# Filter subclones with fewer than min_subclone_size cells
subclone_sizes <- table(subclone_df$subclone)
small_subclones <- names(subclone_sizes[subclone_sizes < min_subclone_size])
if (length(small_subclones) > 0) {
    log("Removing ", length(small_subclones), " subclone(s) with < ",
        min_subclone_size, " cells: ", paste(small_subclones, collapse = ", "))
    subclone_df <- subclone_df[!(subclone_df$subclone %in% small_subclones), ]
}
n_after <- length(unique(subclone_df$subclone))

subclone_out <- file.path(out_dir, "subclone_assignments_filtered.tsv")
write.table(subclone_df, subclone_out, sep = "\t", row.names = FALSE, quote = FALSE)

log("Subclones before size filter: ", n_before, "  |  after: ", n_after)
log("Subclone assignment table -> ", subclone_out)

log("═══════════════════════════════════════════════════")
log("DONE.")
log("Next step (separate tool, not run by this script):")
log("  UPhyloplot2 on this same --out_dir to generate phylogenetic trees")
log("  https://github.com/harbourlab/UPhyloplot2")
log("═══════════════════════════════════════════════════")
