#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript isoquant_read_info_to_qs.R read_assignments.tsv[.gz] output.qs [sample]")
}

suppressPackageStartupMessages({
    library(data.table)
    library(qs)
})

input <- args[1]
output <- args[2]
sample <- if (length(args) >= 3) args[3] else NA_character_

assignments <- if (grepl("\\.gz$", input, ignore.case = TRUE)) {
    fread(cmd = paste("gzip -dc --", shQuote(input)), skip = "read_id",
          sep = "\t", header = TRUE, data.table = TRUE,
          na.strings = c("", ".", "NA"), showProgress = TRUE)
} else {
    fread(file = input, skip = "read_id", sep = "\t", header = TRUE,
          data.table = TRUE, na.strings = c("", ".", "NA"),
          showProgress = TRUE)
}

if ("isoform_assignment_type" %in% names(assignments) &&
    !"assignment_type" %in% names(assignments)) {
    setnames(assignments, "isoform_assignment_type", "assignment_type")
}
required <- c("read_id", "gene_id", "isoform_id", "assignment_type")
missing_columns <- setdiff(required, names(assignments))
if (length(missing_columns) > 0) {
    stop("IsoQuant read_info is missing columns: ",
         paste(missing_columns, collapse = ", "))
}

setnames(assignments,
         c("read_id", "isoform_id", "assignment_type"),
         c("readid", "transcript_id", "type"))

# IsoQuant intentionally emits multiple rows for reads with multiple candidate
# assignments. Keep those rows so the Fullscope all-assignment matrix remains
# lossless; the standard merge script creates a second unique-transcript matrix.
assignments <- assignments[
    !is.na(readid) & nzchar(readid) &
    !is.na(gene_id) & nzchar(gene_id) &
    !is.na(transcript_id) & nzchar(transcript_id)
]
assignments[, sample := sample]
assignments[, assignment_source := "isoquant"]
setcolorder(assignments,
            c("readid", "gene_id", "transcript_id", "type", "sample",
              "assignment_source",
              setdiff(names(assignments),
                      c("readid", "gene_id", "transcript_id", "type", "sample",
                        "assignment_source"))))

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
qsave(assignments, output)

summary_path <- sub("\\.qs$", ".assignment_summary.tsv", output)
summary <- assignments[, .(
    assignment_rows = .N,
    assigned_reads = uniqueN(readid),
    genes = uniqueN(gene_id),
    transcripts = uniqueN(transcript_id)
), by = type][order(-assigned_reads)]
fwrite(summary, summary_path, sep = "\t")

cat("IsoQuant assignment rows:", nrow(assignments), "\n")
cat("IsoQuant assigned reads:", uniqueN(assignments$readid), "\n")
cat("Saved:", output, "\n")
cat("Summary:", summary_path, "\n")
