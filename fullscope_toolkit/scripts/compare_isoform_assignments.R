#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript compare_isoform_assignments.R bambu.qs isoquant.qs output.tsv [sample]")
}

suppressPackageStartupMessages({
    library(data.table)
    library(qs)
})

summarize_assignments <- function(path, annotator, sample) {
    x <- as.data.table(qread(path))
    required <- c("readid", "gene_id", "transcript_id")
    missing_columns <- setdiff(required, names(x))
    if (length(missing_columns) > 0) {
        stop(path, " is missing columns: ", paste(missing_columns, collapse = ", "))
    }
    per_read <- x[, .(
        transcript_candidates = uniqueN(transcript_id, na.rm = TRUE),
        gene_candidates = uniqueN(gene_id, na.rm = TRUE)
    ), by = readid]
    data.table(
        sample = sample,
        annotator = annotator,
        assignment_rows = nrow(x),
        annotated_reads = nrow(per_read),
        unique_transcript_reads = sum(per_read$transcript_candidates == 1L),
        ambiguous_transcript_reads = sum(per_read$transcript_candidates > 1L),
        ambiguous_transcript_fraction = mean(per_read$transcript_candidates > 1L),
        unique_gene_reads = sum(per_read$gene_candidates == 1L),
        ambiguous_gene_reads = sum(per_read$gene_candidates > 1L),
        genes = uniqueN(x$gene_id, na.rm = TRUE),
        transcripts = uniqueN(x$transcript_id, na.rm = TRUE)
    )
}

sample <- if (length(args) >= 4) args[4] else NA_character_
result <- rbindlist(list(
    summarize_assignments(args[1], "bambu", sample),
    summarize_assignments(args[2], "isoquant", sample)
), use.names = TRUE)
fwrite(result, args[3], sep = "\t")
print(result)
