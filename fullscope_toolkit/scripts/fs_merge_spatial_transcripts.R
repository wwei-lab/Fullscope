args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript fs_merge_spatial_transcripts.R <transdf_qs> <cidmap_tsv> <merged_all_qs> <merged_unique_qs> [gtf]")
}

transdf_qs <- args[1]
cidmap_tsv <- args[2]
merged_all_qs <- args[3]
merged_unique_qs <- args[4]
gtf_path <- if (length(args) >= 5) args[5] else NA_character_

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
  library(rtracklayer)
  library(dplyr)
})

MySplit <- function(str, sep, n) {
  unlist(lapply(strsplit(str, sep, fixed = TRUE), `[[`, n))
}

transdf <- qread(transdf_qs)
cidmap <- fread(cidmap_tsv, data.table = FALSE)

if (!all(c("readid", "gene_id", "transcript_id") %in% colnames(transdf))) {
  stop("transdf must contain readid, gene_id, transcript_id")
}
if (!all(c("readid", "cidPos") %in% colnames(cidmap))) {
  stop("cidmap must contain readid, cidPos")
}

if (!all(c("gene_name", "transcript_name") %in% colnames(transdf))) {
  if (is.na(gtf_path) || !file.exists(gtf_path)) {
    stop("A readable GTF path is required when gene_name or transcript_name is absent")
  }
  gtf <- import(gtf_path)
  anno <- as.data.frame(mcols(gtf)) %>%
    select(gene_id, gene_name, transcript_id, transcript_name) %>%
    filter(!is.na(gene_id), !is.na(transcript_id)) %>%
    distinct(gene_id, transcript_id, .keep_all = TRUE)
  transdf <- transdf %>% left_join(anno, by = c("gene_id", "transcript_id"))
}

cidmap_merged <- merge(cidmap, transdf, by = "readid")
cidmap_unique <- cidmap_merged %>%
  group_by(readid) %>%
  filter(dplyr::n_distinct(transcript_id, na.rm = TRUE) == 1) %>%
  ungroup()

cidmap_merged$x <- as.numeric(MySplit(cidmap_merged$cidPos, "_", 1))
cidmap_merged$y <- as.numeric(MySplit(cidmap_merged$cidPos, "_", 2))
cidmap_unique$x <- as.numeric(MySplit(cidmap_unique$cidPos, "_", 1))
cidmap_unique$y <- as.numeric(MySplit(cidmap_unique$cidPos, "_", 2))

qsave(cidmap_merged, file = merged_all_qs)
qsave(cidmap_unique, file = merged_unique_qs)

cat("Merged readids:", length(unique(cidmap_merged$readid)), "\n")
cat("Unique-transcript readids:", length(unique(cidmap_unique$readid)), "\n")
