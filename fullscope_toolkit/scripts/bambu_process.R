#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript bambu_process.R annotations.gtf genome.fa input.bam output_prefix sample")
}

# The packaged R 4.5 environment uses Bioconductor 3.22. Avoid online version
# discovery on compute nodes; both values remain user-overridable.
if (!nzchar(Sys.getenv("R_BIOC_VERSION"))) {
    Sys.setenv(R_BIOC_VERSION = "3.22")
}
if (!nzchar(Sys.getenv("BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS"))) {
    Sys.setenv(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "FALSE")
}

suppressPackageStartupMessages({
    library(bambu)
    library(qs)
})

out_table_from_bambu <- function(read_transcripts, sample) {
    keep <- (lengths(read_transcripts$compatibleMatches) +
             lengths(read_transcripts$equalMatches)) > 0
    read_transcripts <- read_transcripts[keep, ]

    equal_count <- lengths(read_transcripts$equalMatches)
    compatible_count <- lengths(read_transcripts$compatibleMatches)

    equal <- data.frame(
        readid = rep(read_transcripts$readId, equal_count),
        tranid = unlist(read_transcripts$equalMatches),
        type = "equal"
    )
    compatible <- data.frame(
        readid = rep(read_transcripts$readId, compatible_count),
        tranid = unlist(read_transcripts$compatibleMatches),
        type = "compatible"
    )

    result <- unique(rbind(equal, compatible))
    result$sample <- sample
    result
}

bambu_matrix_build <- function(bam, annotations, genome, output, sample) {
    workers <- suppressWarnings(
        as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
    )
    if (is.na(workers) || workers < 1L) {
        workers <- 1L
    }

    se <- bambu(
        reads = bam,
        annotations = annotations,
        genome = genome,
        discovery = FALSE,
        NDR = 0,
        trackReads = TRUE,
        ncore = workers,
        lowMemory = TRUE
    )
    qsave(se, paste0(output, "_se.qs"))

    maps <- lapply(
        metadata(se)$readToTranscriptMaps,
        out_table_from_bambu,
        sample = sample
    )
    read_transcripts <- do.call(rbind, maps)

    annotation <- as.data.frame(rowData(se))[, c("TXNAME", "GENEID", "txid")]
    colnames(annotation) <- c("transcript_id", "gene_id", "txid")
    match_id <- match(read_transcripts$tranid, annotation$txid)
    read_transcripts$transcript_id <- annotation$transcript_id[match_id]
    read_transcripts$gene_id <- annotation$gene_id[match_id]

    message("Annotated reads: ", length(unique(read_transcripts$readid)))
    qsave(read_transcripts, paste0(output, "_trans_total_anno.qs"))
    read_transcripts
}

gtf <- args[1]
genome <- args[2]
bam <- args[3]
output <- args[4]
sample <- args[5]

message("Starting Bambu with offline Bioconductor version ", Sys.getenv("R_BIOC_VERSION"))
annotations <- prepareAnnotations(gtf)
bambu_matrix_build(bam, annotations, genome, output, sample)
