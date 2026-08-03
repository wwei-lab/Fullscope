#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript test_isoquant_conversion.R converter.R fixture.tsv")
}

output <- tempfile(fileext = ".qs")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(args[1], args[2], output, "test_sample"))
stopifnot(status == 0L)

suppressPackageStartupMessages(library(qs))
x <- qread(output)
stopifnot(nrow(x) == 3L)
stopifnot(length(unique(x$readid)) == 2L)
stopifnot(sum(x$readid == "read_ambiguous") == 2L)
stopifnot(all(x$sample == "test_sample"))
stopifnot(all(x$assignment_source == "isoquant"))
cat("IsoQuant conversion test passed\n")
