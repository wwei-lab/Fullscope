args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
    stop("Provide one or more R files to parse.", call. = FALSE)
}

for (path in args) {
    parse(file = path)
    message("Parsed: ", path)
}

message("R syntax test passed.")
