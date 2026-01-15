suppressPackageStartupMessages({
    suppressMessages({
        # Get this script's directory
        this_file <- (function() {
            cmdArgs <- commandArgs(trailingOnly = FALSE)
            needle <- "--file="
            match <- grep(needle, cmdArgs)
            if (length(match) > 0) {
                return(normalizePath(sub(needle, "", cmdArgs[match])))
            } else {
                return(normalizePath(sys.frames()[[1]]$ofile))
            }
        })()
        script_dir <- dirname(this_file)

        # Load required functions from enumerate_clean.R
        # Redirect output to suppress messages
        sink(tempfile())
        source(file.path(script_dir, "..", "enumerate_clean.R"))
        sink()
    })
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 2) {
    cat("Usage: Rscript enumerate_fmatrices.R <n> <outpath>\n")
    cat("  n: number of tips\n")
    cat("  outpath: file location where results will be written\n")
    quit(status = 1)
}

# Get arguments
n <- as.integer(args[1])
outpath <- args[2]

# Validate n
if (is.na(n) || n < 2) {
    cat("Error: n must be an integer >= 2\n")
    quit(status = 1)
}

# Run the functions silently
suppressMessages({
    out <- construct_trees_het_D(n)
    F.list <- Fmat_from_construct_het(out)
})

# Format output as Python-style list of lists of lists
format_matrix_to_python <- function(mat) {
    rows <- apply(mat, 1, function(row) {
        paste0("[", paste(row, collapse = ", "), "]")
    })
    paste0("[", paste(rows, collapse = ", "), "]")
}

# Create the full output string
python_list <- paste0(
    "[",
    paste(sapply(F.list, format_matrix_to_python), collapse = ", "),
    "]"
)

# Write to file silently
writeLines(python_list, outpath)
