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
source(file.path(script_dir, "..", "..", "dependent model", "likelihood for other models", "alpha-beta_likelihood.R"))

wrapper <- function(file_path, alpha, beta) {
    f_mat <- as.matrix(read.csv(file_path, header = FALSE))
    prob <- ab_likelihood(f_mat, alpha, beta)
    outpath <- paste0(file_path, ".R.csv")
    write(prob, outpath)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    cat("Usage: Rscript alpha_beta_likelihood.R <file_path> <alpha> <beta>\n")
    cat("  file_path: csv file location of the F-matrix\n")
    cat("  alpha: alpha of the Beta(alpha, beta)")
    cat("  beta: beta of the Beta(alpha, beta)")
    quit(status = 1)
}
file_path <- as.character(args[1])
alpha <- as.double(args[2])
beta <- as.double(args[3])
if (is.na(alpha) || is.na(beta) || alpha <= 0 || beta <= 0) {
    cat("Error: alpha and beta must be positive\n")
    quit(status = 1)
}

# Run the script
wrapper(file_path, alpha, beta)
