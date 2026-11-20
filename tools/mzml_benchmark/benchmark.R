# Pass CLI arguments as a character vector, for example c("simulate", "calibration", "-o", "synthetic").
mzml_benchmark <- function(args, script, python = "python") {
    output <- system2(python, shQuote(c(normalizePath(script, mustWork = TRUE), args)),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
        stop(paste(output, collapse = "\n"), call. = FALSE)
    }
    invisible(output)
}
