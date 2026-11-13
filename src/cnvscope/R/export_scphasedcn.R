# Write the tables of a scphasedcn result object (call_scn(), saved with saveRDS()) to a folder,
# as write_scphasedcn() does, for cnvscope to read:
#
#   Rscript export_scphasedcn.R COLO829BL_call_scn.rds out_dir/
#
# Uses write_scphasedcn() when the scphasedcn package is installed; otherwise writes the tables
# cnvscope needs with data.table alone (cells, calls, regions, thresholds; tumour mode: clones,
# assignments, clone segments) and the run notes (mode, sex, extraction dir).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript export_scphasedcn.R result.rds out_dir")
suppressPackageStartupMessages(library(data.table))
x <- readRDS(args[1])
out <- args[2]
if (!is.list(x) || is.null(x$cells) || is.null(x$mode)) {
  stop("not a scphasedcn result object (call_scn()): no cells / mode")
}
if (requireNamespace("scphasedcn", quietly = TRUE)) {
  invisible(scphasedcn::write_scphasedcn(x, out))
  quit(save = "no")
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(d, name) {
  d <- if (is.null(d)) data.table() else as.data.table(d)
  fwrite(d, file.path(out, name), sep = "\t")
}
calls <- if (is.null(x$calls)) data.table() else as.data.table(x$calls)
if (nrow(calls)) {
  calls <- calls[level != "pool"]
  if (x$mode == "normal") calls <- calls[called %in% TRUE]
}
w(x$cells, "cells.tsv")
w(calls, "calls.tsv")
w(x$regions, "regions.tsv")
w(x$thresholds, "thresholds.tsv")
w(x$events, "events.tsv")
if (x$mode == "tumour") {
  w(x$clones, "clones.tsv")
  w(x$assignments, "assignments.tsv")
  w(x$segments, "clone_segments.tsv")
}
writeLines(sprintf("mode: %s", x$mode), file.path(out, "params.yaml"))
writeLines(c(sprintf("# scphasedcn result (%s mode)", x$mode), "", "## Run notes", "",
             paste0("- ", c(x$run$inputs, x$notes))),
           file.path(out, "NOTES.md"))
