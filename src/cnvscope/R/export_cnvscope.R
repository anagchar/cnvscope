# Write the small file cnvscope opens, beside an ASCAT.sc result.
#
# Run it once at the end of an ASCAT.sc pipeline, where the result object is saved:
#
#   source(system("cnvscope rscript", intern = TRUE))     # or the path of this file
#   save(res, file = "results/result_object.Rda")
#   export_cnvscope(res, "results/result_object.Rda")
#
# This writes results/result_object.cnvscope.h5ad (about 1-2% of the .Rda's size). cnvscope
# looks for it first whenever it is given results/result_object.Rda (or the folder), so the
# large object never has to be read again:
#
#   cnvscope --mode dna --tool ascatsc --input results/ --outdir viewer/
#
# Without the cnvscope command where R runs, a compact results/result_object.cnvscope.json.gz
# is written instead; cnvscope opens that one too.

.cnvscope_here <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f)) {
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    f <- if (length(a)) sub("^--file=", "", a[1]) else NA
  }
  if (is.na(f)) getwd() else dirname(normalizePath(f))
})

export_cnvscope <- function(res, path, cnvscope = Sys.which("cnvscope"), keep_json = FALSE) {
  if (!is.list(res) || is.null(res$allProfiles)) stop("res is not an ASCAT.sc result (no allProfiles)")
  # results/result_object.Rda -> results/result_object.cnvscope.*
  stem <- sub("\\.(rds|rda|rdata)$", "", path, ignore.case = TRUE)
  json <- paste0(stem, ".cnvscope.json")
  if (!exists("export_for_web", mode = "function"))
    source(file.path(.cnvscope_here, "ascatsc_to_web.R"))
  app <- prepare_explorer_data(res)
  export_for_web(app, output_path = json)
  # gzip it: about 5x smaller, and cnvscope reads .json.gz directly
  raw <- readBin(json, "raw", file.size(json))
  con <- gzfile(paste0(json, ".gz"), "wb", compression = 6)
  writeBin(raw, con)
  close(con)
  unlink(json)
  out <- paste0(json, ".gz")
  if (nzchar(cnvscope)) {
    h5 <- paste0(stem, ".cnvscope.h5ad")
    status <- system2(cnvscope, c("convert", shQuote(out), "-o", shQuote(h5)))
    if (status == 0 && file.exists(h5)) {
      if (!keep_json) unlink(out)
      out <- h5
    } else {
      warning("cnvscope convert failed; keeping ", out)
    }
  }
  message(sprintf("cnvscope file: %s (%.1f MB)", out, file.size(out) / 1e6))
  invisible(out)
}
