# Export an ASCAT.sc result (.rds, or .Rda/.RData holding one) to the web JSON format
# that cnvscope's ascatsc_json reader reads. Uses ascatsc_to_web.R from the ASCAT.sc Explorer.
#
#   Rscript export_ascatsc.R <result.rds|.Rda> <out.json>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: export_ascatsc.R <result.rds|.Rda> <out.json>")
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "ascatsc_to_web.R"))

path <- args[1]
is_result <- function(x) is.list(x) && !is.null(x$allProfiles)
if (grepl("\\.(rda|rdata)$", tolower(path))) {
  env <- new.env()
  message(sprintf("Loading %s ...", path))
  objs <- load(path, envir = env)
  hits <- Filter(function(n) is_result(get(n, envir = env)), objs)
  if (!length(hits)) stop(sprintf("no ASCAT.sc result (a list with allProfiles) among: %s", paste(objs, collapse = ", ")))
  message(sprintf("Using object '%s'", hits[1]))
  res <- get(hits[1], envir = env)
  rm(env)
} else {
  message(sprintf("Loading %s ...", path))
  res <- readRDS(path)
}
if (!is_result(res)) stop("this is not an ASCAT.sc result (no allProfiles)")
app_data <- prepare_explorer_data(res)
export_for_web(app_data, output_path = args[2])
