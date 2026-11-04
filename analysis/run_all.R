# Runs all analysis scripts. Each script writes its figures and statistics to
# output/<script name>/. Then build the supplementary tables with
#   python analysis/make_tables.py
#
# Usage (from anywhere in the repository):  Rscript analysis/run_all.R
# Scripts that need the raw simulation output (sfig_raw, sfig_S2) are skipped
# if data/raw is empty.

# start from the repository root, so that here::here() finds it
f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(f) == 1) setwd(dirname(dirname(normalizePath(f))))

scripts <- c("fig2_analysis.R", "fig3_analysis.R", "fig4_analysis.R",
             "fig5_analysis.R", "fig6_analysis.R", "sfig_mutation_analysis.R",
             "sfig_raw_analysis.R", "sfig_S2_analysis.R")
needs_raw <- c("sfig_raw_analysis.R", "sfig_S2_analysis.R")
raw_available <- length(list.files(here::here("data", "raw"))) > 1   # more than the README

for (s in scripts) {
  if (s %in% needs_raw && !raw_available) {
    message("Skipping ", s, " (no raw simulation output in data/raw)")
    next
  }
  out <- here::here("output", sub("\\.R$", "", s))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  message("Running ", s, " -> ", out)
  old <- setwd(out)
  tryCatch(source(here::here("analysis", s), local = new.env()),
           error = function(e) message("  FAILED: ", conditionMessage(e)),
           finally = setwd(old))
}
