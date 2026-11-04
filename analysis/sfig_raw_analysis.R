# SFig raw (and the initial-vs-final correlation table): developmental
# parameter values at the first and last generation, per initial condition.
#
# Input (in data_dir/<setting>/alt0/Parameters/):
#   parameters_done<replicate>_evo<initial condition>.txt
#   one row per individual and generation; column 'names' ends with the
#   generation number, the other columns are developmental parameters.
# Output: sfig_raw.pdf (all parameters; the figure shows Egr and Mgr),
#   sfig_raw_correlations.csv

source(here::here("analysis", "config.R")); data_dir <- data_raw
settings  <- c("1genome", "2genomes", "2genomeshigh")
labels    <- c("With pleiotropy", "Without pleiotropy",
               "Without pleiotropy, double mutation")
names(labels) <- settings
min_gens  <- 150    # simulations with fewer recorded generations are skipped
# columns that are not analyzed (unused or invariable parameters)
drop_pattern <- "^names$|^INH_Inhibtion|^Nothing|^Rad_Radius|^Dgr_Downward"

# ---- Mean parameter values at the first and last generation -----------------
# One row per simulation: population mean at generation 1 (initial) and at the
# last generation (final).
read_setting <- function(s) {
  dir   <- file.path(data_dir, s, "alt0", "Parameters")
  files <- list.files(dir, pattern = "^parameters_done\\d+_evo\\d+\\.txt$")
  do.call(rbind, lapply(files, function(f) {
    x <- read.table(file.path(dir, f), header = TRUE)
    x$generation <- as.numeric(sub(".*[^0-9]([0-9]+)[^0-9]*$", "\\1", x$names))
    if (length(unique(x$generation)) < min_gens) return(NULL)
    x <- x[, !grepl(drop_pattern, names(x))]
    params <- setdiff(names(x), "generation")
    rbind(
      data.frame(time = "initial", t(colMeans(x[x$generation == 1, params]))),
      data.frame(time = "final",   t(colMeans(x[x$generation == max(x$generation), params])))
    ) |> transform(setting = s,
                   done = as.numeric(sub(".*done(\\d+)_.*", "\\1", f)),
                   evo  = as.numeric(sub(".*evo(\\d+)\\..*", "\\1", f)))
  }))
}
d <- do.call(rbind, lapply(settings, read_setting))
params <- setdiff(names(d), c("time", "setting", "done", "evo"))

# Initial value per setting and initial condition (shared by its replicates)
i0  <- d[d$time == "initial", ]
ini <- aggregate(i0[, params], by = list(setting = i0$setting, evo = i0$evo), FUN = mean)
fin <- d[d$time == "final", ]

# ---- Correlation between initial and final values ---------------------------
# Per setting and parameter: Pearson correlation across simulations between the
# initial value of the simulation's initial condition and its final value.
cors <- do.call(rbind, lapply(settings, function(s) {
  f  <- fin[fin$setting == s, ]
  i0 <- ini[ini$setting == s, ]
  r <- do.call(rbind, lapply(params, function(p) {
    ct <- cor.test(i0[[p]][match(f$evo, i0$evo)], f[[p]])
    data.frame(setting = labels[s], parameter = p, r = ct$estimate, p = ct$p.value)
  }))
  r$p_adj <- p.adjust(r$p, "BH")    # across parameters within a setting
  r
}))
write.csv(cors, "sfig_raw_correlations.csv", row.names = FALSE)

# ---- Plot: initial to final values per initial condition --------------------
ic      <- sort(unique(d$evo))
ic_cols <- hcl.colors(length(ic), "Dark 3")

pdf("sfig_raw.pdf", width = 12, height = 4)
for (p in params) {
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1))
  ylim <- range(c(ini[[p]], fin[[p]]), na.rm = TRUE)
  for (s in settings) {
    f  <- fin[fin$setting == s, ]
    i0 <- ini[ini$setting == s, ]
    plot(NA, xlim = c(0.8, 2.2), ylim = ylim, xaxt = "n", xlab = "",
         ylab = "Parameter value", main = labels[s])
    axis(1, at = 1:2, labels = c("Initial", "Final"))
    for (k in seq_along(ic)) {
      y0 <- i0[[p]][i0$evo == ic[k]]
      y1 <- sort(f[[p]][f$evo == ic[k]])
      if (length(y1) == 0) next
      polygon(c(1, rep(2, length(y1))), c(y0, y1),
              col = adjustcolor(ic_cols[k], 0.1), border = ic_cols[k])
      points(c(1, rep(2, length(y1))), c(y0, y1), pch = 16, col = ic_cols[k])
    }
    cr <- cors[cors$setting == labels[s] & cors$parameter == p, ]
    legend("topleft", bty = "n", cex = 0.8,
           legend = sprintf("r = %.2f, p = %.2g", cr$r, cr$p))
  }
  legend("bottomright", legend = paste("IC", ic), col = ic_cols, pch = 16,
         bty = "n", cex = 0.7, ncol = 2)
  mtext(p, outer = TRUE, font = 2)
}
dev.off()
