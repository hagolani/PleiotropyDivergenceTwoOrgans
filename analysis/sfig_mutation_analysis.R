# SFig mutrate (mutation rates) and SFig mutstr (genome structures):
# final fitness per setting under alternative mutational regimes, with the
# corresponding mutational influx and statistics.
#
# Input (in data_dir):
#   boxfit.dat  final fitness, one column per setting x regime,
#               named X<setting>_alt<k> (setting: 1genome, 2genomes, 2genomeshigh)
# Output: sfig_mutrate.pdf, sfig_mutstr.pdf, sfig_mutation_stats_*.csv

library(car)    # leveneTest
library(coin)   # wilcox_test, kruskal_test

source(here::here("analysis", "config.R")); data_dir <- data_proc
settings <- c("1genome", "2genomes", "2genomeshigh")
labels   <- c("With pleiotropy", "Without pleiotropy",
              "Without pleiotropy,\ndouble mutation")
cols     <- c("red", "blue", "purple")
names(labels) <- names(cols) <- settings

# ---- Mutational regimes -------------------------------------------------------
# mu_A: A_constitutive rate; cis and trans rates are mu_A * f_cis and mu_A * f_trans.
# alt6 and alt7 rescale mu_A so that the total influx equals the baseline.
regimes <- data.frame(
  alt     = c("alt0", "alt1", "alt2", "alt6", "alt7"),
  name    = c("Baseline", "High mutation (x5)", "Low mutation (/5)",
              "All equal", "Cis equal"),
  mu_A    = c(0.025, 0.025 * 5, 0.025 / 5, 0.00958, 0.01403),
  f_cis   = c(0.1, 0.1, 0.1, 1, 1),
  f_trans = c(0.05, 0.05, 0.05, 1, 0.05))

n_par <- 21; pop <- 200; gens <- 200
# Number of mutational targets per element type, as used for the influx in the
# paper (21 per element; A_const duplicated without pleiotropy). Adjust if the
# hox_cis (two teeth) or hox_trans (one value per tooth) targets differ.
loci <- rbind(A = c(n_par, 2 * n_par, 2 * n_par), cis = n_par, trans = n_par)
colnames(loci) <- settings
# Without pleiotropy, mu_A is halved to keep the influx of the with-pleiotropy
# setting; the double-mutation setting keeps the full rate.
mu_scale <- c("1genome" = 1, "2genomes" = 0.5, "2genomeshigh" = 1)

influx <- function(r) {
  sapply(settings, function(s) {
    mu <- c(A = r$mu_A * mu_scale[s], cis = r$mu_A * r$f_cis, trans = r$mu_A * r$f_trans)
    pop * gens * loci[, s] * mu
  })
}

# ---- Data ------------------------------------------------------------------
fit <- read.table(file.path(data_dir, "boxfit.dat"))
get_fit <- function(s, alt) na.omit(fit[[paste0("X", s, "_", alt)]])

# ---- Figures ---------------------------------------------------------------
plot_regimes <- function(alts, file) {
  pdf(file, width = 4 * length(alts), height = 8)
  par(mfrow = c(2, length(alts)), mar = c(9, 5, 3, 1))
  ylim <- range(unlist(lapply(alts, function(a) lapply(settings, get_fit, alt = a))))
  for (a in alts) {
    y <- lapply(settings, get_fit, alt = a)
    boxplot(y, names = labels, col = adjustcolor(cols, 0.3), border = cols,
            outline = FALSE, las = 2, ylim = ylim, ylab = "Fitness",
            main = regimes$name[regimes$alt == a])
    for (i in seq_along(y))
      points(jitter(rep(i, length(y[[i]])), 3), y[[i]],
             pch = 16, cex = 0.4, col = adjustcolor(cols[i], 0.4))
    pw <- stats$p_adj[stats$alt == a & stats$test != "Kruskal-Wallis"]
    mtext(sprintf("p (1-2, 1-3, 2-3) = %s", paste(signif(pw, 2), collapse = ", ")),
          side = 3, line = 0.2, cex = 0.6)
  }
  bar_max <- max(sapply(alts, function(a) colSums(influx(regimes[regimes$alt == a, ]))))
  for (a in alts) {
    barplot(influx(regimes[regimes$alt == a, ]), names.arg = labels, las = 2,
            col = c("forestgreen", "tomato", "steelblue"), ylim = c(0, bar_max * 1.1),
            ylab = "Mutational influx (mutations per simulation)",
            legend.text = c("A_constitutive", "hox_cis", "hox_trans"),
            args.legend = list(x = "topleft", bty = "n", cex = 0.8))
  }
  dev.off()
}

# ---- Statistics --------------------------------------------------------------
pairs <- combn(settings, 2, simplify = FALSE)
long  <- function(a) do.call(rbind, lapply(settings, function(s)
  data.frame(alt = a, setting = factor(s, levels = settings), fitness = get_fit(s, a))))

stats <- do.call(rbind, lapply(regimes$alt, function(a) {
  x  <- long(a)
  kw <- kruskal_test(fitness ~ setting, data = x, distribution = "asymptotic")
  pw <- sapply(pairs, function(p)
    pvalue(wilcox_test(fitness ~ setting, data = droplevels(x[x$setting %in% p, ]),
                       distribution = "asymptotic")))
  data.frame(alt = a,
             test  = c("Kruskal-Wallis", sapply(pairs, paste, collapse = " vs ")),
             p     = c(pvalue(kw), pw),
             p_adj = c(NA, p.adjust(pw, "BH")))
}))

assumptions <- do.call(rbind, lapply(regimes$alt, function(a) {
  x   <- long(a)
  lev <- leveneTest(fitness ~ setting, data = x)
  rbind(do.call(rbind, lapply(settings, function(s) {
          t <- shapiro.test(x$fitness[x$setting == s])
          data.frame(alt = a, test = paste("Shapiro-Wilk", s),
                     statistic = t$statistic, p = t$p.value)
        })),
        data.frame(alt = a, test = "Levene", statistic = lev[["F value"]][1],
                   p = lev[["Pr(>F)"]][1]))
}))

plot_regimes(c("alt0", "alt1", "alt2"), "sfig_mutrate.pdf")
plot_regimes(c("alt0", "alt6", "alt7"), "sfig_mutstr.pdf")

write.csv(stats,       "sfig_mutation_stats_tests.csv",       row.names = FALSE)
write.csv(assumptions, "sfig_mutation_stats_assumptions.csv", row.names = FALSE)
write.csv(merge(regimes, data.frame(
  alt = rep(regimes$alt, each = 3), setting = settings,
  influx = unlist(lapply(seq_len(nrow(regimes)), function(i) colSums(influx(regimes[i, ])))))),
  "sfig_mutation_influx.csv", row.names = FALSE)
