# Figure 3A: correlation of cusp change between upper and lower teeth versus
# final fitness, per setting, plus the corresponding statistics.
#
# Input (in data_dir):
#   alldata_monoTooth.dat  one row per simulation, with columns
#                          exp (1-3), evo (initial condition), fitness, corr_ul
# Output (in the working directory):
#   fig3A_scatter.pdf, sfig_corrIC.pdf, fig3_stats_*.csv

library(car)    # leveneTest
library(coin)   # wilcox_test, kruskal_test



source(here::here("analysis", "config.R")); data_dir <- data_proc
settings <- c("1", "2", "3")
labels   <- c("With pleiotropy", "Without pleiotropy",
              "Without pleiotropy, double mutation")
cols     <- c("firebrick", "steelblue", "purple4")
shapes   <- c(16, 17, 8)
names(labels) <- names(cols) <- names(shapes) <- settings

cor_method <- "spearman"   # correlation reported in Fig 3A and Table S4

# ---- Data ------------------------------------------------------------------
d <- na.omit(read.table(file.path(data_dir, "alldata_monoTooth.dat")))
d <- d[d$exp %in% as.numeric(settings), ]
d$exp <- factor(d$exp, levels = settings)
d$ic  <- factor(d$evo)       # initial condition

vars  <- c("fitness", "corr_ul")
pairs <- combn(settings, 2, simplify = FALSE)

# ---- Fig 3A: scatter plot per setting --------------------------------------
pdf("fig3A_scatter.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1))
for (s in settings) {
  x  <- d[d$exp == s, ]
  ct <- cor.test(x$fitness, x$corr_ul, method = cor_method, exact = FALSE)
  plot(x$fitness, x$corr_ul, pch = shapes[s], col = cols[s], cex = 0.8,
       xlab = "Fitness", ylab = "Correlation of cusp change (upper vs lower)",
       main = labels[s], xlim = range(d$fitness), ylim = range(d$corr_ul))
  abline(lm(corr_ul ~ fitness, data = x), col = cols[s], lwd = 2)
  legend("topleft", bty = "n", text.col = cols[s],
         legend = sprintf("r_s = %.2f\np = %.2g", ct$estimate, ct$p.value))
}
dev.off()

# ---- Statistics --------------------------------------------------------------
# Correlation between fitness and corr_ul per setting (Table S4, Fig 3)
correlation <- do.call(rbind, lapply(settings, function(s) {
  x  <- d[d$exp == s, ]
  ct <- cor.test(x$fitness, x$corr_ul, method = cor_method, exact = FALSE)
  data.frame(setting = labels[s], n = nrow(x), r = ct$estimate, p = ct$p.value)
}))

# Normality per setting (Shapiro-Wilk) and homoscedasticity (Levene)
shapiro <- do.call(rbind, lapply(vars, function(v)
  do.call(rbind, lapply(settings, function(s) {
    t <- shapiro.test(d[d$exp == s, v])
    data.frame(variable = v, setting = labels[s], W = t$statistic, p = t$p.value)
  }))))
levene <- do.call(rbind, lapply(vars, function(v) {
  t <- leveneTest(as.formula(paste(v, "~ exp")), data = d)
  data.frame(variable = v, F = t[["F value"]][1], p = t[["Pr(>F)"]][1])
}))

# Differences between settings: Kruskal-Wallis + pairwise Wilcoxon (BH)
between <- do.call(rbind, lapply(vars, function(v) {
  kw <- kruskal_test(as.formula(paste(v, "~ exp")), data = d,
                     distribution = "asymptotic")
  pw <- sapply(pairs, function(p) {
    s <- droplevels(d[d$exp %in% p, ])
    pvalue(wilcox_test(as.formula(paste(v, "~ exp")), data = s,
                       distribution = "asymptotic"))
  })
  data.frame(variable = v,
             test = c("Kruskal-Wallis",
                      sapply(pairs, function(p) paste(labels[p], collapse = " vs "))),
             p = c(pvalue(kw), pw), p_adj = c(NA, p.adjust(pw, "BH")))
}))

# Effect of initial conditions within each setting (SFig corrIC, Table S1):
#   F-test of lm(fitness ~ corr_ul) vs lm(fitness ~ corr_ul + IC), and
#   Kruskal-Wallis of corr_ul ~ IC
ic_effect <- do.call(rbind, lapply(settings, function(s) {
  x  <- droplevels(d[d$exp == s, ])
  a  <- anova(lm(fitness ~ corr_ul, data = x), lm(fitness ~ corr_ul + ic, data = x))
  kw <- kruskal_test(corr_ul ~ ic, data = x, distribution = "asymptotic")
  data.frame(setting = labels[s],
             F_ic = a$F[2], p_ic = a[["Pr(>F)"]][2],
             KW_chi2_corr = statistic(kw), KW_p_corr = pvalue(kw))
}))

# ---- SFig corrIC: Fig 3A per setting, coloured by initial condition ----------
ic_levels <- levels(d$ic)
ic_cols   <- adjustcolor(hcl.colors(length(ic_levels), "Dark 3"), alpha.f = 0.7)
names(ic_cols) <- ic_levels

pdf("sfig_corrIC.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1))
for (s in settings) {
  x <- d[d$exp == s, ]
  plot(x$fitness, x$corr_ul, pch = 16, cex = 1.2, col = ic_cols[as.character(x$ic)],
       xlab = "Fitness", ylab = "Correlation of cusp change (upper vs lower)",
       main = labels[s], xlim = range(d$fitness), ylim = range(d$corr_ul))
  for (i in ic_levels) {                       # convex hull per initial condition
    xi <- x[x$ic == i, ]
    if (nrow(xi) >= 3) {
      h <- chull(xi$fitness, xi$corr_ul)
      polygon(xi$fitness[h], xi$corr_ul[h], border = ic_cols[i], lty = 2, lwd = 1.5)
    }
  }
  p <- ic_effect$p_ic[ic_effect$setting == labels[s]]
  legend("topleft", bty = "n", legend = sprintf("IC effect: p = %.2g", p))
}
legend("bottomright", legend = paste("IC", ic_levels), col = ic_cols,
       pch = 16, bty = "n", cex = 0.7, ncol = 2)
dev.off()

write.csv(correlation, "fig3_stats_correlation.csv", row.names = FALSE)
write.csv(shapiro,     "fig3_stats_shapiro.csv",     row.names = FALSE)
write.csv(levene,      "fig3_stats_levene.csv",      row.names = FALSE)
write.csv(between,     "fig3_stats_between_settings.csv", row.names = FALSE)
write.csv(ic_effect,   "fig3_stats_initial_conditions.csv", row.names = FALSE)
