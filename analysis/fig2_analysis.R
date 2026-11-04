# Figure 2: fitness, morphological change (EMD) and cusp number in the
# three settings, plus the corresponding statistics.
#
# Input (in data_dir):
#   allexp.dat         one row per individual: the 10 best occluding
#                      individuals of the last generation of each simulation
#   allcuspsallexp.dat cusp numbers of upper and lower teeth (heatmaps)
# Output (in the working directory):
#   fig2_boxplots.pdf, fig2_cusp_number.pdf, fig2_cusp_heatmaps.pdf, sfig_IC.pdf
#   fig2_stats_*.csv

library(car)      # leveneTest
library(coin)     # wilcox_test, kruskal_test
library(ggplot2)
library(ggpubr)   # ggarrange

set.seed(123)     # heatmap subsampling

source(here::here("analysis", "config.R")); data_dir <- data_proc
settings <- c("1genome_alt0", "2genomes_alt0", "2genomeshigh_alt0")
labels   <- c("With pleiotropy", "Without pleiotropy",
              "Without pleiotropy,\ndouble mutation")
cols     <- c("red", "blue", "purple")
names(labels) <- names(cols) <- settings

# ---- Data ------------------------------------------------------------------
d <- read.table(file.path(data_dir, "allexp.dat"))
d <- d[d$exp %in% settings, ]
d$exp      <- factor(d$exp, levels = settings)
d$ic       <- factor(d$evo)          # initial condition
d$emd_up   <- d[, 4]
d$emd_low  <- d[, 5]
d$cusp_up  <- d[, 6]
d$cusp_low <- d[, 7]
d$fitness  <- d[, 8]
d$cusp_all <- d$cusp_up + d$cusp_low

vars <- c(fitness  = "Fitness (occlusion)",
          emd_up   = "EMD upper tooth",
          emd_low  = "EMD lower tooth",
          cusp_all = "Cusp number (upper + lower)")

pairs <- combn(settings, 2, simplify = FALSE)

# Pairwise Wilcoxon rank-sum tests between settings (BH-adjusted)
pairwise_wilcox <- function(v) {
  res <- do.call(rbind, lapply(pairs, function(p) {
    s  <- droplevels(d[d$exp %in% p, ])
    wt <- wilcox_test(as.formula(paste(v, "~ exp")), data = s,
                      distribution = "asymptotic")
    data.frame(variable = v, comparison = paste(p, collapse = " vs "),
               Z = round(statistic(wt), 3), p = pvalue(wt))
  }))
  res$p_adj <- p.adjust(res$p, method = "BH")
  res
}

# ---- Panel A: boxplots of fitness and EMD ----------------------------------
boxpanel <- function(v, ylab, main = "") {
  y <- d[[v]]
  n <- tapply(!is.na(y), d$exp, sum)
  boxplot(y ~ d$exp, names = paste0(labels, " (n = ", n, ")"), las = 2,
          ylab = ylab, main = main, xlab = "", outline = FALSE,
          ylim = range(y, na.rm = TRUE) * c(0.9, 1.2))
  for (i in seq_along(settings)) {
    yi <- y[d$exp == settings[i]]
    points(jitter(rep(i, length(yi)), 3), yi,
           col = adjustcolor(cols[i], 0.5), pch = 16, cex = 0.5)
  }
  # significance brackets
  pw  <- pairwise_wilcox(v)
  top <- max(y, na.rm = TRUE)
  for (k in seq_along(pairs)) {
    x  <- match(pairs[[k]], settings)
    yb <- top * (1 + 0.05 * k)
    segments(x[1], yb, x[2], yb)
    text(mean(x), yb, signif(pw$p_adj[k], 2), pos = 3, cex = 0.8,
         col = ifelse(pw$p_adj[k] < 0.05, "red", "black"))
  }
}

pdf("fig2_boxplots.pdf", width = 10, height = 5)
par(mfrow = c(1, 3), mar = c(14, 5, 3, 1), cex.axis = 0.8)
boxpanel("fitness", "Fitness (occlusion)")
boxpanel("emd_up",  "Morphological distance to initial tooth (EMD)", "Upper teeth")
boxpanel("emd_low", "Morphological distance to initial tooth (EMD)", "Lower teeth")
dev.off()

# ---- Panel B: cusp number distributions ------------------------------------
cusp_props <- function(v) {
  lv <- seq_len(max(d[[v]], na.rm = TRUE))
  sapply(settings, function(s)
    prop.table(table(factor(d[d$exp == s, v], levels = lv))))
}

pdf("fig2_cusp_number.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
for (v in c("cusp_up", "cusp_low")) {
  barplot(t(cusp_props(v)), beside = TRUE, col = cols, ylim = c(0, 1),
          xlab = "Number of cusps", ylab = "Frequency",
          main = ifelse(v == "cusp_up", "Upper teeth", "Lower teeth"))
}
legend("topright", legend = labels, fill = cols, bty = "n", cex = 0.7)
dev.off()

# ---- Panel C: upper vs lower cusp number, all teeth, all generations -------
cusps   <- read.table(file.path(data_dir, "allcuspsallexp.dat"))
max_c   <- max(cusps$cuspUp, cusps$cuspLow, na.rm = TRUE)
n_sample <- 40000                     # random subsample per setting

heatmap_panel <- function(s, title) {
  x  <- cusps[cusps$exp == s, ]
  x  <- x[sample(nrow(x), min(n_sample, nrow(x))), ]
  hm <- as.data.frame(table(up  = factor(x$cuspUp,  levels = 1:max_c),
                            low = factor(x$cuspLow, levels = 1:max_c)))
  hm$Freq[hm$Freq == 0] <- NA        # empty cells left blank on the log scale
  main <- ggplot(hm, aes(up, low, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "log",
                        na.value = "grey90",
                        labels = scales::number_format(accuracy = 1)) +
    labs(title = title, x = "Upper tooth: number of cusps",
         y = "Lower tooth: number of cusps", fill = "Teeth")
  marg_x <- ggplot(aggregate(Freq ~ up,  hm, sum), aes(up,  Freq)) +
    geom_col() + labs(x = NULL, y = NULL)
  marg_y <- ggplot(aggregate(Freq ~ low, hm, sum), aes(low, Freq)) +
    geom_col() + coord_flip() + labs(x = NULL, y = NULL)
  ggarrange(marg_x, NULL, main + theme(legend.position = "left"), marg_y,
            nrow = 2, ncol = 2, widths = c(3, 1), heights = c(1, 3))
}

pdf("fig2_cusp_heatmaps.pdf", width = 6, height = 6)
for (s in settings) print(heatmap_panel(s, gsub("\n", " ", labels[s])))
dev.off()

# ---- Statistics --------------------------------------------------------------
# Normality per setting (Shapiro-Wilk)
shapiro <- do.call(rbind, lapply(names(vars), function(v)
  do.call(rbind, lapply(settings, function(s) {
    t <- shapiro.test(d[d$exp == s, v])
    data.frame(variable = v, setting = s, W = t$statistic, p = t$p.value)
  }))))

# Homoscedasticity across settings (Levene)
levene <- do.call(rbind, lapply(names(vars), function(v) {
  t <- leveneTest(as.formula(paste(v, "~ exp")), data = d)
  data.frame(variable = v, F = t[["F value"]][1], p = t[["Pr(>F)"]][1])
}))

# Kruskal-Wallis across settings + pairwise Wilcoxon (reported in Table S4)
kw <- do.call(rbind, lapply(names(vars), function(v) {
  t <- kruskal_test(as.formula(paste(v, "~ exp")), data = d,
                    distribution = "asymptotic")
  data.frame(variable = v, chi2 = statistic(t), p = pvalue(t))
}))
wilcox <- do.call(rbind, lapply(names(vars), pairwise_wilcox))

# ---- SFig IC and Table S1: effect of initial conditions ---------------------
# Unit of the test: "simulation" (mean of the 10 best individuals of each
# simulation, i.e. 10 values per initial condition) or "individual" (all rows).
ic_unit <- "simulation"
ic_data <- if (ic_unit == "simulation") {
  aggregate(cbind(fitness, cusp_all) ~ exp + ic + done, data = d, FUN = mean)
} else d
ic_vars <- c(cusp_all = "Number of cusps (upper + lower)", fitness = "Fitness")

ic_effect <- do.call(rbind, lapply(settings, function(s) {
  x <- droplevels(ic_data[ic_data$exp == s, ])
  do.call(rbind, lapply(names(ic_vars), function(v) {
    t <- kruskal_test(as.formula(paste(v, "~ ic")), data = x,
                      distribution = "asymptotic")
    data.frame(setting = s, variable = v, chi2 = statistic(t),
               df = nlevels(x$ic) - 1, p = pvalue(t))
  }))
}))

# boxplots by initial condition: rows = variables, columns = settings
pdf("sfig_IC.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(5, 5, 3, 1))
for (v in names(ic_vars)) for (s in settings) {
  x <- ic_data[ic_data$exp == s, ]
  boxplot(x[[v]] ~ x$ic, outline = FALSE, xlab = "Initial condition",
          ylab = ic_vars[v], main = gsub("\n", " ", labels[s]),
          ylim = range(ic_data[[v]], na.rm = TRUE))
  points(jitter(as.numeric(x$ic), 0.7), x[[v]],
         col = adjustcolor(cols[s], 0.5), pch = 16, cex = 0.6)
  p <- ic_effect$p[ic_effect$setting == s & ic_effect$variable == v]
  legend("topleft", bty = "n", legend = sprintf("Kruskal-Wallis p = %.2g", p))
}
dev.off()

write.csv(shapiro,   "fig2_stats_shapiro.csv",   row.names = FALSE)
write.csv(levene,    "fig2_stats_levene.csv",    row.names = FALSE)
write.csv(kw,        "fig2_stats_kruskal.csv",   row.names = FALSE)
write.csv(wilcox,    "fig2_stats_wilcoxon.csv",  row.names = FALSE)
write.csv(ic_effect, "fig2_stats_initial_conditions.csv", row.names = FALSE)
