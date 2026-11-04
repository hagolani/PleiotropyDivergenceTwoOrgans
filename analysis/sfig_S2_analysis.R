# SFig S2: adaptive changes per parameter. For each simulation, the proportion
# of 50-generation windows in which a parameter's evolution correlates with
# fitness improvement (r > 0.7), compared between settings at equal fitness.
#
# Input (in data_dir/<setting folder>/):
#   Dataevos/data<done>_evo<evo>_best50.txt        fitness of the top 50%
#                                                  occluding individuals
#   Paras/paras<1|2>_done<done>_evo<evo>_best50.dat parameters of the same
#                                                  individuals (1 upper, 2 lower)
# Output: sfig_S2.pdf, sfig_S2_stats.csv

library(multcomp)   # Tukey contrasts

source(here::here("analysis", "config.R")); data_dir <- data_raw
folders  <- c(with = "1genome/fromstart", without = "2genomes/fromstart2g",
              double = "2genomes/highmuta")
labels   <- c(with = "With pleiotropy", without = "Without pleiotropy",
              double = "Without pleiotropy, double mutation")
window   <- 50
r_min    <- 0.7          # a window is adaptive if r > r_min
fit_lim  <- c(20, 50)    # simulations outside this fitness range are excluded
exclude  <- "INH_Inhibtion_of_SHH_over_BMP4"

# ---- Rolling-window correlations per simulation ------------------------------
roll_cor <- function(x, y) sapply(seq_len(length(x) - window + 1), function(i) {
  k <- i:(i + window - 1); r <- suppressWarnings(cor(x[k], y[k]))
  ifelse(is.na(r), 0, r)
})

one_sim <- function(s, done, evo) {
  base <- file.path(data_dir, folders[s])
  f_fit <- file.path(base, "Dataevos", sprintf("data%d_evo%d_best50.txt", done, evo))
  f_par <- file.path(base, "Paras", sprintf("paras%d_done%d_evo%d_best50.dat", 1:2, done, evo))
  if (!all(file.exists(c(f_fit, f_par)))) return(NULL)

  fit <- read.table(f_fit)
  names(fit)[5:6] <- c("fitness", "generation")
  fit$fitness <- 1 / fit$fitness                       # stored as inverse fitness
  mean_fit    <- tapply(fit$fitness, fit$generation, mean)
  final_fit   <- mean(head(sort(fit$fitness, decreasing = TRUE), 10))  # 10 best

  do.call(rbind, lapply(1:2, function(tooth) {
    p   <- read.table(f_par[tooth])
    gen <- p$generation
    p   <- p[, setdiff(names(p)[seq_len(ncol(p) - 5)], "generation"), drop = FALSE]
    p   <- p[, !grepl(exclude, names(p)), drop = FALSE]
    mean_par <- aggregate(p, by = list(generation = gen), FUN = mean)
    n   <- min(nrow(mean_par), length(mean_fit))
    do.call(rbind, lapply(names(p), function(v) {
      r <- roll_cor(mean_fit[seq_len(n)], mean_par[[v]][seq_len(n)])
      data.frame(setting = s, done = done, evo = evo, fitness = final_fit,
                 tooth = c("upper", "lower")[tooth],
                 element = ifelse(v == "hn", "hox_trans",
                                  ifelse(grepl("hox", v), "hox_cis", "A_constitutive")),
                 parameter = sub("_hox$", "", v),
                 adaptive = sum(r > r_min), windows = length(r))
    }))
  }))
}

d <- do.call(rbind, lapply(names(folders), function(s)
  do.call(rbind, lapply(1:12, function(done)
    do.call(rbind, lapply(1:10, function(evo) one_sim(s, done, evo)))))))
if (is.null(d)) stop("No simulations found: check data_dir (", normalizePath(data_dir, mustWork = FALSE),
                     ") and the folder names in 'folders'.")
message("Simulations read per setting: ",
        paste(names(folders), table(factor(unique(d[, c("setting", "done", "evo")])$setting,
                                           levels = names(folders))), collapse = ", "),
        "; final fitness range: ", paste(round(range(d$fitness), 2), collapse = " - "))
d <- d[d$fitness > fit_lim[1] & d$fitness < fit_lim[2], ]
if (nrow(d) == 0) stop("All simulations fall outside fit_lim; check the fitness scale.")
d$setting <- factor(d$setting, levels = names(folders))

# ---- Comparison between settings at equal fitness ----------------------------
# Effect size: binomial GLM of adaptive windows ~ fitness + setting, predicted
# proportion at the lowest fitness. Significance: Tukey contrasts in a linear
# model of the number of adaptive windows ~ setting + fitness.
groups <- unique(d[, c("element", "tooth", "parameter")])
stats <- do.call(rbind, lapply(seq_len(nrow(groups)), function(g) {
  x  <- merge(d, groups[g, ])
  bm <- glm(cbind(adaptive, windows - adaptive) ~ fitness + setting,
            family = binomial, data = x)
  pr <- predict(bm, data.frame(fitness = min(x$fitness), setting = levels(d$setting)),
                type = "response")
  tk <- tryCatch(summary(glht(lm(adaptive ~ setting + fitness, data = x),
                              linfct = mcp(setting = "Tukey")))$test$pvalues,
                 error = function(e) rep(NA, 3))   # e.g. no adaptive windows at all
  cbind(groups[g, ],
        delta_with_without = pr[1] - pr[2], delta_with_double = pr[1] - pr[3],
        p_with_without = tk[1], p_with_double = tk[2], p_without_double = tk[3])
}))
rownames(stats) <- NULL
write.csv(stats, "sfig_S2_stats.csv", row.names = FALSE)

# ---- SFig S2: difference in adaptive windows per parameter -------------------
pdf("sfig_S2.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(12, 5, 3, 1))
for (e in c("A_constitutive", "hox_cis")) {
  s <- stats[stats$element == e, ]
  s <- s[order(s$parameter, s$tooth), ]
  barplot(s$delta_with_without, las = 2, cex.names = 0.6,
          names.arg = paste(sub("_.*", "", s$parameter), s$tooth),
          col = ifelse(is.na(s$p_with_without) | s$p_with_without >= 0.05, "grey",
                       ifelse(s$tooth == "upper", "darkgreen", "orange")),
          ylab = "Difference in proportion of adaptive windows\n(with - without pleiotropy)",
          main = e)
  abline(h = 0)
}
dev.off()
