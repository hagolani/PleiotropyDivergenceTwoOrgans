# Figure 6 (and SFig S4): proportion of single mutations that improve fitness
# and are passed on to the next generation, per genetic element and setting.
#
# Input (in data_dir), one file per setting (exp1 = with pleiotropy,
# exp2 = without pleiotropy, exp3 = without pleiotropy, double mutation):
#   ori_norm_allevos_propsWFIT_exp{1,2,3}.dat
#   one row per simulation; columns X<t>.<c>.<a> give the proportion of
#   adaptive mutations among mutation events with t hox_trans, c hox_cis and
#   a A_constitutive mutations (t = 1/2: upper/lower hox_trans).
#   Only single-mutation events are analyzed.
# Output: fig6_boxplots.pdf (with pleiotropy vs without),
#   sfig_S4_boxplots.pdf (all three settings), fig6_stats.csv

library(ggplot2)

source(here::here("analysis", "config.R")); data_dir <- data_proc
settings <- c(with = 1, without = 2, double = 3)
labels   <- c(with = "With pleiotropy", without = "Without pleiotropy",
              double = "Without pleiotropy, double mutation")
cols     <- c(with = "red", without = "blue", double = "purple")
elements <- c("A_constitutive", "hox_cis", "hox_trans")

# ---- Data (long format: one row per simulation x element) -------------------
d <- do.call(rbind, lapply(names(settings), function(s) {
  x <- read.table(file.path(data_dir,
         sprintf("ori_norm_allevos_propsWFIT_exp%d.dat", settings[s])))
  id <- seq_len(nrow(x))
  rbind(
    data.frame(sim = id, setting = s, element = "A_constitutive", prop = x$X0.0.1),
    data.frame(sim = id, setting = s, element = "hox_cis",        prop = x$X0.1.0),
    # upper and lower hox_trans combined
    data.frame(sim = id, setting = s, element = "hox_trans",
               prop = c(x$X1.0.0, x$X2.0.0)))
}))
d$setting <- factor(d$setting, levels = names(settings))
d$element <- factor(d$element, levels = elements)

# ---- Plots ------------------------------------------------------------------
plot_props <- function(sub, file) {
  p <- ggplot(sub, aes(element, prop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = setting), width = 0.2, height = 0, alpha = 0.5) +
    facet_wrap(~ setting, labeller = as_labeller(labels)) +
    scale_color_manual(values = cols, guide = "none") +
    labs(x = NULL, y = "Adaptive mutations / total number of mutations") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file, p, width = 3 * nlevels(droplevels(sub$setting)), height = 5)
}
plot_props(droplevels(subset(d, setting != "double")), "fig6_boxplots.pdf")
plot_props(d, "sfig_S4_boxplots.pdf")

# ---- Statistics --------------------------------------------------------------
# Within a setting, hox_cis and A_constitutive come from the same simulations:
# paired tests. Between settings: unpaired (Welch) tests.
# Wilcoxon tests are given alongside as a non-parametric check.
within <- function(s) {
  a <- subset(d, setting == s & element == "A_constitutive")
  h <- subset(d, setting == s & element == "hox_cis")
  h <- h[match(a$sim, h$sim), ]
  data.frame(comparison = paste0("hox_cis vs A_constitutive, ", labels[s]),
             mean_diff  = mean(h$prop - a$prop, na.rm = TRUE),
             p_t        = t.test(h$prop, a$prop, paired = TRUE)$p.value,
             p_wilcox   = wilcox.test(h$prop, a$prop, paired = TRUE,
                                      exact = FALSE)$p.value)
}
between <- function(e, s1, s2) {
  x <- d$prop[d$element == e & d$setting == s1]
  y <- d$prop[d$element == e & d$setting == s2]
  data.frame(comparison = paste0(e, ", ", labels[s1], " vs ", labels[s2]),
             mean_diff  = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
             p_t        = t.test(x, y)$p.value,
             p_wilcox   = wilcox.test(x, y, exact = FALSE)$p.value)
}

stats <- rbind(
  within("with"), within("without"), within("double"),
  between("A_constitutive", "with", "without"),
  between("hox_cis",        "with", "without"),
  between("hox_trans",      "with", "without"),
  between("A_constitutive", "with", "double"),
  between("hox_cis",        "with", "double"),
  between("hox_trans",      "with", "double"))
write.csv(stats, "fig6_stats.csv", row.names = FALSE)
