# Figure 4: total acceleration of fitness versus final fitness, per setting,
# plus the corresponding statistics.
#
# Input (in data_dir):
#   infocategories_threshol0.1_50pbest_delta1_method2.dat
#     one row per simulation, with columns experiment, evo (initial condition),
#     mean_fitness (inverse fitness) and sum.abs.differences.. (total
#     acceleration: sum of the absolute second derivative of fitness)
# Output (in the working directory):
#   fig4_acceleration.pdf, fig4_stats_*.csv

library(car)    # leveneTest
library(coin)   # wilcox_test, kruskal_test



source(here::here("analysis", "config.R")); data_dir <- data_proc
settings <- c("1genome", "2genomes", "2genomesHigh")
labels   <- c("With pleiotropy", "Without pleiotropy",
              "Without pleiotropy, double mutation")
cols     <- adjustcolor(c("red", "blue", "purple"), alpha.f = 0.5)
shapes   <- c(16, 17, 8)
names(labels) <- names(cols) <- names(shapes) <- settings

# ---- Data ------------------------------------------------------------------
d <- read.table(file.path(data_dir,
                          "infocategories_threshol0.1_50pbest_delta1_method2.dat"))
d <- d[d$experiment %in% settings, ]
d$exp          <- factor(d$experiment, levels = settings)
d$fitness      <- 1 / d$mean_fitness
d$acceleration <- d$sum.abs.differences..
if ("evo" %in% names(d)) d$ic <- factor(d$evo)   # initial condition

pairs <- combn(settings, 2, simplify = FALSE)

# ---- Fig 4: acceleration vs final fitness ----------------------------------
pdf("fig4_acceleration.pdf", width = 6, height = 6)
plot(acceleration ~ fitness, data = d, type = "n",
     xlim = c(0, 52), ylim = c(0, 20),
     xlab = "Fitness", ylab = "Total acceleration of fitness")
for (s in settings)
  points(acceleration ~ fitness, data = d[d$exp == s, ],
         col = cols[s], pch = shapes[s], cex = 1.2)
legend("topleft", legend = labels, col = cols, pch = shapes, bty = "n")
dev.off()

# ---- Statistics --------------------------------------------------------------
# Normality per setting (Shapiro-Wilk) and homoscedasticity (Levene)
shapiro <- do.call(rbind, lapply(settings, function(s) {
  t <- shapiro.test(d$acceleration[d$exp == s])
  data.frame(setting = labels[s], W = t$statistic, p = t$p.value)
}))
lev    <- leveneTest(acceleration ~ exp, data = d)
levene <- data.frame(F = lev[["F value"]][1], p = lev[["Pr(>F)"]][1])

# Setting effect at equal fitness (reported in Fig 4 and Table S4):
#   acceleration ~ fitness   vs   acceleration ~ fitness + setting
# and test of equal slopes (fitness x setting interaction)
m_fit <- lm(acceleration ~ fitness,       data = d)
m_add <- lm(acceleration ~ fitness + exp, data = d)
m_int <- lm(acceleration ~ fitness * exp, data = d)
cmp <- function(a, b, label) {
  av <- anova(a, b)
  data.frame(test = label, df = av$Df[2], F = av$F[2], p = av[["Pr(>F)"]][2])
}
ancova <- rbind(cmp(m_fit, m_add, "setting (fitness as covariate)"),
                cmp(m_add, m_int, "fitness x setting interaction"))
coefs  <- as.data.frame(summary(m_add)$coefficients)   # setting differences

# Differences between settings ignoring fitness: Kruskal-Wallis + Wilcoxon
kw <- kruskal_test(acceleration ~ exp, data = d, distribution = "asymptotic")
pw <- sapply(pairs, function(p)
  pvalue(wilcox_test(acceleration ~ exp, data = droplevels(d[d$exp %in% p, ]),
                     distribution = "asymptotic")))
between <- data.frame(
  test  = c("Kruskal-Wallis",
            sapply(pairs, function(p) paste(labels[p], collapse = " vs "))),
  p     = c(pvalue(kw), pw),
  p_adj = c(NA, p.adjust(pw, "BH")))

# Effect of initial conditions on acceleration within each setting
if (!is.null(d$ic)) {
  ic_effect <- do.call(rbind, lapply(settings, function(s) {
    t <- kruskal_test(acceleration ~ ic, data = droplevels(d[d$exp == s, ]),
                      distribution = "asymptotic")
    data.frame(setting = labels[s], chi2 = statistic(t), p = pvalue(t))
  }))
  write.csv(ic_effect, "fig4_stats_initial_conditions.csv", row.names = FALSE)
}

write.csv(shapiro, "fig4_stats_shapiro.csv", row.names = FALSE)
write.csv(levene,  "fig4_stats_levene.csv",  row.names = FALSE)
write.csv(ancova,  "fig4_stats_ancova.csv",  row.names = FALSE)
write.csv(coefs,   "fig4_stats_ancova_coefficients.csv")
write.csv(between, "fig4_stats_between_settings.csv", row.names = FALSE)
