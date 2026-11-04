# Figure 5 (and SFigs S1, S3, Table S2): fitness-adjusted differences in
# parameter ranges between settings.
#
# Input (in data_dir), one row per simulation, one column per parameter plus
# 'fit' (final fitness) and 'type' (setting x element x tooth, see below):
#   nohoxsiniend.dat, hoxsiniend.dat   range between first and last generation
#   nohoxsminmax.dat, hoxsminmax.dat   range between min and max over evolution
#   'nohox' files hold the A_constitutive parameters and hn (hox_trans),
#   'hox' files the hox_cis parameters.
# type codes: 1/2 = A_const upper/lower, 3/4 = hox_cis upper/lower in the
#   with-pleiotropy setting; +4 without pleiotropy; +8 double mutation.
# Output: fig5A_example.pdf, fig5B_violins.pdf, sfig_S1_bars.pdf,
#   sfig_S3_violins.pdf, table_S2.csv

library(multcomp)   # Tukey contrasts
library(vioplot)

source(here::here("analysis", "config.R")); data_dir <- data_proc
settings <- c("with", "without", "double")
labels   <- c(with = "With pleiotropy", without = "Without pleiotropy",
              double = "Without pleiotropy, double mutation")
cols     <- adjustcolor(c("red", "blue", "purple"), 0.5)
names(cols) <- settings
exclude  <- "INH_Inhibtion_of_SHH_over_BMP4"   # parameter not analyzed

normalize <- function(x) (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))

# ---- Read and reshape to long format ---------------------------------------
read_ranges <- function(file) {
  x <- read.table(file.path(data_dir, file))
  params <- setdiff(names(x), c("fit", "type", "generation", "row_num"))
  params <- params[!grepl(exclude, params)]
  do.call(rbind, lapply(params, function(p) data.frame(
    parameter = sub("(_hox)?\\.\\d+$", "", p),
    value     = x[[p]],
    fitness   = x$fit,
    type      = x$type)))
}

decode <- function(d) {
  t <- (d$type - 1) %% 4                       # 0-3 within a setting
  d$setting <- factor(settings[(d$type - 1) %/% 4 + 1], levels = settings)
  d$element <- ifelse(t < 2, "A_constitutive", "hox_cis")
  d$element[d$parameter == "hn"] <- "hox_trans"
  d$tooth   <- ifelse(t %% 2 == 0, "upper", "lower")
  d
}

load_range <- function(range_type) {
  a <- read_ranges(paste0("nohoxs", range_type, ".dat"))
  h <- read_ranges(paste0("hoxs",   range_type, ".dat"))
  d <- decode(rbind(a, h))
  # min-max normalization per parameter, jointly across settings and elements
  # (hn, which exists only as hox_trans, is normalized on its own)
  for (p in unique(d$parameter)) {
    i <- d$parameter == p
    d$value[i] <- normalize(d$value[i])
  }
  d$range_type <- range_type
  d
}

d <- rbind(load_range("iniend"), load_range("minmax"))

# ---- Fitness-adjusted comparison per parameter, element, tooth -------------
# range ~ fitness + setting (common fitness slope). The setting coefficients
# are the differences in intercepts, i.e. the difference in range at equal
# fitness; significance from Tukey's all-pairwise comparisons.
groups <- unique(d[, c("range_type", "element", "tooth", "parameter")])
table_S2 <- do.call(rbind, lapply(seq_len(nrow(groups)), function(g) {
  sub <- merge(d, groups[g, ])
  m   <- lm(value ~ fitness + setting, data = sub)
  tk  <- summary(glht(m, linfct = mcp(setting = "Tukey")))
  cf  <- coef(m)
  cbind(groups[g, ],
        delta_with_without = -cf["settingwithout"],
        delta_with_double  = -cf["settingdouble"],
        p_with_without     = tk$test$pvalues[1],
        p_with_double      = tk$test$pvalues[2],
        p_without_double   = tk$test$pvalues[3])
}))
rownames(table_S2) <- NULL
write.csv(table_S2, "table_S2.csv", row.names = FALSE)

# ---- Fig 5A: example parameter ---------------------------------------------
example <- subset(d, range_type == "iniend" & element == "A_constitutive" &
                    tooth == "upper" &
                    parameter == "Egr_Epithelial_proliferation_rate")
m  <- lm(value ~ fitness + setting, data = example)
fx <- seq(min(example$fitness), max(example$fitness), length.out = 100)

pdf("fig5A_example.pdf", width = 6, height = 5)
plot(value ~ fitness, data = example, col = cols[example$setting], pch = 16,
     xlab = "Fitness", ylab = "Parameter range (normalized 0 to 1)",
     main = "Epithelial proliferation, A_constitutive")
for (s in settings) {
  pr <- predict(m, data.frame(fitness = fx, setting = s), interval = "confidence")
  polygon(c(fx, rev(fx)), c(pr[, "lwr"], rev(pr[, "upr"])), col = cols[s], border = NA)
  lines(fx, pr[, "fit"], col = cols[s], lwd = 2)
  points(min(fx), pr[1, "fit"], pch = 4, cex = 2, lwd = 2)   # intercept (cross)
}
legend("topleft", legend = labels, col = cols, pch = 16, bty = "n")
dev.off()

# ---- Fig 5B and SFig S3: violins of intercept differences ------------------
# White: significantly larger range with pleiotropy; black: significantly larger
# without; grey: no significant difference. Upper and lower teeth pooled.
violins <- function(delta, p, file) {
  pdf(file, width = 8, height = 5)
  panels <- expand.grid(element = c("A_constitutive", "hox_cis"),
                        range_type = c("iniend", "minmax"),
                        stringsAsFactors = FALSE)
  par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
  for (k in seq_len(nrow(panels))) {
    s   <- merge(table_S2, panels[k, ])
    y   <- s[[delta]]
    col <- ifelse(s[[p]] >= 0.05, "grey", ifelse(y > 0, "white", "black"))
    vioplot(y, drawRect = FALSE, col = "white", ylim = c(-0.35, 0.35),
            names = panels$element[k],
            main = ifelse(panels$range_type[k] == "iniend",
                          "First to last generation", "Min. to max. value"))
    abline(h = 0, lty = 2)
    points(jitter(rep(1, length(y)), 5), y, pch = 21, bg = col, cex = 1.5)
  }
  dev.off()
}
violins("delta_with_without", "p_with_without", "fig5B_violins.pdf")
violins("delta_with_double",  "p_with_double",  "sfig_S3_violins.pdf")

# ---- SFig S1: intercept difference per parameter ---------------------------
pdf("sfig_S1_bars.pdf", width = 12, height = 8)
par(mfrow = c(2, 1), mar = c(12, 4, 3, 1))
for (r in c("iniend", "minmax")) {
  s <- table_S2[table_S2$range_type == r & table_S2$element != "hox_trans", ]
  s <- s[order(s$element, s$parameter, s$tooth), ]
  barplot(s$delta_with_without, las = 2, cex.names = 0.5,
          names.arg = paste(sub("_.*", "", s$parameter), s$element),
          col = ifelse(s$p_with_without >= 0.05, "grey",
                       ifelse(s$tooth == "upper", "darkgreen", "orange")),
          ylab = "Intercept difference (with - without)",
          main = ifelse(r == "iniend", "First to last generation",
                        "Min. to max. value"))
}
dev.off()
