# 110_fig_p1_hybrid_gains.R -- P1 Figure 5, the six hybrids' gains under each reading, from saved tables only:
# naive, naive/h and instrumented with 100-cluster joint bootstrap intervals (p1_table1_joint_bootstrap.tsv, k = 100),
# and the gain adjusted for h and measured kappa with gene-bootstrap intervals (p1_arm_error_covariance.tsv, script 109).
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); fig <- file.path(root, "doc/04_figures/final")
source(file.path(root, "scripts/figures/theme.R")); source(file.path(root, "scripts/figures/fig_style.R"))
T1 <- fread(file.path(root, "data/02_processed/crossspecies/p1_table1_joint_bootstrap.tsv"))[k == "100"]
K <- dcast(fread(file.path(root, "data/02_processed/hsf1/p1_arm_error_covariance.tsv"))[time == "all" & stat == "gamma_hk"],
           hybrid ~ stat, value.var = c("est", "lo", "hi"))
setnames(K, c("group", "hk", "hk_lo", "hk_hi")); D <- merge(T1, K, by = "group")
stopifnot(nrow(D) == 6L)
D <- D[order(-naive)]
LAB <- c(ScxSk = "S. cerevisiae × S. kudriavzevii", ScxSp = "S. cerevisiae × S. paradoxus", ScxSu = "S. cerevisiae × S. uvarum",
         SkxSu = "S. kudriavzevii × S. uvarum", SpxSk = "S. paradoxus × S. kudriavzevii", SpxSu = "S. paradoxus × S. uvarum")
draw <- function() {
  par(mar = c(3.6, 13.2, 5.6, 0.8)); y <- seq_len(nrow(D)); off <- c(0.27, 0.09, -0.09, -0.27)
  plot(NA, xlim = c(0.6, 1.25), ylim = c(0.5, nrow(D) + 0.5), yaxt = "n", xlab = "gain (1 = no divergence)", ylab = "")
  axis(2, at = y, labels = FALSE); text(par("usr")[1] - 0.012, y, LAB[D$group], adj = 1, xpd = NA, cex = CEX_TXT, font = 3)
  abline(v = 1, lty = 2, col = CHARCOAL, lwd = LW$ref)
  segments(D$rel_lo, y + off[2], D$rel_hi, y + off[2], col = PAL$rel, lwd = LW$ci)
  segments(D$iv_lo, y + off[3], D$iv_hi, y + off[3], col = PAL$corrected, lwd = LW$ci)
  segments(D$hk_lo, y + off[4], D$hk_hi, y + off[4], col = PAL$hk, lwd = LW$ci)
  points(D$naive, y + off[1], pch = 22, bg = COL_SCI$naive, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  points(D$rel, y + off[2], pch = 21, bg = PAL$rel, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  points(D$iv, y + off[3], pch = 21, bg = COL_SCI$corrected, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  points(D$hk, y + off[4], pch = 23, bg = PAL$hk, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  legend_top(c("naive gain", "naive/h", "instrumented gain", "adjusted for h and measured κ"),
             pch = c(22, 21, 21, 23), pt.bg = c(COL_SCI$naive, PAL$rel, COL_SCI$corrected, PAL$hk),
             col = COL_SCI$truth, pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2)
}
fig_both("FP1_5_hybrid_gains", 108, 86, draw, fig)
cat("FP1_5_hybrid_gains ok\n")
