#!/usr/bin/env Rscript
## 84_p2_rho2_sensitivity.R -- sensitivity of P2's presence and 37/42 ratio verdicts (script 80) to the
## correlation model. Script 80 inflated these amplitude-type (quadratic) statistics with VIF = 1 + (n - 1) rho;
## for a quadratic statistic the governing quantity is the bias-corrected mean squared pairwise correlation
## (the mean-rho^2 form of script 85). P2's verdicts were fixed before this analysis was run,
## so this is reported as a sensitivity only.
## Same inputs as 80.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
CM <- fread(file.path(root, "data/02_processed/count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(root, "data/01_metadata/sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + factor(SM$grp)); v <- voom(dge, des, plot = FALSE)
R <- v$E - v$E %*% des %*% solve(crossprod(des)) %*% t(des); df <- ncol(R) - ncol(des)
set.seed(20260926L); G <- nrow(R); pr <- cbind(sample.int(G, 40000, TRUE), sample.int(G, 40000, TRUE)); pr <- pr[pr[, 1] != pr[, 2], ]
r <- vapply(seq_len(nrow(pr)), function(k) cor(R[pr[k, 1], ], R[pr[k, 2], ]), 0)
rho <- mean(r); rho2 <- max(0, mean(r^2) - 1 / df)
V1 <- 1 + (G - 1) * rho; V2 <- 1 + (G - 1) * rho2
BP <- fread(file.path(hsd, "p2_detectability_bypairing.tsv"))
for (nm in c("1", "2")) { V <- if (nm == "1") V1 else V2
  BP[, (paste0("lo_", nm)) := pmax(0, A - sqrt(V) * (A - lo))][, (paste0("hi_", nm)) := A + sqrt(V) * (hi - A)] }
Q1 <- BP[, .(zero_pairings = sum(A == 0), present_rho = all(A > 0) & all(lo_1 > 0), present_rho2 = all(A > 0) & all(lo_2 > 0),
             min_lo_rho = min(lo_1), min_lo_rho2 = min(lo_2)), by = .(scale, condition)]
Q2 <- rbindlist(lapply(c("level", "response"), function(sc) rbindlist(lapply(c("10", "30"), function(tp) {
  a <- BP[scale == sc & condition == paste0("37_", tp)]; b <- BP[scale == sc & condition == paste0("42_", tp)]
  m <- merge(a[, .(labels_25, pairing, l1 = lo_1, l2 = lo_2)], b[, .(labels_25, pairing, h1 = hi_1, h2 = hi_2)], by = c("labels_25", "pairing"))
  f <- function(l, h) ifelse(h > 0, l / h, ifelse(l > 0, Inf, 0))
  data.table(scale = sc, timepoint = tp, holds_rho = all(f(m$l1, m$h1) > 1), holds_rho2 = all(f(m$l2, m$h2) > 1),
             min_lo_rho = min(f(m$l1, m$h1)), min_lo_rho2 = min(f(m$l2, m$h2))) }))))
## check: validate the rho-form results against script 80
LQ <- fread(file.path(hsd, "p2_presence_adjusted.tsv")); LR <- fread(file.path(hsd, "p2_ratio_adjusted.tsv"))
c1 <- merge(Q1, LQ[, .(scale, condition, present)], by = c("scale", "condition")); c2 <- merge(Q2, LR[, .(scale, timepoint = as.character(timepoint), holds)], by = c("scale", "timepoint"))
cat(sprintf("genome-wide residual rho %.4f (VIF %.1f); bias-corrected mean rho^2 %.5f (VIF %.1f); residual df %d\n", rho, V1, rho2, V2, df))
cat(sprintf("rho-form reproduces the saved results: presence %s, ratio %s\n", all(c1$present_rho == c1$present), all(c2$holds_rho == c2$holds)))
fwrite(Q1, file.path(hsd, "p2_presence_rho2_sensitivity.tsv"), sep = "\t"); fwrite(Q2, file.path(hsd, "p2_ratio_rho2_sensitivity.tsv"), sep = "\t")
print(Q1); print(Q2)
