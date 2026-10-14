#!/usr/bin/env Rscript
## 107_p2_per_library_score.R -- per-library class score for P2. For the Msn2/4-only, redundant and Hsf1-only classes,
## each library's median log-CPM over the class minus that over the not-induced reference; I(37/30), I(37/30) - I(42/30)
## and (H-W) 37/30 - (H-W) 37/10 are estimated from the group means with a pooled-variance t interval. The estimands
## were fixed before this analysis was run.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT"); hsd <- file.path(root, "data/02_processed/hsf1"); pro <- file.path(root, "data/02_processed")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(root, "data/01_metadata/sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); dge <- calcNormFactors(dge[filterByExpr(dge, group = SM$grp), , keep.lib.sizes = FALSE])
E <- cpm(dge, log = TRUE, prior.count = 0.5)
CLS <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
U <- CLS[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][k == 1L & nm == 3L]
ref <- intersect(rownames(E), U[cls == "not_induced", gene])
G4 <- c("Double.KDKO", "HSF1.KD", "MSN24.KO", "Wildtype")
OUT <- rbindlist(lapply(c("Msn24_only", "redundant", "Hsf1_only"), function(cl) {
  g <- intersect(rownames(E), U[cls == cl, gene])
  sc <- apply(E[g, ], 2, median) - apply(E[ref, ], 2, median)
  D <- data.table(SM[, .(knockout, cond, replicate)], score = sc)
  mu <- D[, .(m = mean(score)), by = .(knockout, cond)]; s2 <- D[, var(score) * (.N - 1), by = .(knockout, cond)][, sum(V1)] / (nrow(D) - nrow(mu))
  df <- nrow(D) - nrow(mu); m <- function(k, c) mu[knockout == k & cond == c, m]
  est <- function(w) { v <- sum(w$coef * mapply(m, w$k, w$c)); se <- sqrt(s2 * sum(w$coef^2) / 3); c(v, v - qt(.975, df) * se, v + qt(.975, df) * se) }
  I <- function(c) list(k = G4, c = rep(c, 4), coef = c(1, -1, -1, 1))
  cat2 <- function(a, b) list(k = c(a$k, b$k), c = c(a$c, b$c), coef = c(a$coef, -b$coef))
  HW <- function(c) list(k = c("HSF1.KD", "Wildtype"), c = rep(c, 2), coef = c(1, -1))
  r <- rbind(`I(37/30)` = est(I("37_30")), `I(37/30) - I(42/30)` = est(cat2(I("37_30"), I("42_30"))),
             `(H-W) 37/30 - (H-W) 37/10` = est(cat2(HW("37_30"), HW("37_10"))))
  data.table(cls = cl, n_genes = length(g), estimand = rownames(r), est = r[, 1], lo = r[, 2], hi = r[, 3], df = df) }))
fwrite(OUT, file.path(hsd, "p2_per_library_score.tsv"), sep = "\t")
print(OUT[, .(cls, n_genes, estimand, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3), df)])
