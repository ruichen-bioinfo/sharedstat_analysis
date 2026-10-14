#!/usr/bin/env Rscript
## 87c_p3_composition_adjusted.R -- Figure 2A: between-/within-species response similarity for the Hsf1 core and
## other induced genes (script 03's ratio_stat), a parental quantity. Uncertainty as in script 85 for a parental
## quadratic statistic: interval widened by the larger of the rho-bar form (gene-bootstrap SE on the log scale x
## sqrt(VIF), rho-bar from parental log fold changes over species x block x time, species x time fit) and the
## leave-one-species-out jackknife (multiplier 1.96, fixed before this analysis was run);
## t(3) and mean-rho^2 forms beside.
suppressPackageStartupMessages({ library(data.table); library(limma) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
EX <- parse(file.path(root, "scripts/hsf1/03_fay_parental_dimensions.R"))
for (e in EX) { s <- paste(deparse(e), collapse = " "); if (grepl("^CMP <- list\\(\\)", s)) break
  if (grepl("^sink\\(|^fwrite\\(|quit\\(|^say\\(", s)) next; eval(e, envir = globalenv()) }
stopifnot(is.function(ratio_stat), exists("X"), exists("SPECIES"))
NB <- 1000L; T3 <- qt(0.975, 3)
rs <- function(w, idx, sp) { g <- function(a, bk) w[[paste0(a, "_", bk)]][idx]
  wc <- sapply(sp, function(a) cor(g(a, "block1"), g(a, "block2"), use = "pairwise.complete.obs"))
  bc <- sapply(combn(sp, 2, simplify = FALSE), function(p) mean(c(cor(g(p[1], "block1"), g(p[2], "block1"), use = "pairwise.complete.obs"),
                                                                  cor(g(p[1], "block2"), g(p[2], "block2"), use = "pairwise.complete.obs"))))
  median(bc) / median(wc) }
rho_par <- function(s) { W <- dcast(X[stratum == s], hog ~ species + block + time, value.var = "lfc"); M <- as.matrix(W[, -1]); M <- M[complete.cases(M), ]
  des <- model.matrix(~ 0 + factor(sub("_block[12]_", "_", colnames(M)))); set.seed(873)
  R <- M - M %*% des %*% solve(crossprod(des)) %*% t(des); g <- rownames(M) <- W$hog[complete.cases(as.matrix(W[, -1]))]
  pr <- cbind(sample(nrow(M), 8000, TRUE), sample(nrow(M), 8000, TRUE)); pr <- pr[pr[, 1] != pr[, 2], ]
  r <- vapply(seq_len(nrow(pr)), function(k) cor(R[pr[k, 1], ], R[pr[k, 2], ]), 0)
  c(rho = max(0, interGeneCorrelation(M, des)$correlation), rho2 = max(0, mean(r^2) - 1 / (ncol(M) - ncol(des)))) }
LC <- fread(file.path(root, "data/02_processed/hsf1/fay_parental_composition_37_30.tsv"))
OUT <- rbindlist(lapply(c("Hsf1_core", "other_induced"), function(s) { W <- dcast(X[stratum == s], hog + time ~ species + block, value.var = "lfc"); p <- rho_par(s)
  rbindlist(lapply(TP, function(t) { w <- W[time == t]; n <- nrow(w); r <- rs(w, seq_len(n), SPECIES); set.seed(20260905 + t)
    b <- replicate(NB, rs(w, sample.int(n, n, TRUE), SPECIES)); sb <- sd(log(b[b > 0]), na.rm = TRUE)
    jk <- sapply(SPECIES, function(a) log(rs(w, seq_len(n), setdiff(SPECIES, a)))); sej <- sqrt(3 / 4 * sum((jk - mean(jk))^2))
    se_rho <- sb * sqrt(1 + (n - 1) * p[["rho"]]); se2 <- sb * sqrt(1 + (n - 1) * p[["rho2"]]); se <- max(se_rho, sej)
    data.table(stratum = s, time = t, n = n, ratio = r, lo = exp(log(r) - 1.96 * se), hi = exp(log(r) + 1.96 * se), se_rho = se_rho, se_jk = sej,
               rho = p[["rho"]], lo_jk_t3 = exp(log(r) - T3 * sej), hi_jk_t3 = exp(log(r) + T3 * sej), lo_rho2 = exp(log(r) - 1.96 * se2),
               hi_rho2 = exp(log(r) + 1.96 * se2), lo_boot = quantile(b, .025, na.rm = TRUE), hi_boot = quantile(b, .975, na.rm = TRUE)) })) }))
ck <- merge(OUT[, .(stratum, time, ratio)], LC[, .(stratum, time, landed = ratio)], by = c("stratum", "time")); stopifnot(max(abs(ck$ratio - ck$landed)) < 1e-9)
cat(sprintf("fidelity: %d saved composition ratios reproduced\n", nrow(ck)))
fwrite(OUT, file.path(root, "data/02_processed/hsf1/p3_adjusted_composition.tsv"), sep = "\t")
print(OUT[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
