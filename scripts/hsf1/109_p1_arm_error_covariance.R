#!/usr/bin/env Rscript
## 109_p1_arm_error_covariance.R -- P1 arm error covariance, with the estimand fixed before this analysis was run.
## Measures kappa = cov(u, v)/var(R-hat) from the two replicate blocks, on script 91's log-CPM block responses.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>"); set.seed(20261005L)
ext <- file.path(root, "data/03_external/fay2023"); hsd <- file.path(root, "data/02_processed/hsf1")
HYB <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu"); TP <- c(15L, 30L, 60L)
dat <- read.csv(file.path(ext, "htseq.ortho.csv")); design <- read.csv(file.path(ext, "design.csv"))
one <- dat$one == 1; cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
D <- as.data.table(design); D[, col := Name]; D[, num := as.integer(sub("^[A-Za-z]+", "", Name))]
D[, block := fifelse(num <= 32L, "block1", "block2")]
H <- D[Strain %in% HYB & Condition == "heat"]; H[, Time := as.integer(Time)]
L <- list()
for (g in HYB) {
  hs <- H[Strain == g]; a1 <- substr(g, 1, 2); a2 <- substr(g, 4, 5)
  m <- cts[, hs$col, drop = FALSE]; y <- calcNormFactors(DGEList(counts = m), method = "TMM")
  lc <- cpm(y, log = TRUE, prior.count = 2); cr <- cpm(y)
  keep <- rowSums(cr[, hs[Allele == a1, col]] >= 1) >= 4L & rowSums(cr[, hs[Allele == a2, col]] >= 1) >= 4L
  gv <- function(al, bk, t) lc[keep, hs[Allele == al & block == bk & Time == t, col]]
  for (t in TP) L[[length(L) + 1L]] <- data.table(hog = rownames(lc)[keep], hybrid = g, time = t,
    rA1 = gv(a1, "block1", t) - gv(a1, "block1", 0L), rA2 = gv(a1, "block2", t) - gv(a1, "block2", 0L),
    rB1 = gv(a2, "block1", t) - gv(a2, "block1", 0L), rB2 = gv(a2, "block2", t) - gv(a2, "block2", 0L))
}
X <- rbindlist(L)
X[, `:=`(R = (rA1 + rA2) / 2, A = (rB1 + rB2) / 2, nR = (rA2 - rA1) / 2, nA = (rB2 - rB1) / 2)]
## centre within hybrid x time so that stacking timepoints does not add between-time variance
X[, `:=`(R = R - mean(R), A = A - mean(A), nR = nR - mean(nR), nA = nA - mean(nA)), by = .(hybrid, time)]
est <- function(d) { vR <- mean(d$R^2); h <- 1 - mean(d$nR^2) / vR; k <- mean(d$nR * d$nA) / vR
  nv <- sum(d$A * d$R) / sum(d$R^2)
  c(n = nrow(d), naive = nv, h = h, kappa = k, gamma_h = nv / h, gamma_hk = (nv - k) / h,
    rho_err = mean(d$nR * d$nA) / sqrt(mean(d$nR^2) * mean(d$nA^2)),
    ## also fixed in advance: reverse orientation and amplitude ratio
    gamma_rev = { vA <- mean(d$A^2); hA <- 1 - mean(d$nA^2) / vA; kA <- mean(d$nR * d$nA) / vA
                  (sum(d$A * d$R) / sum(d$A^2) - kA) / hA },
    amp = sqrt((mean(d$A^2) - mean(d$nA^2)) / (mean(d$R^2) - mean(d$nR^2))),
    r_true = (mean(d$A * d$R) - mean(d$nR * d$nA)) /
             sqrt((mean(d$R^2) - mean(d$nR^2)) * (mean(d$A^2) - mean(d$nA^2)))) }
B <- 1000L; out <- list()
for (g in HYB) for (tt in list(TP, 15L, 30L, 60L)) {
  d <- X[hybrid == g & time %in% tt]; e0 <- est(d); genes <- unique(d$hog); setkey(d, hog)
  bs <- t(replicate(B, est(d[.(sample(genes, replace = TRUE)), allow.cartesian = TRUE])))
  ci <- apply(bs, 2, quantile, c(0.025, 0.975))
  out[[length(out) + 1L]] <- data.table(hybrid = g, time = if (length(tt) > 1) "all" else as.character(tt),
    stat = names(e0), est = e0, lo = ci[1, ], hi = ci[2, ])
}
O <- rbindlist(out)[stat != "n" | TRUE]
fwrite(O, file.path(hsd, "p1_arm_error_covariance.tsv"), sep = "\t")
W <- dcast(O[time == "all"], hybrid ~ stat, value.var = "est")
print(W[, .(hybrid, n, naive = round(naive, 3), h = round(h, 3), kappa = round(kappa, 4), rho_err = round(rho_err, 3),
            gamma_h = round(gamma_h, 3), gamma_hk = round(gamma_hk, 3))])
print(O[time == "all" & stat %in% c("kappa", "gamma_hk", "gamma_rev", "amp", "r_true"), .(hybrid, stat, est = round(est, 4), lo = round(lo, 4), hi = round(hi, 4))])
print(dcast(O[stat == "kappa"], hybrid ~ time, value.var = "est")[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)])
