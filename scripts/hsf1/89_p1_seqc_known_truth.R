#!/usr/bin/env Rscript
## 89_p1_seqc_known_truth.R -- the SEQC known-truth test of P1, with the decision criteria fixed before this analysis was run.
## Run once per site.
## Usage: Rscript 89_p1_seqc_known_truth.R BGI   (then MAY)
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
source(file.path(root, "pkg/sharedstat/R/reliability.R"))
site <- commandArgs(TRUE)[1]; stopifnot(site %in% c("BGI", "MAY"))
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
out <- file.path(root, "data/02_processed/hsf1"); pre <- sprintf("p1_seqc_%s_", tolower(site))

## ---- data: sum 16 lane columns per library ----------------------------------------------------------------
X <- fread(cmd = sprintf("zcat %s/data/03_external/seqc/GSE47774_SEQC_ILM_%s.txt.gz", root, site))
ids <- X[[1]]; M <- as.matrix(X[, -1]); storage.mode(M) <- "double"
lab <- sub(sprintf("^SEQC_ILM_%s_([A-F])_([0-9])_.*$", site), "\\1\\2", colnames(M))
keepc <- grepl("^[A-E][1-5]$", lab)
L <- sapply(sort(unique(lab[keepc])), function(k) rowSums(M[, lab == k, drop = FALSE]))
## Deviation from the pre-specified procedure, noted before any value was computed: the layout assumed in advance says A-E have five
## libraries each; in the deposit E has two (E1, E2). E enters no computation except the duplicate-vector filter and
## thinning, so the filter uses the 22 A-E libraries that exist.
stopifnot(all(sapply(c("A", "B", "C", "D"), function(s) sum(grepl(paste0("^", s), colnames(L))) == 5L)),
          all(L == round(L)))
isE <- grepl("^ERCC", ids); nuc <- !isE
say("libraries: %d; features %d (ERCC %d)", ncol(L), nrow(L), sum(isE))

## ---- truth from the spike-ins --------------------------------------------------------------------------------
mix <- fread(file.path(root, "data/03_external/seqc/ERCC_mix_cms_095046.txt"))
setnames(mix, c("resort", "id", "subgroup", "c1", "c2", "ratio", "log2r"))
mix[, id := sub("-", "_", id)]
rho <- setNames(mix$c1 / mix$c2, mix$id); sub_of <- setNames(mix$subgroup, mix$id)
gam_ercc <- function(libs_drop = integer(0)) {
  li <- setdiff(1:5, libs_drop)
  lsz <- colSums(L[nuc, ])
  E <- sweep(L[isE, ], 2, lsz, "/")
  Cc <- paste0("C", li); Dc <- paste0("D", li)
  eC <- rowMeans(E[, Cc, drop = FALSE]); eD <- rowMeans(E[, Dc, drop = FALSE])
  mC <- rowMeans(L[isE, Cc, drop = FALSE]); mD <- rowMeans(L[isE, Dc, drop = FALSE])
  eid <- ids[isE]; r <- rho[eid]
  g <- 3 * (eD / eC) * (0.75 * r + 0.25) / (0.25 * r + 0.75)
  ok <- mC >= 100 & mD >= 100 & is.finite(g)
  data.table(id = eid, subgroup = sub_of[eid], rho = r, gamma_j = g, ok = ok)
}
GE <- gam_ercc(); gE <- GE[ok == TRUE, median(gamma_j)]
loo <- sapply(1:5, function(k) gam_ercc(k)[ok == TRUE, median(gamma_j)])
seE <- sqrt(4 / 5 * sum((loo - mean(loo))^2))
subm <- GE[ok == TRUE, .(median = median(gamma_j), n = .N), by = subgroup][order(subgroup)]
spread <- (max(subm$median) - min(subm$median)) / gE
ercc_usable <- sum(GE$ok) >= 20 && spread <= 0.25
say("gamma_ERCC = %.3f [%.3f, %.3f] over %d controls; subgroup medians %s; spread %.1f%%; usable %s",
    gE, gE - qt(.975, 4) * seE, gE + qt(.975, 4) * seE, sum(GE$ok),
    paste(sprintf("%s %.3f (n=%d)", subm$subgroup, subm$median, subm$n), collapse = ", "), 100 * spread, ercc_usable)
fwrite(GE, file.path(out, paste0(pre, "ercc.tsv")), sep = "\t")

## ---- features (fixed at full depth) --------------------------------------------------------------------------
nm <- grepl("^NM_", ids)
key <- apply(L, 1, paste, collapse = ",")
dup <- duplicated(key) & nm
cpm0 <- sweep(L, 2, colSums(L[nuc, ]) / 1e6, "/")
AD <- grep("^[ABCD]", colnames(L), value = TRUE)
feat <- nm & !dup & apply(cpm0[, AD] >= 1, 1, all)
say("features: NM_ %d, identical-count duplicates dropped %d, passing CPM >= 1 in all 20 A-D libraries %d",
    sum(nm), sum(dup), sum(feat))

## ---- one construction ------------------------------------------------------------------------------------------
rot_libs <- function(k) { s <- ((k - 1 + 0:4) %% 5) + 1; list(ref = s[1:2], rsp = s[3:4], ins = s[5]) }
build <- function(Lx, k, logscale = FALSE, t4 = FALSE) {
  lsz <- colSums(Lx[nuc, ]) / 1e6
  cpm <- sweep(Lx, 2, lsz, "/")[feat, ]
  if (logscale) V <- log2(cpm + 1) else V <- cpm / sqrt(rowMeans(cpm[, paste0("B", 1:5)]) + 1)
  r <- rot_libs(k); g <- function(s, b) V[, paste0(s, b)]
  c1 <- g("C", r$ref[1]) - g("A", r$ref[1]); c2 <- g("C", r$ref[2]) - g("A", r$ref[2])
  R <- (c1 + c2) / 2
  rsp <- if (!t4) ((g("D", r$rsp[1]) - g("A", r$rsp[1])) + (g("D", r$rsp[2]) - g("A", r$rsp[2]))) / 2 else
         ((g("D", r$rsp[1]) - g("A", r$ref[1])) + (g("D", r$rsp[2]) - g("A", r$ref[2]))) / 2
  e <- rsp - R; S <- g("C", r$ins) - g("A", r$ins)
  h <- reliability_from_blocks(c1, c2, center = TRUE)$h
  kap <- var((g("A", r$ref[1]) - g("A", r$ref[2])) / 2) / var(R)
  list(gN = 1 + cov(e, R) / var(R), gIV = 1 + cov(e, S) / cov(R, S), fs = cor(R, S), h = h, kappa = kap,
       e = e, R = R, S = S, c1 = c1, c2 = c2)
}
boot_gap <- function(b, gT, nb = 1000L) {
  set.seed(1L); n <- length(b$R)
  v <- replicate(nb, { i <- sample.int(n, n, TRUE)
    h <- reliability_from_blocks(b$c1[i], b$c2[i], center = TRUE)$h
    (1 + cov(b$e[i], b$R[i]) / var(b$R[i])) - gT * h })
  quantile(v, c(.025, .975))
}
## residual correlation between features (reported): mean off-diagonal correlation of per-library residuals
rho_bar <- function(Lx) {
  cpm <- sweep(Lx, 2, colSums(Lx[nuc, ]) / 1e6, "/")[feat, ]
  V <- cpm / sqrt(rowMeans(cpm[, paste0("B", 1:5)]) + 1)
  rr <- sapply(c("A", "C", "D"), function(s) { Z <- V[, paste0(s, 1:5)]; Z <- Z - rowMeans(Z)
    Z <- Z / sqrt(rowSums(Z^2) / 4); Z[!is.finite(Z)] <- 0; n <- nrow(Z)
    (sum(colSums(Z)^2) / 4 - n) / (n * (n - 1)) })
  mean(rr)
}

gT0 <- if (ercc_usable) gE else NA_real_
thin <- function(p, lev) { if (p == 1) return(L); set.seed(20260927 + lev)
  matrix(rbinom(length(L), size = as.integer(L), prob = p), nrow(L), dimnames = dimnames(L)) }
rows <- list(); lev <- 0L; p <- 1
repeat {
  Lx <- thin(p, lev)
  B <- lapply(1:5, function(k) build(Lx, k)); b1 <- B[[1]]
  if (lev == 0L) { gStar <- b1$gIV; gT <- if (is.na(gT0)) gStar else gT0 }
  rotv <- function(f) sapply(B, f)
  gap <- rotv(function(b) b$gN - gT * b$h); giv <- rotv(function(b) b$gIV - gT)
  tq <- qt(.975, 4)
  lb <- log(b1$gN > 0)
  T4 <- build(Lx, 1, t4 = TRUE); T5 <- build(Lx, 1, logscale = TRUE)
  bq <- boot_gap(b1, gT)
  medc <- median(apply(Lx[feat, ], 2, median))
  rows[[length(rows) + 1]] <- data.table(site = site, level = lev, p = p, median_count = medc,
    h = b1$h, gamma_naive = b1$gN, gamma_h = gT * b1$h, gap = gap[1], gap_rel = gap[1] / (gT * b1$h),
    gap_lo = gap[1] - tq * sd(gap), gap_hi = gap[1] + tq * sd(gap), gap_boot_lo = bq[1], gap_boot_hi = bq[2],
    gamma_iv = b1$gIV, first_stage = b1$fs, iv_dev = giv[1], iv_lo = giv[1] - tq * sd(giv), iv_hi = giv[1] + tq * sd(giv),
    h_rot_sd = sd(rotv(function(b) b$h)),
    t4_gamma_naive = T4$gN, t4_kappa = T4$kappa, t4_pred = gT * T4$h + T4$kappa,
    t5_h = T5$h, t5_gamma_naive = T5$gN, t5_gamma_iv = T5$gIV,
    rho_bar = if (lev %% 3 == 0) rho_bar(Lx) else NA_real_)
  r <- rows[[length(rows)]]
  say("level %2d p=2^-%-2d  median count %9.1f  h %.3f  naive %.3f  gamma*h %.3f  gap %+.3f [%+.3f, %+.3f]  IV %.3f (fs %.2f)",
      lev, lev, medc, r$h, r$gamma_naive, r$gamma_h, r$gap, r$gap_lo, r$gap_hi, r$gamma_iv, r$first_stage)
  if (b1$h < 1 / gT - 0.10 || medc < 2 || p < 2^-20) break
  lev <- lev + 1L; p <- p / 2
}
RES <- rbindlist(rows)
fwrite(RES, file.path(out, paste0(pre, "levels.tsv")), sep = "\t")

## ---- verdicts --------------------------------------------------------------------------------------------------
say("\ntruth: %s = %.3f (gamma_ERCC %.3f, usable %s; gamma* %.3f)", if (is.na(gT0)) "gamma*" else "gamma_ERCC",
    gT, gE, ercc_usable, gStar)
T0 <- gT > 1 && gT < 9 && RES$first_stage[1] >= 0.5 && abs(gStar - gT) / gT <= 0.25
J <- RES[h >= 0.10]
t1lev <- J[, abs(gap_rel) <= 0.10 & gap_lo <= 0 & gap_hi >= 0]
T1 <- length(t1lev) > 0 && all(t1lev)
low <- RES[h < 1 / gT - 0.05]; high <- J[h > 1 / gT + 0.05]
T2 <- if (!nrow(low)) NA else (nrow(high) > 0 && all(low$gamma_naive < 1) && all(high$gamma_naive > 1))
J3 <- J[first_stage >= 0.3]
t3lev <- J3[, abs(iv_dev) / gT <= 0.10 & iv_lo <= 0 & iv_hi >= 0]
T3 <- length(t3lev) > 0 && all(t3lev)
V <- data.table(site = site, truth = if (is.na(gT0)) "gamma_star" else "gamma_ERCC", gamma_truth = gT,
                gamma_ercc = gE, gamma_ercc_lo = gE - qt(.975, 4) * seE, gamma_ercc_hi = gE + qt(.975, 4) * seE,
                ercc_n = sum(GE$ok), ercc_spread = spread, gamma_star = gStar,
                n_features = sum(feat), n_levels = nrow(RES), n_judged = nrow(J),
                T0 = T0, T1 = T1, T1_levels_pass = sum(t1lev), T2 = T2, T3 = T3, T3_levels_pass = sum(t3lev),
                T3_levels_judged = length(t3lev))
fwrite(V, file.path(out, paste0(pre, "verdicts.tsv")), sep = "\t")
say("T0 %s | T1 %s (%d of %d judged levels) | T2 %s | T3 %s (%d of %d)", T0, T1, sum(t1lev), length(t1lev),
    if (is.na(T2)) "not reached" else T2, T3, sum(t3lev), length(t3lev))
