#!/usr/bin/env Rscript
## 82_p1_h_correlation.R -- inter-gene correlation in P1's allele-paired design, and the reliability intervals and
## reversal-regime counts with correlation-adjusted width. Reads stored tables and the deposited counts only.
## (1) Per hybrid, the 16 allele-level libraries (2 alleles x 4 times x 2 blocks) are normalised as script 04 does
##     (per-hybrid TMM, logCPM, prior.count from 04), a group-means model on allele x time (8 groups, 8 residual df) is
##     fitted, and for each P1 stratum the mean residual correlation rho is taken with limma::interGeneCorrelation.
##     Also reported: rho after centring each library's residuals across genes (the part a centred estimator sees) and
##     the bias-corrected mean squared correlation (mean r^2 - 1/df over random gene pairs), the quantity that governs
##     the variance of a mean of squares.
## (2) The centred reliability intervals of P1 (p1x_b5_reliability_centred_ci.tsv, the Figure 3B source) are widened
##     about the point estimate by sqrt(VIF), VIF = 1 + (n - 1) rho (the convention of script 80 for P2, conservative
##     for a variance ratio), and by sqrt(VIF2), VIF2 = 1 + (n - 1) mean_rho2 (the variance-of-squares form). The
##     reversal-regime counts (whole interval below 1/gamma) are recounted under both; point-estimate counts unchanged.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
ext <- file.path(root, "data/03_external/fay2023")
say <- function(...) cat(sprintf(...), "\n", sep = "")
src <- readLines(file.path(root, "scripts/hsf1/04_divergence_by_dependency.R"))
PRIOR <- eval(parse(text = sub(".*<-", "", grep("^PRIOR *<-", src, value = TRUE)[1])))
stopifnot(is.numeric(PRIOR))
dat <- read.csv(file.path(ext, "htseq.ortho.csv")); design <- read.csv(file.path(ext, "design.csv"))
one <- dat$one == 1; cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
D <- as.data.table(design); D[, col := Name][, num := as.integer(sub("^[A-Za-z]+", "", Name))]
D[, block := fifelse(num <= 32L, "block1", "block2")]
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu")
H <- D[Strain %in% HYB & Condition == "heat"][, Time := as.integer(Time)]
X <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a5_cis_dod_37_30.tsv.gz"))))
set.seed(20260926L)
## Note: the reliability is of the allele contrast inside one library, and library-level noise shared by both
## alleles cancels in it. The relevant correlation is therefore that of the
## per-library allele difference (a2 - a1 in the same library), 8 libraries per hybrid (4 times x 2 blocks), with
## residuals from a time-means fit (4 groups, 4 residual df). The allele-level version is kept as a column, because it
## is what a reader would compute first and it overstates the correlation.
RH <- rbindlist(lapply(HYB, function(hy) {
  hs <- H[Strain == hy]; m <- cts[, hs$col]
  y <- calcNormFactors(DGEList(counts = m), method = "TMM"); lc <- cpm(y, log = TRUE, prior.count = PRIOR)
  a1 <- substr(hy, 1, 2); a2 <- substr(hy, 4, 5)
  L1 <- hs[Allele == a1][order(num)]; L2 <- hs[Allele == a2][order(num)]
  stopifnot(identical(L1$num, L2$num), identical(L1$Time, L2$Time), nrow(L1) == 8L)
  Dm <- lc[, L2$col] - lc[, L1$col]; colnames(Dm) <- paste0("lib", L1$num)
  des <- model.matrix(~ 0 + factor(L1$Time))
  R <- Dm - Dm %*% des %*% solve(crossprod(des)) %*% t(des)          # residuals, 4 df
  grpA <- factor(paste(hs$Allele, hs$Time)); desA <- model.matrix(~ 0 + grpA)
  rbindlist(lapply(unique(X$stratum), function(st) {
    g <- intersect(unique(X[hybrid == hy & stratum == st, hog]), rownames(lc)); if (length(g) < 10) return(NULL)
    gs <- if (length(g) > 2000) sample(g, 2000) else g
    rho <- interGeneCorrelation(Dm[gs, ], des)$correlation
    rho_allele <- interGeneCorrelation(lc[gs, ], desA)$correlation
    pr <- matrix(sample(gs, 2 * min(20000, choose(length(gs), 2)), replace = TRUE), ncol = 2); pr <- pr[pr[, 1] != pr[, 2], ]
    r <- vapply(seq_len(nrow(pr)), function(k) cor(R[pr[k, 1], ], R[pr[k, 2], ]), 0)
    data.table(hybrid = hy, stratum = st, n_genes = length(g), rho = rho, rho_allele_level = rho_allele,
               mean_r2 = mean(r^2, na.rm = TRUE), rho2_bc = max(0, mean(r^2, na.rm = TRUE) - 1 / 4))
  }))
}))
fwrite(RH, file.path(hsd, "p1_intergene_correlation.tsv"), sep = "\t")
B5 <- fread(file.path(hsd, "p1x_b5_reliability_centred_ci.tsv"))[dataset == "Fay 2023"]
B5[, c("stratum", "hybrid", "time") := tstrsplit(arm, "|", fixed = TRUE)]
B5 <- merge(B5, RH[, .(stratum, hybrid, rho, rho2_bc)], by = c("stratum", "hybrid"))
stopifnot(nrow(B5) == 108L)
B5[, `:=`(vif = 1 + (n_genes - 1) * pmax(rho, 0), vif2 = 1 + (n_genes - 1) * rho2_bc)]
B5[, `:=`(hi_adj = pmin(1, h_centred + sqrt(vif) * (hi - h_centred)), lo_adj = pmax(0, h_centred - sqrt(vif) * (h_centred - lo)),
          hi_adj2 = pmin(1, h_centred + sqrt(vif2) * (hi - h_centred)), lo_adj2 = pmax(0, h_centred - sqrt(vif2) * (h_centred - lo)))]
fwrite(B5, file.path(hsd, "p1_h_adjusted.tsv"), sep = "\t")
G <- c(1.1, 1.2, 1.3, 1.5, 2.0)
CNT <- rbindlist(lapply(G, function(g) data.table(gamma = g, point = sum(B5$h_centred < 1/g), interval_gene = sum(B5$hi < 1/g),
                                                  interval_vif_rho = sum(B5$hi_adj < 1/g), interval_vif_rho2 = sum(B5$hi_adj2 < 1/g))))
fwrite(CNT, file.path(hsd, "p1_reversal_counts_adjusted.tsv"), sep = "\t")
say("=== 82_p1_h_correlation | %s ===", format(Sys.time()))
say("residual correlation by stratum (median over hybrids):")
print(RH[, .(rho = round(median(rho), 4), rho_allele_level = round(median(rho_allele_level), 4), mean_r2 = round(median(mean_r2), 4),
             rho2_bc = round(median(rho2_bc), 4), n = as.integer(median(n_genes))), by = stratum])
say("\nVIF over the 108 cells: rho form median %.1f (range %.1f-%.1f); rho^2 form median %.2f (range %.2f-%.2f)",
    median(B5$vif), min(B5$vif), max(B5$vif), median(B5$vif2), min(B5$vif2), max(B5$vif2))
say("\nreversal-regime counts (of 108):"); print(CNT)
