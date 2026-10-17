#!/usr/bin/env Rscript
## 102d_p1_joint_cluster_bootstrap.R -- joint cluster bootstrap for P1 Table 1. Genes are clustered by k-means on their
## standardised limma R across the six hybrids (k = 100, 50, 200; k = 0 is a plain gene bootstrap), clusters are
## resampled (2000 draws), and naive, IV, h, naive/h and the closure naive - IV x h are recomputed jointly per hybrid.
## The procedure was fixed before this analysis was run.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT"); csd <- file.path(root, "data/02_processed/crossspecies"); ext <- file.path(root, "data/03_external/fay2023")
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu"); METHODS <- c("DESeq2","edgeR","limma"); PRIOR <- 2
partner <- c(ScxSp="SkxSu",SkxSu="ScxSp",ScxSu="SpxSk",SpxSk="ScxSu",ScxSk="SpxSu",SpxSu="ScxSk")
args <- commandArgs(TRUE); KS <- if (length(args)) as.integer(args) else c(100L, 50L, 200L, 0L); NB <- 2000L
X <- fread(cmd = paste("zcat", shQuote(file.path(csd, "fay_multimethod_p1order.tsv.gz"))))[condition == "heat" & group %in% HYB]   # SpxSk refitted with the first-named allele as reference
Rtab <- X[model == "A" & term == "time", .(gene, group, method, R = eff)]
inst <- function(g, mt) { keep <- setdiff(HYB, c(g, partner[[g]]))
  Rtab[group %in% keep & method == mt, .(S = mean(R, na.rm = TRUE), nS = sum(is.finite(R))), by = gene][nS >= 3] }
## block responses for h (script 94 h_arm_on, 15 min)
D <- as.data.table(read.csv(file.path(ext, "design.csv"))); D[, col := Name]
D[, num := as.integer(sub("^[A-Za-z]+", "", Name))]; D[, block := fifelse(num <= 32L, "block1", "block2")]
Hd <- D[Strain %in% HYB & Condition == "heat"]; Hd[, Time := as.integer(Time)]
dat <- read.csv(file.path(ext, "htseq.ortho.csv")); one <- dat$one == 1
cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
blocks <- function(g) { hs <- Hd[Strain == g][order(Allele, block, Time)]; a1 <- substr(g, 1, 2)
  lc <- cpm(calcNormFactors(DGEList(counts = cts[, hs$col, drop = FALSE]), method = "TMM"), log = TRUE, prior.count = PRIOR)
  gv <- function(bk, t) lc[, hs[Allele == a1 & block == bk & Time == t, col]]
  data.table(gene = rownames(lc), b1 = gv("block1", 15L) - gv("block1", 0L), b2 = gv("block2", 15L) - gv("block2", 0L)) }
## per hybrid data: one row per gene, method-specific e, R, S
DAT <- lapply(setNames(HYB, HYB), function(g) {
  L <- lapply(METHODS, function(mt) { ii <- X[group == g & method == mt & model == "A" & term == "inter", .(gene, e = eff)]
    d <- merge(merge(ii, Rtab[group == g & method == mt, .(gene, R)], by = "gene"), inst(g, mt), by = "gene")
    d <- d[is.finite(e) & is.finite(R) & is.finite(S), .(gene, e, R, S)]; setnames(d, c("e","R","S"), paste0(c("e","R","S"), "_", mt)); d })
  d <- Reduce(function(a, b) merge(a, b, by = "gene"), L); d <- merge(d, blocks(g), by = "gene"); d[is.finite(b1) & is.finite(b2)] })
## clusters on R across hybrids (limma)
W <- dcast(Rtab[method == "limma"], gene ~ group, value.var = "R")
for (h in HYB) { v <- W[[h]]; v[!is.finite(v)] <- mean(v, na.rm = TRUE); W[[h]] <- (v - mean(v)) / sd(v) }
genes_all <- unique(unlist(lapply(DAT, `[[`, "gene"))); W <- W[gene %in% genes_all]
wcov <- function(x, y, w) { mx <- sum(w * x) / sum(w); my <- sum(w * y) / sum(w); sum(w * (x - mx) * (y - my)) / sum(w) }
stat <- function(d, w) {
  sig <- (d$b1 + d$b2) / 2; nul <- (d$b2 - d$b1) / 2; h <- 1 - wcov(nul, nul, w) / wcov(sig, sig, w)
  per <- sapply(METHODS, function(mt) { e <- d[[paste0("e_", mt)]]; R <- d[[paste0("R_", mt)]]; S <- d[[paste0("S_", mt)]]
    c(naive = 1 + wcov(e, R, w) / wcov(R, R, w), iv = 1 + wcov(e, S, w) / wcov(R, S, w)) })
  nv <- median(per["naive", ]); iv <- median(per["iv", ]); c(naive = nv, iv = iv, h = h, rel = nv / h, closure = nv - iv * h) }
OUT <- list()
for (K in KS) { set.seed(20261004L)
  cl <- if (K == 0L) setNames(seq_len(nrow(W)), W$gene) else setNames(kmeans(as.matrix(W[, ..HYB]), centers = K, nstart = 25, iter.max = 100)$cluster, W$gene)
  ncl <- max(cl)
  BOOT <- replicate(NB, { m <- tabulate(sample.int(ncl, ncl, TRUE), ncl)
    sapply(HYB, function(g) { d <- DAT[[g]]; stat(d, m[cl[d$gene]]) }) }, simplify = "array")   # 5 x 6 x NB
  for (g in HYB) { d <- DAT[[g]]; p <- stat(d, rep(1, nrow(d)))
    q <- apply(BOOT[, g, ], 1, quantile, c(.025, .975), na.rm = TRUE)
    OUT[[length(OUT) + 1]] <- data.table(k = if (K == 0L) "genes" else as.character(K), group = g, n = nrow(d),
      naive = p["naive"], iv = p["iv"], iv_lo = q[1, "iv"], iv_hi = q[2, "iv"], h = p["h"], h_lo = q[1, "h"], h_hi = q[2, "h"],
      rel = p["rel"], rel_lo = q[1, "rel"], rel_hi = q[2, "rel"], closure = p["closure"], cl_lo = q[1, "closure"], cl_hi = q[2, "closure"]) }
  cat("done k =", K, "\n") }
O <- rbindlist(OUT); fwrite(O, file.path(csd, "p1_table1_joint_bootstrap.tsv"), sep = "\t")
print(O[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
for (K in unique(O$k)) { z <- O[k == K]; cat(sprintf("k=%s: IV resolved below 1 in %d of 6; naive/h resolved below 1 in %d of 6; closure excludes 0 in %d of 6\n",
  K, sum(z$iv_hi < 1), sum(z$rel_hi < 1), sum(z$cl_lo > 0 | z$cl_hi < 0))) }
