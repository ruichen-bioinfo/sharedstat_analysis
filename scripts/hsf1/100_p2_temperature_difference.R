#!/usr/bin/env Rscript
## 100_p2_temperature_difference.R -- temperature difference of the P2 interaction, with the decision
## fixed before this analysis was run.
## D(t) = I(37,t) - I(42,t). B1 genome-wide amplitude by replicate pairing (scripts 78/80 construction, 42 C labels
## permuted 6 ways); B2 class medians vs not-induced with the correlation-adjusted interval of script 80.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
NB <- 2000L; set.seed(20261002L)
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + factor(SM$grp)); colnames(des) <- levels(factor(SM$grp))
v <- voom(dge, des, plot = FALSE); E <- v$E; G <- nrow(E)
G4 <- c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"); SG <- c(1, -1, -1, 1)
col_of <- function(g, c, r) which(SM$knockout == g & SM$cond == c & SM$replicate == r)
I_rep <- function(c, r) as.vector(E[, vapply(G4, function(g) col_of(g, c, r), 0L)] %*% SG)
IR <- lapply(c("37_10","37_30","42_10","42_30"), function(c_) sapply(1:3, function(r) I_rep(c_, r))); names(IR) <- c("37_10","37_30","42_10","42_30")
PR <- combn(1:3, 2); perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
cols <- list(); MT <- list()
for (tp in c("10","30")) for (pk in perms) for (j in 1:3) { a <- PR[1, j]; b <- PR[2, j]
  X <- IR[[paste0("37_", tp)]]; Z <- IR[[paste0("42_", tp)]][, pk]
  s <- (X[, a] + X[, b]) / 2 - (Z[, a] + Z[, b]) / 2; n <- (X[, b] - X[, a]) / 2 - (Z[, b] - Z[, a]) / 2
  cols[[length(cols) + 1]] <- s^2 - n^2
  MT[[length(MT) + 1]] <- data.table(timepoint = tp, labels_42 = paste(pk, collapse = ""), pairing = sprintf("r%d-r%d", a, b)) }
MT <- rbindlist(MT); Dm <- do.call(cbind, cols)
amp <- function(w) sqrt(pmax(0, colSums(Dm * w) / sum(w)))
A0 <- amp(rep(1, G)); BOOT <- t(replicate(NB, amp(tabulate(sample.int(G, G, replace = TRUE), G))))
rho_g <- fread(file.path(hsd, "p2_genomewide_adjusted.tsv"))[regime == "all genes (as published)" & contrast == "I(37_30)", rho]
VG <- 1 + (G - 1) * rho_g
MT[, `:=`(A = A0, lo = apply(BOOT, 2, quantile, 0.025), hi = apply(BOOT, 2, quantile, 0.975))][, lo_adj := pmax(0, A - sqrt(VG) * (A - lo))]
B1 <- MT[, .(n_pairings = .N, zero_pairings = sum(A == 0), min_A = min(A), median_A = median(A), min_lo_adj = min(lo_adj),
             present = all(A > 0) & all(lo_adj > 0)), by = timepoint]
## B2: class medians vs not-induced
CLS3 <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
U <- CLS3[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][, cls := ifelse(k == 1L & nm == 3L, cls, "unassigned")]
mk <- function(c_) { ct <- setNames(rep(0, ncol(des)), colnames(des)); ct[paste(G4, c_, sep = "__")] <- SG; ct }
CT <- cbind(`D(10)` = mk("37_10") - mk("42_10"), `D(30)` = mk("37_30") - mk("42_30"))
cf <- contrasts.fit(lmFit(v, des), CT)
rho_of <- function(g) { g <- intersect(g, rownames(v)); if (length(g) > 2000) { set.seed(1); g <- sample(g, 2000) }
  interGeneCorrelation(v$E[g, ], des)$correlation }
RH <- rbindlist(lapply(unique(U$cls), function(cl) data.table(cls = cl, rho = rho_of(U[cls == cl, gene]))))
B2 <- rbindlist(lapply(colnames(CT), function(ct) {
  X <- merge(data.table(gene = rownames(cf$coefficients), D = cf$coefficients[, ct], se = cf$stdev.unscaled[, ct] * cf$sigma), U[, .(gene, cls)], by = "gene")
  S <- merge(X[, .(n = .N, med = median(D), se_gene = median(se) / sqrt(.N)), by = cls], RH, by = "cls")[, vif := 1 + (n - 1) * pmax(rho, 0)]
  r <- S[cls == "not_induced"]
  S[, `:=`(contrast = ct, vs_ref = med - r$med, se_adj = sqrt(vif * se_gene^2 + r$vif * r$se_gene^2))]
  S[, `:=`(z_adj = vs_ref / se_adj, lo_adj = vs_ref - 1.96 * se_adj, hi_adj = vs_ref + 1.96 * se_adj)][cls != "not_induced"] }))
fwrite(MT, file.path(hsd, "p2_tempdiff_amplitude_bypairing.tsv"), sep = "\t"); fwrite(B1, file.path(hsd, "p2_tempdiff_amplitude.tsv"), sep = "\t")
fwrite(B2, file.path(hsd, "p2_tempdiff_classes.tsv"), sep = "\t")
cat(sprintf("genome-wide rho %.4f VIF %.1f\n", rho_g, VG)); print(B1)
print(B2[order(contrast, -n), .(contrast, cls, n, vs_ref = round(vs_ref, 3), z_adj = round(z_adj, 2), lo = round(lo_adj, 3), hi = round(hi_adj, 3))])
d <- B2[contrast == "D(30)" & cls == "Msn24_only"]
cat(sprintf("RESULT: Msn2/4-only D(30) interval [%.3f, %.3f] -> %s\n", d$lo_adj, d$hi_adj,
            if (d$lo_adj > 0 || d$hi_adj < 0) "excludes 0: temperature-dependent difference resolved" else "includes 0: difference not resolved"))
