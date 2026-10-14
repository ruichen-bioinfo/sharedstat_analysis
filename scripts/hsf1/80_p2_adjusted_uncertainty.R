#!/usr/bin/env Rscript
## 80_p2_adjusted_uncertainty.R -- correlation-adjusted uncertainty for P2, fixed before this analysis was run.
## Point estimates are read from the saved tables and validated before uncertainty is estimated; only uncertainty changes.
## (1) class vs reference: SE^2 = VIF_c se_c^2 + VIF_r se_r^2, se = median(gene se)/sqrt(n) as in Table 1, VIF = 1 +
##    (n - 1) rho, rho = limma::interGeneCorrelation of the group-means model residuals within the class (<= 2000 genes
##    drawn when larger), per normalisation regime; camera p for each class (all genes as the competing set, the
##    contrast as the test) beside it; the delete-one-replicate jackknife of 79 is joined as a sensitivity column.
##    The same adjustment for the genome-wide median (rho over all retained genes).
## (2) presence: no pairing at 0, every leave-one-replicate-out amplitude > 0 (with three replicates, leaving one out
##    leaves exactly one pairing, so this is the same set of numbers as the pairings), and the bootstrap lower bound
##    inflated as A - sqrt(VIF_g) (A - lo) > 0, VIF_g from the genome-wide mean residual correlation.
## (3) 37/42 ratio interval inflated the same way on the ratio scale.  (4) h3 = SNR^2 / (SNR^2 + 2/3).
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
say <- function(...) cat(sprintf(...), "\n", sep = "")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
CLS3 <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
U <- CLS3[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][, cls := ifelse(k == 1L & nm == 3L, cls, "unassigned")]
RESP <- unique(CLS3[class %in% c("Hsf1_only", "Msn24_only", "redundant", "both_required"), gene])   # as script 20
G4 <- c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"); SG <- c(1, -1, -1, 1)
fit_regime <- function(norm_genes) {
  dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  sub <- dge[rownames(dge) %in% norm_genes, , keep.lib.sizes = TRUE]
  dge$samples$norm.factors <- calcNormFactors(sub)$samples$norm.factors
  grp <- factor(SM$grp); des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
  v <- voom(dge, des, plot = FALSE); fit <- lmFit(v, des)
  mk <- function(c_) { ct <- setNames(rep(0, ncol(des)), colnames(des)); ct[paste(G4, c_, sep = "__")] <- SG; ct }
  CT <- cbind(`I(37_30)` = mk("37_30"), `I(25)` = mk("25_0"), `I(37_30) - I(25)` = mk("37_30") - mk("25_0"))
  list(v = v, des = des, CT = CT, cf = contrasts.fit(fit, CT))
}
ALLG <- rownames(M)
REG <- list(`all genes (as published)` = fit_regime(ALLG),
            `response classes excluded from normalisation` = fit_regime(setdiff(ALLG, RESP)))
rho_of <- function(v, des, g) { g <- intersect(g, rownames(v)); if (length(g) > 2000) { set.seed(1); g <- sample(g, 2000) }
  interGeneCorrelation(v$E[g, ], des)$correlation }
OUT <- list(); GW <- list(); CAM <- list()
for (rg in names(REG)) { R <- REG[[rg]]
  RH <- rbindlist(lapply(unique(U$cls), function(cl) data.table(cls = cl, rho = rho_of(R$v, R$des, U[cls == cl, gene]))))
  rho_g <- rho_of(R$v, R$des, rownames(R$v))
  for (ct in colnames(R$CT)) {
    X <- data.table(gene = rownames(R$cf$coefficients), I = R$cf$coefficients[, ct], se = R$cf$stdev.unscaled[, ct] * R$cf$sigma)
    GW[[length(GW) + 1]] <- data.table(regime = rg, contrast = ct, n = nrow(X), med = median(X$I), se_gene = median(X$se) / sqrt(nrow(X)),
                                       rho = rho_g, vif = 1 + (nrow(X) - 1) * rho_g)
    X <- merge(X, U[, .(gene, cls)], by = "gene")
    S <- X[, .(n = .N, med = median(I), se_gene = median(se) / sqrt(.N)), by = cls]
    S <- merge(S, RH, by = "cls")[, vif := 1 + (n - 1) * pmax(rho, 0)]
    r <- S[cls == "not_induced"]
    S[, `:=`(regime = rg, contrast = ct, vs_ref = med - r$med,
             z_gene = (med - r$med) / se_gene,
             se_adj = sqrt(vif * se_gene^2 + r$vif * r$se_gene^2))]
    S[cls == "not_induced", se_adj := NA_real_]
    S[, `:=`(z0_gene = med / se_gene, z0_adj = med / (sqrt(vif) * se_gene))]
    S[, `:=`(z_adj = vs_ref / se_adj, lo_adj = vs_ref - 1.96 * se_adj, hi_adj = vs_ref + 1.96 * se_adj)]
    idx <- lapply(setdiff(unique(U$cls), "not_induced"), function(cl) which(rownames(R$v) %in% U[cls == cl, gene]))
    names(idx) <- setdiff(unique(U$cls), "not_induced")
    ## inter.gene.cor = NA: camera estimates each set's own correlation; its default is a fixed 0.01, which is not ours
    cam <- camera(R$v, idx, R$des, contrast = R$CT[, ct], inter.gene.cor = NA, use.ranks = FALSE)
    CAM[[length(CAM) + 1]] <- data.table(regime = rg, contrast = ct, cls = rownames(cam), camera_direction = cam$Direction, camera_p = cam$PValue)
    OUT[[length(OUT) + 1]] <- S
  }
}
S <- rbindlist(OUT); CAMT <- rbindlist(CAM); S <- merge(S, CAMT, by = c("regime", "contrast", "cls"), all.x = TRUE)
GWT <- rbindlist(GW)[, `:=`(z_gene = med / se_gene, se_adj = sqrt(vif) * se_gene)][, z_adj := med / se_adj]
## validation: the published-regime, anchor, level-scale point estimates and gene z equal Table 1 (p2_class_exclusive)
T1 <- fread(file.path(hsd, "p2_class_exclusive.tsv"))
fz <- merge(S[contrast == "I(37_30)", .(regime, cls, vs_ref, z_gene)], T1[, .(regime, cls, vsref, zref)], by = c("regime", "cls"))
stopifnot(nrow(fz) == 14L, max(abs(fz$vs_ref - fz$vsref)) < 1e-10, max(abs(fz$z_gene - fz$zref)) < 1e-8)
LS <- fread(file.path(hsd, "p2_offset_origin_summary.tsv"))[condition == "37_30"]
stopifnot(max(abs(GWT[contrast == "I(37_30)"][order(regime), med] - LS[order(regime), med_I])) < 1e-10)
say("validation: 14 Table 1 point estimates and gene-independence z reproduced; 2 genome-wide anchor medians reproduced")
## jackknife sensitivity from script 79 under the primary normalization regime
JK <- fread(file.path(hsd, "p2_replicate_jackknife.tsv"))
S <- merge(S, rbind(JK[, .(regime = "all genes (as published)", contrast = "I(37_30)", cls, z_jackknife = z_lvl_jk)],
                    JK[, .(regime = "all genes (as published)", contrast = "I(37_30) - I(25)", cls, z_jackknife = z_rsp_jk)]),
           by = c("regime", "contrast", "cls"), all.x = TRUE)
fwrite(S, file.path(hsd, "p2_class_adjusted.tsv"), sep = "\t"); fwrite(GWT, file.path(hsd, "p2_genomewide_adjusted.tsv"), sep = "\t")
## (2)-(4) on the amplitude tables of 78
BP <- fread(file.path(hsd, "p2_detectability_bypairing.tsv"))
rho_g <- GWT[regime == "all genes (as published)" & contrast == "I(37_30)", rho]; VG <- 1 + (nrow(REG[[1]]$v) - 1) * rho_g
BP[, lo_adj := pmax(0, A - sqrt(VG) * (A - lo))]
Q1 <- BP[, .(n_pairings = .N, zero_pairings = sum(A == 0), loro_all_positive = all(A > 0), min_lo_adj = min(lo_adj),
             present = all(A > 0) & all(lo_adj > 0)), by = .(scale, condition)]
## (3) A ratio with a zero denominator has no interval to inflate, so the lower bound of the ratio is taken
## conservatively as the inflated lower bound of the 37 °C amplitude over the inflated upper bound of the 42 °C
## amplitude, same pairing (same labels): lo_ratio = lo_adj(37) / hi_adj(42), hi_adj = A + sqrt(VIF) (hi - A).
## This interval construction it is more conservative than inflating the ratio.
BP[, hi_adj := A + sqrt(VG) * (hi - A)]
RB <- rbindlist(lapply(c("level", "response"), function(sc) rbindlist(lapply(c("10", "30"), function(tp) {
  a <- BP[scale == sc & condition == paste0("37_", tp)]; b <- BP[scale == sc & condition == paste0("42_", tp)]
  m <- merge(a[, .(labels_25, pairing, A37 = A, lo37 = lo_adj)], b[, .(labels_25, pairing, A42 = A, hi42 = hi_adj)], by = c("labels_25", "pairing"))
  stopifnot(nrow(m) == nrow(a))
  m[, `:=`(scale = sc, timepoint = tp, lo_ratio = ifelse(hi42 > 0, lo37 / hi42, ifelse(lo37 > 0, Inf, 0)))] }))))
Q2 <- RB[, .(n_pairings = .N, min_lo_ratio = min(lo_ratio), holds = all(lo_ratio > 1)), by = .(scale, timepoint)]
fwrite(RB, file.path(hsd, "p2_ratio_adjusted_bypairing.tsv"), sep = "\t")
Q3 <- fread(file.path(hsd, "p2_detectability_readability.tsv"))[, readable_h3 := h3 >= 0.5]
fwrite(Q1, file.path(hsd, "p2_presence_adjusted.tsv"), sep = "\t"); fwrite(Q2, file.path(hsd, "p2_ratio_adjusted.tsv"), sep = "\t")
say("=== 80_p2_adjusted_uncertainty | %s ===", format(Sys.time()))
say("genome-wide residual correlation %.4f, VIF %.1f", rho_g, VG)
say("\nGENOME-WIDE MEDIANS"); print(GWT[, .(regime = substr(regime, 1, 12), contrast, med = round(med, 4), z_gene = round(z_gene, 1), vif = round(vif, 1), z_adj = round(z_adj, 2))])
say("\nCLASSES vs not-induced"); print(S[cls != "not_induced"][order(regime, contrast, -n), .(regime = substr(regime, 1, 12), contrast, cls, n, vs_ref = round(vs_ref, 3), z_gene = round(z_gene, 1), vif = round(vif, 1), z_adj = round(z_adj, 2), lo = round(lo_adj, 3), hi = round(hi_adj, 3), camera_p = signif(camera_p, 2), z_jk = round(z_jackknife, 1))])
say("\nU2 PRESENCE"); print(Q1); say("\nU3 RATIO"); print(Q2); say("\nU4 READABILITY"); print(Q3[, .(scale, condition, snr = round(snr, 3), h3 = round(h3, 3), readable_h3)])
