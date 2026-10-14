#!/usr/bin/env Rscript
## 79_p2_replicate_jackknife.R -- a diagnostic independent of the preceding diagnostic: the uncertainty of a class reading when
## the unit of replication is the replicate, not the gene. Delete-one-replicate jackknife: drop replicate r from every
## genotype x condition group (r = 1, 2, 3), refit the same voom group-means model (published normalisation), and
## recompute each class median of I(37_30) and of I(37_30) - I(25) minus the not-induced median. Jackknife SE =
## sqrt((k - 1)/k * sum((theta_r - mean)^2)), k = 3, which is crude with three units and is reported as such. Compared
## with the gene-independence SE of Table 1 (median(se)/sqrt(n)).
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
keep0 <- filterByExpr(DGEList(M, group = factor(SM$grp)), group = SM$grp)       # same gene set in every refit
CLS3 <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
U <- CLS3[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][, cls := ifelse(k == 1L & nm == 3L, cls, "unassigned")]
G4 <- c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"); SG <- c(1, -1, -1, 1)
fitI <- function(use) {
  d <- calcNormFactors(DGEList(M[keep0, use], group = factor(SM$grp[use])))
  grp <- factor(SM$grp[use]); des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
  fit <- lmFit(voom(d, des, plot = FALSE), des)
  mk <- function(c_) { ct <- setNames(rep(0, ncol(des)), colnames(des)); ct[paste(G4, c_, sep = "__")] <- SG; ct }
  cf <- contrasts.fit(fit, cbind(lvl = mk("37_30"), rsp = mk("37_30") - mk("25_0")))
  X <- data.table(gene = rownames(cf$coefficients), lvl = cf$coefficients[, "lvl"], rsp = cf$coefficients[, "rsp"],
                  se_lvl = cf$stdev.unscaled[, "lvl"] * cf$sigma, se_rsp = cf$stdev.unscaled[, "rsp"] * cf$sigma)
  X <- merge(X, U[, .(gene, cls)], by = "gene")
  S <- X[, .(n = .N, lvl = median(lvl), rsp = median(rsp), se_lvl = median(se_lvl) / sqrt(.N), se_rsp = median(se_rsp) / sqrt(.N)), by = cls]
  S[, `:=`(lvl_ref = lvl - S[cls == "not_induced", lvl], rsp_ref = rsp - S[cls == "not_induced", rsp])]; S }
FULL <- fitI(rep(TRUE, nrow(SM)))
## fidelity: the full fit reproduces Table 1's published vs-ref values
T1 <- fread(file.path(hsd, "p2_class_exclusive.tsv"))[regime == "all genes (as published)"]
chk <- merge(FULL[, .(cls, lvl_ref)], T1[, .(cls, vsref)], by = "cls"); stopifnot(max(abs(chk$lvl_ref - chk$vsref)) < 1e-8)
cat("fidelity: full refit reproduces Table 1 vs-ref values\n")
JK <- rbindlist(lapply(1:3, function(r) fitI(SM$replicate != r)[, drop := r]))
S <- JK[, .(jk_se_lvl = sqrt(2/3 * sum((lvl_ref - mean(lvl_ref))^2)), jk_se_rsp = sqrt(2/3 * sum((rsp_ref - mean(rsp_ref))^2)),
            lvl_ref_range = paste(sprintf("%+.3f", range(lvl_ref)), collapse = " to "),
            rsp_ref_range = paste(sprintf("%+.3f", range(rsp_ref)), collapse = " to ")), by = cls]
S <- merge(FULL[, .(cls, n, lvl_ref, rsp_ref, se_lvl, se_rsp)], S, by = "cls")
S[, `:=`(z_lvl_gene = lvl_ref / se_lvl, z_lvl_jk = lvl_ref / jk_se_lvl, z_rsp_gene = rsp_ref / se_rsp, z_rsp_jk = rsp_ref / jk_se_rsp,
         se_ratio_lvl = jk_se_lvl / se_lvl)]
fwrite(S, file.path(hsd, "p2_replicate_jackknife.tsv"), sep = "\t")
print(S[order(-n), .(cls, n, lvl_ref = round(lvl_ref, 3), z_gene = round(z_lvl_gene, 1), z_jk = round(z_lvl_jk, 1), lvl_range = lvl_ref_range,
                     rsp_ref = round(rsp_ref, 3), z_rsp_gene = round(z_rsp_gene, 1), z_rsp_jk = round(z_rsp_jk, 1), rsp_range = rsp_ref_range,
                     se_ratio = round(se_ratio_lvl, 1))])
