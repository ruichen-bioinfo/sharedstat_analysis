#!/usr/bin/env Rscript
## 74_p2_joint_candidates.R -- Supplementary Table S12 of P2, the joint-requirement candidates.
## A gene is listed when any of the three differential-expression frameworks places it in the redundantly controlled
## class with both single-perturbation effects negative (sub = red_joint, the rule of 01_epistasis_decomposition.R).
## Per method: the class call, the sub-label, and the two single-perturbation log2 fold changes and adjusted p. Support
## that never entered the class definition: HSE and STRE in the promoter (the full scan, p2x_full_scan_per_gene.tsv,
## the one P2 reports) and YEASTRACT records for Hsf1 and Msn2/4, direct and any evidence, as in
## 46_p2_sequence_orthogonal_validation.R. Every row is a candidate; nothing here demonstrates joint requirement.
## Reads stored tables only; no statistic is recomputed.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
stopifnot(all(c("gene", "class", "sub", "method", "lfc_kd", "padj_kd", "lfc_ko", "padj_ko", "lfc_db", "padj_db") %in% names(CL)))
## verify the per-method partition counts
cnt <- CL[class == "redundant", .N, by = .(method, sub)]
want <- data.table(method = rep(c("DESeq2", "edgeR", "limma"), each = 3),
                   sub = rep(c("red_Msn24_subthreshold", "red_joint", "red_unclassified"), 3),
                   N = c(96L, 43L, 18L, 89L, 45L, 13L, 96L, 51L, 13L))
chk <- merge(want, cnt, by = c("method", "sub"), suffixes = c("_want", "_got")); stopifnot(nrow(chk) == 9L, all(chk$N_want == chk$N_got))
cat("Class partition reproduced in all 9 method x sub-label cells\n")
cand <- unique(CL[sub == "red_joint", gene])
W <- CL[gene %in% cand, .(gene, method, class, sub, lfc_kd, padj_kd, lfc_ko, padj_ko, lfc_db, padj_db)]
Wd <- dcast(W, gene ~ method, value.var = c("sub", "lfc_kd", "padj_kd", "lfc_ko", "padj_ko", "lfc_db", "padj_db"))
Wd[, n_methods_joint := rowSums(as.matrix(.SD) == "red_joint", na.rm = TRUE), .SDcols = patterns("^sub_")]
## unanimous exclusive class label: one class across all three methods, otherwise "unassigned"
U <- CL[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][, exclusive_class := ifelse(k == 1L & nm == 3L, cls, "unassigned")]
Wd <- merge(Wd, U[, .(gene, exclusive_class)], by = "gene", all.x = TRUE)
## promoter support
P <- fread(file.path(hsd, "p2x_full_scan_per_gene.tsv"))[, .(gene = sys, hse_classic = n_classic, hse_gapped = n_gapped, hse_any, n_stre, stre_per_kb)]
Wd <- merge(Wd, P, by = "gene", all.x = TRUE)
## curated targets
Y <- fread(file.path(root, "data/03_reference/yeastract2019-flat-file.tsv"), header = FALSE, quote = "", sep = "\t", select = c(2, 3, 9),
           col.names = c("tf", "tgt", "evidence"))[, tf := toupper(tf)]
tg <- function(tfs, ev = NULL) unique(Y[tf %in% tfs & (is.null(ev) | evidence %in% ev), tgt])
Wd[, `:=`(yeastract_hsf1_direct = gene %in% tg("HSF1", "Direct"), yeastract_hsf1_any = gene %in% tg("HSF1"),
          yeastract_msn24_direct = gene %in% tg(c("MSN2", "MSN4"), "Direct"), yeastract_msn24_any = gene %in% tg(c("MSN2", "MSN4")))]
Wd[, status := "candidate for joint requirement, not demonstrated"]
setorder(Wd, -n_methods_joint, gene)
fwrite(Wd, file.path(hsd, "p2_joint_requirement_candidates.tsv"), sep = "\t")
cat(sprintf("%d candidate genes: %d called red_joint by all three methods, %d by two, %d by one\n", nrow(Wd),
            sum(Wd$n_methods_joint == 3), sum(Wd$n_methods_joint == 2), sum(Wd$n_methods_joint == 1)))
cat(sprintf("promoter scanned: %d; HSE in %d; STRE >= 1 in %d; YEASTRACT Hsf1 direct %d, Msn2/4 direct %d; exclusive class 'redundant' %d\n",
            sum(!is.na(Wd$hse_any)), sum(Wd$hse_any %in% TRUE), sum(Wd$n_stre > 0, na.rm = TRUE),
            sum(Wd$yeastract_hsf1_direct), sum(Wd$yeastract_msn24_direct), sum(Wd$exclusive_class == "redundant", na.rm = TRUE)))
