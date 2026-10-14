#!/usr/bin/env Rscript
## 103b_p2_torc1_panel.R -- follow-up TORC1 check at 25 C, fixed before this analysis was run:
## limma log2 fold changes with 95% intervals for the panel GAP1, MEP2, DAL80, GDH2, DAL5, CIT2, and a per-library
## ribosomal-protein score (Welch t interval), for H-W 25, D-M 25 and wild type 37/10 vs 25.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT"); hsd <- file.path(root, "data/02_processed/hsf1"); pro <- file.path(root, "data/02_processed")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(root, "data/01_metadata/sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + factor(SM$grp)); colnames(des) <- levels(factor(SM$grp))
v <- voom(dge, des, plot = FALSE); fit <- lmFit(v, des)
ct <- function(a, b) { x <- setNames(rep(0, ncol(des)), colnames(des)); x[a] <- 1; x[b] <- -1; x }
CT <- cbind(`H-W 25` = ct("HSF1.KD__25_0", "Wildtype__25_0"), `D-M 25` = ct("Double.KDKO__25_0", "MSN24.KO__25_0"),
            `WT 37/10 vs 25` = ct("Wildtype__37_10", "Wildtype__25_0"))
f2 <- eBayes(contrasts.fit(fit, CT))
gl <- grep("\tgene\t", readLines(file.path(root, "data/03_reference/Saccharomyces_cerevisiae.R64-1-1.75.gtf"), warn = FALSE), value = TRUE)
nm <- data.table(sys = sub('.*gene_id "([^"]+)".*', "\\1", gl), std = ifelse(grepl('gene_name "', gl), sub('.*gene_name "([^"]+)".*', "\\1", gl), NA))
key <- rownames(f2); std <- nm$std[match(key, nm$sys)]
PANEL <- c("GAP1", "MEP2", "DAL80", "GDH2", "DAL5", "CIT2")
OUT <- rbindlist(lapply(colnames(CT), function(cn) { tt <- topTable(f2, coef = cn, number = Inf, sort.by = "none", confint = TRUE)
  i <- which(std %in% PANEL); data.table(gene = std[i], contrast = cn, log2fc = tt$logFC[i], lo = tt$CI.L[i], hi = tt$CI.R[i], adj_p = tt$adj.P.Val[i]) }))
cat("panel genes found:", paste(sort(unique(OUT$gene)), collapse = " "), "\n")
## per-library RP score
E <- v$E; rp <- key[!is.na(std) & grepl("^RP[LS][0-9]", std)]
score <- colMeans(E[rp, ]) - colMeans(E)
sc <- function(a, b) { x <- score[SM$grp == a]; y <- score[SM$grp == b]; t <- t.test(x, y); c(diff = mean(x) - mean(y), lo = t$conf.int[1], hi = t$conf.int[2]) }
RPS <- rbind(`H-W 25` = sc("HSF1.KD__25_0", "Wildtype__25_0"), `D-M 25` = sc("Double.KDKO__25_0", "MSN24.KO__25_0"), `WT 37/10 vs 25` = sc("Wildtype__37_10", "Wildtype__25_0"))
RPS <- data.table(gene = "RP library score", contrast = rownames(RPS), log2fc = RPS[, 1], lo = RPS[, 2], hi = RPS[, 3], adj_p = NA_real_)
OUT <- rbind(OUT, RPS); fwrite(OUT, file.path(hsd, "p2_torc1_panel.tsv"), sep = "\t")
print(OUT[order(gene, contrast), .(gene, contrast, l2fc = round(log2fc, 2), lo = round(lo, 2), hi = round(hi, 2), adj_p = signif(adj_p, 2))])
hw <- OUT[contrast == "H-W 25" & gene %in% PANEL]; rpc <- RPS[contrast == "H-W 25"]
closed <- !any(hw$log2fc >= 1 & hw$lo > 0) && RPS[contrast == "H-W 25", all(lo > -0.58)]   # the criterion uses the H-W contrast only
cat("DECISION:", if (closed) "TORC1 route closed is supported" else "not supported", "\n")
