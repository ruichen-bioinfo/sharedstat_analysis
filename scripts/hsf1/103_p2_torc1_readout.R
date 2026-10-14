#!/usr/bin/env Rscript
## 103_p2_torc1_readout.R -- TORC1 readout at 25 C: median log2 fold change, with a 2000-draw gene bootstrap interval,
## of ribosomal-protein genes (RPL/RPS), ribosome biogenesis (GO:0042254), the Msn2/4-only class and all genes, for
## H-W 25, D-M 25 and wild type 37/10 and 37/30 vs 25. The decision thresholds (no signature: |median| <= 0.25 and
## lower bound > -0.58 in both RP contrasts; signature: median < -0.58 with upper bound < 0) were fixed before this analysis was run.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata"); ref <- file.path(root, "data/03_reference")
set.seed(20261003L); NB <- 2000L
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + factor(SM$grp)); colnames(des) <- levels(factor(SM$grp))
fit <- lmFit(voom(dge, des, plot = FALSE), des)
ct <- function(a, b) { v <- setNames(rep(0, ncol(des)), colnames(des)); v[a] <- 1; v[b] <- -1; v }
CT <- cbind(`H-W 25` = ct("HSF1.KD__25_0", "Wildtype__25_0"), `D-M 25` = ct("Double.KDKO__25_0", "MSN24.KO__25_0"),
            `WT 37/10 vs 25` = ct("Wildtype__37_10", "Wildtype__25_0"), `WT 37/30 vs 25` = ct("Wildtype__37_30", "Wildtype__25_0"))
cf <- contrasts.fit(fit, CT)$coefficients
## gene sets
gl <- grep("\tgene\t", readLines(file.path(ref, "Saccharomyces_cerevisiae.R64-1-1.75.gtf"), warn = FALSE), value = TRUE)
nm <- data.table(sys = sub('.*gene_id "([^"]+)".*', "\\1", gl), std = ifelse(grepl('gene_name "', gl), sub('.*gene_name "([^"]+)".*', "\\1", gl), NA))
key <- rownames(cf); sysk <- if (mean(key %in% nm$sys) > 0.5) key else nm$sys[match(toupper(key), toupper(nm$std))]
RP <- key[sysk %in% nm[grepl("^RP[LS][0-9]", std), sys]]
GO <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "go_annotations_propagated_universe.tsv.gz"))))
gocol <- intersect(names(GO), c("go", "go_id", "GO", "term"))[1]; scol <- intersect(names(GO), c("sys", "gene", "systematic"))[1]
RIBI <- key[sysk %in% GO[get(gocol) == "GO:0042254", get(scol)]]
CLS <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
U <- CLS[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][k == 1L & nm == 3L]
MSN <- intersect(key, U[cls == "Msn24_only", gene])
SETS <- list(`ribosomal protein (RPL/RPS)` = RP, `ribosome biogenesis (GO:0042254)` = RIBI, `Msn2/4-only class` = MSN, `all genes` = key)
med_ci <- function(x) { b <- replicate(NB, median(sample(x, length(x), TRUE))); c(median(x), quantile(b, c(.025, .975))) }
OUT <- rbindlist(lapply(names(SETS), function(s) rbindlist(lapply(colnames(CT), function(cn) {
  v <- med_ci(cf[SETS[[s]], cn]); data.table(set = s, n = length(SETS[[s]]), contrast = cn, median_log2fc = v[1], lo = v[2], hi = v[3]) }))))
fwrite(OUT, file.path(hsd, "p2_torc1_readout.tsv"), sep = "\t")
print(OUT[, .(set, n, contrast, med = round(median_log2fc, 3), lo = round(lo, 3), hi = round(hi, 3))])
rp <- OUT[set == "ribosomal protein (RPL/RPS)" & contrast %in% c("H-W 25", "D-M 25")]
verdict <- if (all(abs(rp$median_log2fc) <= 0.25 & rp$lo > -0.58)) "no detectable TORC1 signature" else
           if (any(rp$median_log2fc < -0.58 & rp$hi < 0)) "TORC1 signature present" else "inconclusive"
cat("DECISION:", verdict, "\n")
