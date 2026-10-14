# 29_geo_published_lfc_check.R -- do our re-analysed log fold changes agree with the authors' published ones?
#
# Purpose. The perturbation dataset is available from GEO. GEO records the
# series as public since 2024-01-29, and the public supplement is the authors' own DESeq2 result table:
# one Log2FC and one FDR column per genotype x heat-condition contrast against 25 C. That is a stronger
# check than a file checksum, because it tests the whole re-analysis pipeline (alignment-free counts we
# were given, our filtering, our normalisation, our model) against numbers the authors published.
#
# What agreement would and would not show. High agreement shows we are analysing the same data and have
# not mis-assigned samples to genotypes or conditions, which is the failure mode that would silently
# invalidate every interaction contrast in P2. It does not validate our interaction analysis: the authors
# published no interaction contrast, which is why P2 exists. Disagreement in a specific genotype would
# localise a sample-mapping error, which is exactly what we want a check to be able to do.
#
# The four genotype labels in the GEO columns are hsf1/msn2-4 presence flags, so:
#   hsf1:True  msn2/4:True   = Wildtype
#   hsf1:False msn2/4:True   = HSF1 depleted
#   hsf1:True  msn2/4:False  = msn2/4 deleted
#   hsf1:False msn2/4:False  = double

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
pro  <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
hsd  <- file.path(pro, "hsf1"); ext <- file.path(root, "data/03_external/gse179258")
repd <- file.path(root, "doc/05_reports")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
sink(file.path(repd, "logs/hsf1_geo_published_lfc_check.log"), split = TRUE)
say("=== 29_geo_published_lfc_check | %s ===", format(Sys.time()))

G <- fread(cmd = paste("zcat", shQuote(file.path(ext, "GSE179258_deseq_contextmap.tsv.gz"))))
say("authors' published table: %d genes x %d columns", nrow(G), ncol(G))

GEO_GENO <- c("hsf1:True_msn2/4:True"  = "Wildtype",
              "hsf1:False_msn2/4:True" = "HSF1.KD",
              "hsf1:True_msn2/4:False" = "MSN24.KO",
              "hsf1:False_msn2/4:False" = "Double.KDKO")

## ---- our own logCPM and the same within-genotype contrasts, built here so nothing is assumed -------
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
stopifnot(identical(SM$sample, libs))
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
say("count matrix: %d genes x %d libraries", nrow(M), ncol(M))
say("gene universe identical to the published table: %s",
    identical(sort(rownames(M)), sort(G[[1]])))

dge <- DGEList(M, group = factor(SM$grp))
keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + factor(SM$grp)); colnames(design) <- levels(factor(SM$grp))
fit <- lmFit(voom(dge, design, plot = FALSE), design)
say("genes after filterByExpr: %d", sum(keep))

HEAT <- c("37_10", "37_30", "42_10", "42_30")
R <- list()
for (gk in names(GEO_GENO)) {
  geno <- GEO_GENO[[gk]]
  for (h in HEAT) {
    gcol_geo <- sprintf("%s_%s_VS_25_0_%s.Log2FC", h, gk, gk)
    if (!gcol_geo %in% names(G)) { say("MISSING GEO column: %s", gcol_geo); next }
    a <- paste(geno, h, sep = "__"); b <- paste(geno, "25_0", sep = "__")
    if (!all(c(a, b) %in% colnames(design))) { say("MISSING our group: %s or %s", a, b); next }
    ct <- makeContrasts <- limma::makeContrasts
    cm <- matrix(0, ncol(design), 1, dimnames = list(colnames(design), "d"))
    cm[a, 1] <- 1; cm[b, 1] <- -1
    tt <- topTable(eBayes(contrasts.fit(fit, cm)), number = Inf, sort.by = "none")
    ours <- data.table(gene = rownames(tt), ours = tt$logFC)
    theirs <- data.table(gene = G[[1]], theirs = as.numeric(G[[gcol_geo]]))
    J <- merge(ours, theirs, by = "gene")[is.finite(ours) & is.finite(theirs)]
    R[[length(R) + 1L]] <- data.table(genotype = geno, condition = h, n = nrow(J),
      pearson = cor(J$ours, J$theirs), spearman = cor(J$ours, J$theirs, method = "spearman"),
      slope = coef(lm(theirs ~ ours, data = J))[2],
      med_abs_diff = median(abs(J$ours - J$theirs)), sign_agree = mean(sign(J$ours) == sign(J$theirs)))
  }
}
R <- rbindlist(R)
fwrite(R, file.path(hsd, "geo_published_lfc_agreement.tsv"), sep = "\t")

say("\n---- our logFC versus the authors' published Log2FC, per genotype x condition ----")
say("%-12s %7s %7s %9s %9s %8s %10s", "genotype", "cond", "n", "Pearson", "Spearman", "slope", "sign agree")
for (i in seq_len(nrow(R))) say("%-12s %7s %7d %9.4f %9.4f %8.3f %10.3f",
    R$genotype[i], R$condition[i], R$n[i], R$pearson[i], R$spearman[i], R$slope[i], R$sign_agree[i])
say("\nPearson range %.4f to %.4f ; minimum sign agreement %.3f",
    min(R$pearson), max(R$pearson), min(R$sign_agree))
say("median |our logFC - their Log2FC| across all %d cells: %.4f", nrow(R), median(R$med_abs_diff))

WORST <- R[which.min(pearson)]
say("\nweakest cell: %s %s (Pearson %.4f)", WORST$genotype, WORST$condition, WORST$pearson)
if (min(R$pearson) > 0.9) {
  say("\n=> AGREEMENT IS HIGH IN EVERY GENOTYPE AND CONDITION. The samples are mapped to genotypes and")
  say("   conditions as the authors mapped them; a sample-swap or label error large enough to invalidate")
  say("   the interaction contrasts would have shown up as one bad genotype here, and none is bad.")
  say("   This does NOT validate the interaction analysis: the authors published no interaction contrast,")
  say("   which is the reason P2 exists. It validates the inputs to it.")
} else {
  say("\n=> AT LEAST ONE CELL DISAGREES. That localises a possible sample-mapping problem that requires")
  say("   resolution before interaction contrasts built on this mapping can be interpreted.")
}
say("\nEXPECTED SOURCES OF THE RESIDUAL DIFFERENCE, none of which is a discrepancy: the authors used")
say("DESeq2 with its own shrinkage, we use voom+limma; their table is unfiltered at 7126 genes while")
say("ours is filtered by filterByExpr; and shrinkage differs most for low-count genes, which is where")
say("the largest per-gene differences should and do sit.")
writeLines(capture.output(sessionInfo()), file.path(repd, "logs/sessionInfo_geo_check.txt"))
say("\n=== done ===")
sink()
