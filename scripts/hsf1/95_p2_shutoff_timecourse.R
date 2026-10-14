# 95_p2_shutoff_timecourse.R -- descriptive check fixed before this analysis was run.
# Class medians of published-pipeline limma log2FC over time; classes exactly as 21_p2_exclusive_classes.R.
suppressPackageStartupMessages(library(data.table))
CLS <- fread("data/02_processed/hsf1/hsf1_epistasis_classes_37_30.tsv")
CL <- unique(CLS[, .(gene, cls = class, method)])
U <- CL[, .(k = uniqueN(cls), n_meth = uniqueN(method), cls = cls[1]), by = gene]
U[, final := ifelse(k == 1L & n_meth == uniqueN(CL$method), cls, "unassigned")]
stopifnot(!any(duplicated(U$gene)))
conds <- c("37_10", "37_30", "42_10", "42_30"); L <- "data/02_processed/limma/"
rd <- function(f) fread(paste0(L, f))[, .(gene, log2FC)]
out <- rbindlist(lapply(conds, function(c) rbindlist(lapply(list(
  WT = paste0("Wildtype_", c, "_vs_25_0.tsv"), KD = paste0("HSF1.KD_", c, "_vs_25_0.tsv"),
  KD_minus_WT = paste0("HSF1.KD_vs_Wildtype_at_", c, ".tsv")), function(f) rd(f)), idcol = "quantity")[, cond := c]))
M <- merge(out, U[, .(gene, cls = final)], by = "gene")
R <- M[cls %in% c("Msn24_only", "redundant", "Hsf1_only", "not_induced"),
       .(n = uniqueN(gene), median_log2FC = round(median(log2FC, na.rm = TRUE), 3)), by = .(cls, quantity, cond)]
fwrite(R, "data/02_processed/hsf1/p2_shutoff_timecourse.tsv", sep = "\t")
W <- dcast(R, cls + quantity ~ cond, value.var = "median_log2FC"); print(W)
m <- W[cls == "Msn24_only"]; g <- function(q, c) m[quantity == q][[c]]
cat("\nP1 WT 37: 30 <", "10 ->", g("WT","37_30") < g("WT","37_10"), "\n")
cat("P2 WT 42: 30 >=", "10 ->", g("WT","42_30") >= g("WT","42_10"), "\n")
k <- sapply(conds, function(c) g("KD_minus_WT", c))
cat("P3 KD-WT largest at 37_30 ->", which.max(k) == 2, " values", k, "\n")
h <- W[cls == "Hsf1_only" & quantity == "KD_minus_WT", ..conds]
cat("Control Hsf1_only KD-WT < 0 at all ->", all(unlist(h) < 0), "\n")
