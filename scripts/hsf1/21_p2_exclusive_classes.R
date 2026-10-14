# 21_p2_exclusive_classes.R -- method-unanimous, mutually exclusive class assignment for the P2
# class table. The class file carries one assignment per gene per DE method and the three methods
# disagree for 312 genes; taking the union of gene-class pairs double-counts them and contaminates
# the reference class. A gene is assigned only where all three methods agree; the rest are labelled
# "unassigned" and reported as their own row. Exclusivity is checked with stopifnot().
suppressPackageStartupMessages(library(data.table))
CLS <- fread("data/02_processed/hsf1/hsf1_epistasis_classes_37_30.tsv")
gcol <- names(CLS)[1]; cc <- intersect(c("class","epistasis_class","cls"), names(CLS))[1]
mcol <- intersect(c("method","meth"), names(CLS))[1]
cat("gene col:", gcol, "| class col:", cc, "| method col:", mcol, "\n")
stopifnot(!is.na(mcol))
CL <- unique(CLS[, .(gene=get(gcol), cls=get(cc), method=get(mcol))])
cat("methods:", paste(sort(unique(CL$method)), collapse=", "), "\n")
# unanimous across the three methods, else 'unassigned'
U <- CL[, .(k=uniqueN(cls), n_meth=uniqueN(method), cls=cls[1]), by=gene]
U[, final := ifelse(k==1L & n_meth==uniqueN(CL$method), cls, "unassigned")]
stopifnot(!any(duplicated(U$gene)))          # one class per gene, checked
cat("\nunanimous:", U[final!="unassigned", .N], "| non-unanimous -> unassigned:", U[final=="unassigned", .N], "\n")
print(U[, .N, by=final][order(-N)])
IX <- fread(cmd="zcat data/02_processed/hsf1/p2_offset_origin_genewise.tsv.gz")
M <- merge(IX, U[, .(gene, cls=final)], by="gene")
cat("\nmerged rows per regime:", M[, .N, by=regime]$N[1], "(must equal model genes with a label)\n")
R <- M[, .(n=.N, med_I=median(I), se_of_med=median(se)/sqrt(.N)), by=.(regime, cls)]
ref <- R[cls=="not_induced", .(regime, refv=med_I)]
R <- merge(R, ref, by="regime")[, `:=`(z0=med_I/se_of_med, vsref=med_I-refv, zref=(med_I-refv)/se_of_med)]
fwrite(R, "data/02_processed/hsf1/p2_class_exclusive.tsv", sep="\t")
for (rg in unique(R$regime)) {
  cat("\n=== ", rg, " ===\n", sep="")
  z <- R[regime==rg][order(-n)]
  cat(sprintf("%-16s %6s %9s %8s %10s %9s\n","class","n","med_I","z vs 0","vs ref","z vs ref"))
  for (i in seq_len(nrow(z))) cat(sprintf("%-16s %6d %+9.4f %8.1f %+10.4f %9.1f\n",
    z$cls[i], z$n[i], z$med_I[i], z$z0[i], z$vsref[i], z$zref[i]))
}
