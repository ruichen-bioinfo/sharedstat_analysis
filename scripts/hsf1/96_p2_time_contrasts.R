# 96_p2_time_contrasts.R -- P2 time contrasts, with the decision criteria fixed before this analysis was run.
# Gene-level time contrasts in the script-17 group-means model; limma treat and BH/|LFC| double
# threshold at t = 1 and t = log2(1.5); peak classes; class medians with correlation-adjusted SE.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); pro <- file.path(root, "data/02_processed"); hsd <- file.path(pro, "hsf1")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(root, "data/01_metadata/sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")]; SM[, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
stopifnot(identical(SM$sample, libs), all(SM[, .N, by = grp]$N == 3L))
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(counts = M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge)
grp <- factor(SM$grp); des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
v <- voom(dge, des, plot = FALSE); fit <- lmFit(v, des)
cat("genes", nrow(v), "(must be 6070)\n"); stopifnot(nrow(v) == 6070)
G <- c("Wildtype", "HSF1.KD", "MSN24.KO", "Double.KDKO")
e <- function(...) { a <- list(...); ct <- setNames(rep(0, ncol(des)), colnames(des)); for (i in seq(1, length(a), 2)) ct[a[[i]]] <- ct[a[[i]]] + a[[i + 1]]; ct }
n_ <- function(g, c) paste(g, c, sep = "__")
CT <- list()
for (g in G) {
  CT[[paste0("E37|", g)]] <- e(n_(g, "37_10"), 1, n_(g, "25_0"), -1)
  CT[[paste0("L37|", g)]] <- e(n_(g, "37_30"), 1, n_(g, "37_10"), -1)
  CT[[paste0("E42|", g)]] <- e(n_(g, "42_10"), 1, n_(g, "25_0"), -1)
  CT[[paste0("L42|", g)]] <- e(n_(g, "42_30"), 1, n_(g, "42_10"), -1)
  CT[[paste0("F37|", g)]] <- e(n_(g, "37_30"), 1, n_(g, "25_0"), -1)
  CT[[paste0("F42|", g)]] <- e(n_(g, "42_30"), 1, n_(g, "25_0"), -1)
}
for (T in c("37", "42")) for (g in c("HSF1.KD", "MSN24.KO"))
  CT[[paste0("TxG", T, "|", g)]] <- CT[[paste0("L", T, "|", g)]] - CT[[paste0("L", T, "|Wildtype")]]
CTm <- do.call(cbind, CT); cf <- contrasts.fit(fit, CTm); eb <- eBayes(cf)
TS <- c(`1` = 1, `0.585` = log2(1.5))
R <- rbindlist(lapply(colnames(CTm), function(k) rbindlist(lapply(names(TS), function(tn) {
  tr <- treat(cf, fc = 2^TS[[tn]]); j <- k
  data.table(gene = rownames(cf), contrast = k, t = tn, lfc = cf$coefficients[, j],
             se = cf$stdev.unscaled[, j] * cf$sigma,
             padj_treat = p.adjust(tr$p.value[, j], "BH"), padj_eb = p.adjust(eb$p.value[, j], "BH"))
}))))
R[, `:=`(sig_treat = padj_treat < 0.05, sig_dbl = padj_eb < 0.05 & abs(lfc) >= TS[t])]
fwrite(R, file.path(hsd, "p2_time_contrasts_genewise.tsv.gz"), sep = "\t", compress = "gzip")
# classes
CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))[, .(gene, cls = class, method)]
U <- unique(CL)[, .(k = uniqueN(cls), m = uniqueN(method), cls = cls[1]), by = gene]
U[, cls := ifelse(k == 1L & m == 3L, cls, "unassigned")]
R <- merge(R, U[, .(gene, cls)], by = "gene")
CLS <- c("Msn24_only", "redundant", "Hsf1_only", "TF_independent", "both_required", "not_induced")
# class medians with VIF-adjusted SE (script 80 method)
rho_of <- function(g) { g <- intersect(g, rownames(v)); if (length(g) > 2000) { set.seed(1); g <- sample(g, 2000) }; interGeneCorrelation(v$E[g, ], des)$correlation }
RH <- data.table(cls = CLS)[, rho := sapply(cls, function(c) rho_of(U[cls == c, gene]))]
S <- R[t == "1" & cls %in% CLS, .(n = .N, med = median(lfc), se_gene = median(se) / sqrt(.N)), by = .(cls, contrast)]
S <- merge(S, RH, by = "cls")[, vif := 1 + (n - 1) * pmax(rho, 0)][, se_adj := sqrt(vif) * se_gene]
S[, `:=`(lo = med - 1.96 * se_adj, hi = med + 1.96 * se_adj)]
fwrite(S, file.path(hsd, "p2_time_class_medians.tsv"), sep = "\t")
# significance counts per class, both tests, both thresholds, with noise floor
C <- R[cls %in% CLS, .(n = .N, down_treat = sum(sig_treat & lfc < 0), up_treat = sum(sig_treat & lfc > 0),
                       down_dbl = sum(sig_dbl & lfc < 0), up_dbl = sum(sig_dbl & lfc > 0), noise_floor = median(se)), by = .(cls, contrast, t)]
fwrite(C, file.path(hsd, "p2_time_sig_counts.tsv"), sep = "\t")
# peak classes
PK <- rbindlist(lapply(c("37", "42"), function(T) rbindlist(lapply(G, function(g) {
  E <- R[contrast == paste0("E", T, "|", g), .(gene, t, cls, e_sig = sig_treat & lfc > 0, e_sd = sig_dbl & lfc > 0)]
  L <- R[contrast == paste0("L", T, "|", g), .(gene, t, l_sig = sig_treat, l_lfc = lfc, l_sd = sig_dbl)]
  X <- merge(E, L, by = c("gene", "t"))
  rbind(X[e_sig == TRUE, .(test = "treat", n_induced = .N, transient = sum(l_sig & l_lfc < 0), sustained = sum(!l_sig), rising = sum(l_sig & l_lfc > 0)), by = .(cls, t)],
        X[e_sd == TRUE, .(test = "BH+LFC", n_induced = .N, transient = sum(l_sd & l_lfc < 0), sustained = sum(!l_sd), rising = sum(l_sd & l_lfc > 0)), by = .(cls, t)])[, `:=`(temp = T, genotype = g)]
}))))
PK[, frac_transient := round(transient / n_induced, 3)]
fwrite(PK, file.path(hsd, "p2_time_peak_classes.tsv"), sep = "\t")
# report
p <- function(x) cat(x, "\n")
m <- S[cls == "Msn24_only"]; h <- S[cls == "Hsf1_only"]; f <- function(D, k) D[contrast == k]
for (k in c("L37|Wildtype", "L42|Wildtype", "L37|HSF1.KD", "L42|HSF1.KD", "TxG37|HSF1.KD", "TxG42|HSF1.KD")) { z <- f(m, k); cat(sprintf("Msn24_only %-15s med %+.3f [%+.3f, %+.3f] n %d vif %.1f\n", k, z$med, z$lo, z$hi, z$n, z$vif)) }
for (k in c("TxG37|HSF1.KD", "TxG42|HSF1.KD", "L37|Wildtype")) { z <- f(h, k); cat(sprintf("Hsf1_only  %-15s med %+.3f [%+.3f, %+.3f]\n", k, z$med, z$lo, z$hi)) }
q1 <- f(m, "L37|Wildtype")$hi < 0; q2 <- !(f(m, "L42|Wildtype")$hi < 0); q3 <- f(m, "TxG37|HSF1.KD")$lo > 0
pk <- PK[cls == "Msn24_only" & temp == "37" & genotype %in% c("Wildtype", "HSF1.KD")]
q4 <- all(sapply(c("1", "0.585"), function(tt) { a <- pk[test == "treat" & t == tt]; a[genotype == "Wildtype", frac_transient] > a[genotype == "HSF1.KD", frac_transient] }))
ctl <- f(h, "TxG37|HSF1.KD")$med < 0
cat(sprintf("\nMsn2/4-only wild-type 37 C decline %s  no wild-type 42 C decline %s  positive HSF1-depletion interaction at 37 C %s  peak-class criterion %s  Hsf1-only control %s\n", q1, q2, q3, q4, ctl))
print(PK[cls %in% c("Msn24_only", "redundant", "Hsf1_only") & genotype %in% c("Wildtype", "HSF1.KD")][order(test, t, temp, cls, genotype)])
