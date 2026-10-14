#!/usr/bin/env Rscript
## 104b_p2_class_split.R -- held-out class split, fixed before this analysis was run:
## no leak (filter + TMM inside each fold), selection-bias contrast.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT"); hsd <- file.path(root, "data/02_processed/hsf1"); pro <- file.path(root, "data/02_processed")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(root, "data/01_metadata/sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
G4 <- c("Double.KDKO", "HSF1.KD", "MSN24.KO", "Wildtype"); SG <- c(1, -1, -1, 1); A <- 0.05
classify <- function(reps) {
  s <- SM$replicate %in% reps; g <- factor(SM$grp[s]); d <- DGEList(M[, s], group = g)
  d <- calcNormFactors(d[filterByExpr(d, group = g), , keep.lib.sizes = FALSE])
  des <- model.matrix(~ 0 + g); colnames(des) <- levels(g); fit <- lmFit(voom(d, des, plot = FALSE), des)
  ct <- function(a, b) { v <- setNames(rep(0, ncol(des)), colnames(des)); v[a] <- 1; v[b] <- -1; v }
  C <- cbind(wt = ct("Wildtype__37_30", "Wildtype__25_0"), kd = ct("HSF1.KD__37_30", "Wildtype__37_30"),
             ko = ct("MSN24.KO__37_30", "Wildtype__37_30"), db = ct("Double.KDKO__37_30", "Wildtype__37_30"))
  f2 <- eBayes(contrasts.fit(fit, C)); r <- lapply(colnames(C), function(cn) topTable(f2, coef = cn, number = Inf, sort.by = "none")); names(r) <- colnames(C)
  ind <- r$wt$adj.P.Val < A & r$wt$logFC > 0; req <- function(x) ind & x$adj.P.Val < A & x$logFC < 0
  h <- req(r$kd); m <- req(r$ko); b <- req(r$db)
  data.table(gene = rownames(f2), cls = ifelse(!ind, "not_induced", ifelse(h & !m, "Hsf1_only", ifelse(!h & m, "Msn24_only",
    ifelse(h & m, "both_required", ifelse(b, "redundant", "TF_independent")))))) }
## held-out replicate k: its own filter/TMM on the 4 libraries of 37/30 plus 25 libs of rep k are not enough, so the
## held-out log-CPM uses TMM on all 20 libraries of replicate k only (no information from the defining replicates).
I_rep <- function(k) { s <- SM$replicate == k; d <- calcNormFactors(DGEList(M[, s])); E <- cpm(d, log = TRUE, prior.count = 0.5)
  ss <- SM[s]; cols <- vapply(G4, function(g) which(ss$knockout == g & ss$cond == "37_30"), 0L); setNames(as.vector(E[, cols] %*% SG), rownames(E)) }
med_diff <- function(C, I) { C <- C[gene %in% names(I)]; C[, I := I[gene]]; r <- C[cls == "not_induced", median(I)]; C[, .(n = .N, d = median(I) - r), by = cls] }
CALL <- classify(1:3)
OUT <- rbindlist(lapply(1:3, function(k) { Ik <- I_rep(k)
  Ck <- classify(setdiff(1:3, k)); a <- med_diff(Ck, Ik)[, .(cls, n_out = n, held_out = d)]
  Itr <- (I_rep(setdiff(1:3, k)[1]) + I_rep(setdiff(1:3, k)[2])) / 2       # in-sample estimate on the defining replicates
  a <- merge(a, med_diff(Ck, Itr)[, .(cls, in_sample = d)], by = "cls")
  b <- med_diff(CALL, Ik)[, .(cls, n_all = n, with_k = d)]
  merge(a, b, by = "cls", all = TRUE)[, rep := k] }))
OUT[, selection_bias := with_k - held_out]
fwrite(OUT, file.path(hsd, "p2_class_split_b.tsv"), sep = "\t")
S <- OUT[, .(held_out_mean = mean(held_out), lo_t2 = mean(held_out) - qt(.975, 2) * sd(held_out) / sqrt(3), hi_t2 = mean(held_out) + qt(.975, 2) * sd(held_out) / sqrt(3),
             with_k_mean = mean(with_k), in_sample_mean = mean(in_sample), selection_bias = mean(selection_bias)), by = cls]
print(OUT[order(cls, rep), lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); print(S[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
pub <- c(Msn24_only = -1.517, redundant = -1.788, Hsf1_only = -0.336, TF_independent = -0.662)
for (cl in names(pub)) { o <- OUT[cls == cl]; s <- S[cls == cl]
  ok <- all(sign(o$held_out) == sign(pub[[cl]])) && abs(s$held_out_mean) >= 0.5 * abs(s$in_sample_mean)   # as in script 104, the comparator is the in-sample estimate on the defining replicates
  cat(sprintf("RESULT %s: %s\n", cl, if (ok) "not a selection artefact" else "not shown to survive")) }
