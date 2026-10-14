#!/usr/bin/env Rscript
## 104_p2_class_split.R -- selection check for the P2 epistasis classes: classes are re-defined on two replicates and
## the class interaction I(37/30) minus the not-induced reference is measured on the held-out third replicate, beside
## the in-sample estimate and the overlap with the published classes. The criterion (held-out sign equals the published
## sign in every fold, and the held-out mean is at least half the in-sample mean) was fixed before this analysis was run.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
E_all <- cpm(dge, log = TRUE, prior.count = 0.5)   # voom's log-CPM scale
G4 <- c("Double.KDKO", "HSF1.KD", "MSN24.KO", "Wildtype"); SG <- c(1, -1, -1, 1); A <- 0.05
classify <- function(reps) {
  s <- SM$replicate %in% reps; d <- dge[, s]; g <- factor(SM$grp[s])
  des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
  fit <- lmFit(voom(d, des, plot = FALSE), des)
  ct <- function(a, b) { v <- setNames(rep(0, ncol(des)), colnames(des)); v[a] <- 1; v[b] <- -1; v }
  C <- cbind(wt = ct("Wildtype__37_30", "Wildtype__25_0"), kd = ct("HSF1.KD__37_30", "Wildtype__37_30"),
             ko = ct("MSN24.KO__37_30", "Wildtype__37_30"), db = ct("Double.KDKO__37_30", "Wildtype__37_30"))
  iv <- setNames(rep(0, ncol(des)), colnames(des)); iv[paste(G4, "37_30", sep = "__")] <- SG
  f2 <- eBayes(contrasts.fit(fit, cbind(C, I = iv)))
  r <- lapply(c("wt", "kd", "ko", "db"), function(cn) topTable(f2, coef = cn, number = Inf, sort.by = "none")[, c("logFC", "adj.P.Val")])
  names(r) <- c("wt", "kd", "ko", "db")
  ind <- r$wt$adj.P.Val < A & r$wt$logFC > 0
  req <- function(x) ind & x$adj.P.Val < A & x$logFC < 0
  h <- req(r$kd); m <- req(r$ko); b <- req(r$db)
  cls <- ifelse(!ind, "not_induced", ifelse(h & !m, "Hsf1_only", ifelse(!h & m, "Msn24_only", ifelse(h & m, "both_required",
          ifelse(b, "redundant", "TF_independent")))))
  data.table(gene = rownames(f2), cls = cls, I_in = f2$coefficients[, "I"])
}
I_rep <- function(k) { cols <- vapply(G4, function(g) which(SM$knockout == g & SM$cond == "37_30" & SM$replicate == k), 0L)
  as.vector(E_all[, cols] %*% SG) }
PUB <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
PU <- PUB[, .(k = uniqueN(class), nm = uniqueN(method), pcls = class[1]), by = gene][, pcls := ifelse(k == 1L & nm == 3L, pcls, "unassigned")]
OUT <- rbindlist(lapply(1:3, function(k) {
  C <- classify(setdiff(1:3, k)); C[, I_out := I_rep(k)[match(gene, rownames(E_all))]]
  ref_in <- C[cls == "not_induced", median(I_in)]; ref_out <- C[cls == "not_induced", median(I_out)]
  ov <- merge(C, PU[, .(gene, pcls)], by = "gene")
  C[, .(n = .N, in_sample = median(I_in) - ref_in, held_out = median(I_out) - ref_out,
        overlap_with_published = mean(ov[cls == .BY$cls, pcls] == .BY$cls)), by = cls][, held_out_rep := k] }))
fwrite(OUT, file.path(hsd, "p2_class_split.tsv"), sep = "\t")
print(OUT[order(cls, held_out_rep), .(cls, rep = held_out_rep, n, in_sample = round(in_sample, 3), held_out = round(held_out, 3), overlap = round(overlap_with_published, 2))])
S <- OUT[, .(n_mean = mean(n), in_mean = mean(in_sample), out_mean = mean(held_out), out_min = min(held_out), out_max = max(held_out)), by = cls]
print(S[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
pub <- c(Msn24_only = -1.517, redundant = -1.788, Hsf1_only = -0.336)
for (cl in names(pub)) { o <- OUT[cls == cl]; s <- S[cls == cl]
  ok <- all(sign(o$held_out) == sign(pub[[cl]])) && abs(s$out_mean) >= 0.5 * abs(s$in_mean)
  cat(sprintf("RESULT %s: %s (held-out mean %.3f, in-sample %.3f)\n", cl, if (ok) "not a selection artefact" else "not shown to survive held-out estimation", s$out_mean, s$in_mean)) }
