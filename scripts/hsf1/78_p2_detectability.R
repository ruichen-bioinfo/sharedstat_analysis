#!/usr/bin/env Rscript
## 78_p2_detectability.R -- P2 detectability, with the decision criteria fixed before this analysis was run.
## Three analyses are performed:
##   Presence    rms_tru with a 95% gene-resampling bootstrap interval (2000 resamples) for every admissible pairing:
##                  level scale I(c), 3 replicate pairings per condition; response scale I(c) - I(25), the 3 pairings
##                  under each of the 6 relabellings of the 25 °C replicates (18). Present if every lower bound > 0 and
##                  no pairing collapses to 0; otherwise not detected at this replication.
##   Magnitude   rms_tru, and the 37 °C / 42 °C ratio at each timepoint, same resamples, same pairing; the contrast
##                  holds if every pairing's interval excludes 1.
##   Readability h = SNR^2 / (1 + SNR^2) from the mean-over-pairings amplitude and floor; readable if h >= 0.5.
## Two post hoc diagnostics are reported separately and do not affect the analysis:
##   Replicate-adjusted h for the three-replicate estimate the paper reports, by Spearman-Brown: the planned SNR compares the signal
##      with the noise of a two-replicate average (null = (b - a)/2), whereas the per-gene I(c) of the model averages
##      three replicates, so h3 = SNR^2 / (SNR^2 + 2/3).
##   Inter-gene correlation of the model residuals per class (limma::interGeneCorrelation, as in camera, Wu and
##      Smyth 2012), the variance inflation 1 + (n - 1) rho of a class-median-type statistic, and the class z values of
##      Table 1 divided by its square root. The presence bootstrap and the SE of a class median both treat genes as
##      independent; the residual-correlation diagnostic measures how far that fails.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
say <- function(...) cat(sprintf(...), "\n", sep = "")
NB <- 2000L; set.seed(20260926L)
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + factor(SM$grp)); colnames(design) <- levels(factor(SM$grp))
v <- voom(dge, design, plot = FALSE); E <- v$E
G4 <- c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"); SG <- c(1, -1, -1, 1)
col_of <- function(g, c, r) which(SM$knockout == g & SM$cond == c & SM$replicate == r)
I_rep <- function(c, r) as.vector(E[, vapply(G4, function(g) col_of(g, c, r), 0L)] %*% SG)
HEAT <- c("37_10", "37_30", "42_10", "42_30"); PR <- combn(1:3, 2)
perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
IR <- lapply(c("25_0", HEAT), function(c_) sapply(1:3, function(r) I_rep(c_, r))); names(IR) <- c("25_0", HEAT)
## per-pairing per-gene s^2 - n^2 and n^2 (the bootstrap works on these columns)
cols <- list(); meta_ <- list()
add <- function(scale, cond, lab, j, s, n) { k <- length(cols) + 1
  cols[[k]] <<- cbind(d = s^2 - n^2, n2 = n^2)
  meta_[[k]] <<- data.table(scale = scale, condition = cond, labels_25 = lab, pairing = sprintf("r%d-r%d", PR[1, j], PR[2, j])) }
for (c_ in c("25_0", HEAT)) for (j in 1:3) { a <- PR[1, j]; b <- PR[2, j]; X <- IR[[c_]]
  add("level", c_, "-", j, (X[, a] + X[, b]) / 2, (X[, b] - X[, a]) / 2) }
for (c_ in HEAT) for (pk in perms) for (j in 1:3) { a <- PR[1, j]; b <- PR[2, j]; X <- IR[[c_]]; Z <- IR[["25_0"]][, pk]
  add("response", c_, paste(pk, collapse = ""), j, (X[, a] + X[, b]) / 2 - (Z[, a] + Z[, b]) / 2,
      (X[, b] - X[, a]) / 2 - (Z[, b] - Z[, a]) / 2) }
MT <- rbindlist(meta_); K <- length(cols); G <- nrow(E)
Dm <- sapply(cols, function(x) x[, "d"]); Nm <- sapply(cols, function(x) x[, "n2"])
amp <- function(w) sqrt(pmax(0, colSums(Dm * w) / sum(w)))           # w = resample multiplicities
A0 <- amp(rep(1, G))
BOOT <- matrix(NA_real_, NB, K)
for (b in seq_len(NB)) BOOT[b, ] <- amp(tabulate(sample.int(G, G, replace = TRUE), G))
MT[, `:=`(A = A0, lo = apply(BOOT, 2, quantile, 0.025), hi = apply(BOOT, 2, quantile, 0.975),
          N = sqrt(colSums(Nm) / G))]
## fidelity: the deposited-label mean-over-pairings amplitudes equal the stored values
LA <- fread(file.path(hsd, "p2_i25_amplitude.tsv"))
fl <- MT[scale == "level" & condition %in% HEAT, .(A = mean(A)), by = condition]
fr <- MT[scale == "response" & labels_25 == "123", .(A = mean(A)), by = condition]
stopifnot(max(abs(fl$A - LA[match(sprintf("I(%s)", fl$condition), quantity), rms_tru])) < 1e-10,
          max(abs(fr$A - LA[match(sprintf("I(%s) - I(25)", fr$condition), quantity), rms_tru])) < 1e-10,
          abs(MT[scale == "level" & condition == "25_0", mean(A)] - LA[quantity == "I(25)", rms_tru]) < 1e-10)
say("fidelity: saved amplitudes reproduced (level 5 conditions, response 4 conditions)")
## presence analysis
Q1 <- MT[, .(n_pairings = .N, min_A = min(A), min_lo = min(lo), zero_pairings = sum(A == 0),
             present = all(lo > 0) & all(A > 0)), by = .(scale, condition)]
## magnitude analysis: ratio 37/42 at each timepoint, per matching pairing (same labels, same pairing), same resamples
Q2 <- rbindlist(lapply(c("level", "response"), function(sc) rbindlist(lapply(c("10", "30"), function(tp) {
  a <- MT[, which(scale == sc & condition == paste0("37_", tp))]; b <- MT[, which(scale == sc & condition == paste0("42_", tp))]
  stopifnot(length(a) == length(b), all(MT$pairing[a] == MT$pairing[b]), all(MT$labels_25[a] == MT$labels_25[b]))
  rbindlist(lapply(seq_along(a), function(i) { r <- BOOT[, a[i]] / BOOT[, b[i]]
    r[BOOT[, a[i]] == 0 & BOOT[, b[i]] == 0] <- 1          # both collapsed: no contrast; x/0 stays Inf (37 °C > 0)
    data.table(scale = sc, timepoint = tp, labels_25 = MT$labels_25[a[i]], pairing = MT$pairing[a[i]],
               ratio = if (A0[a[i]] == 0 && A0[b[i]] == 0) 1 else A0[a[i]] / A0[b[i]], lo = quantile(r, 0.025), hi = quantile(r, 0.975)) })) }))))
Q2s <- Q2[, .(n_pairings = .N, min_lo = min(lo), ratio_median = median(ratio), holds = all(lo > 1)), by = .(scale, timepoint)]
## readability and replicate-adjusted reliability
Q3 <- MT[labels_25 %in% c("-", "123"), .(A = mean(A), N = mean(N)), by = .(scale, condition)]
Q3[, snr := A / N][, `:=`(h2 = snr^2 / (1 + snr^2), h3 = snr^2 / (snr^2 + 2/3))][, readable := h2 >= 0.5]
## residual-correlation diagnostic: inter-gene correlation per exclusive class, from the residuals of the group-means model
CLS3 <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
U <- CLS3[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene][, cls := ifelse(k == 1L & nm == 3L, cls, "unassigned")]
T1 <- fread(file.path(hsd, "p2_class_exclusive.tsv"))[regime == "all genes (as published)"]
D2 <- rbindlist(lapply(unique(T1$cls), function(cl) { g <- intersect(U[cls == cl, gene], rownames(v))
  if (length(g) > 2000) { set.seed(1); g <- sample(g, 2000) }                      # correlation estimated on <= 2000 genes
  ic <- interGeneCorrelation(v$E[g, ], design)$correlation
  n <- T1[cls == cl, n]; vif <- 1 + (n - 1) * ic
  data.table(cls = cl, n = n, inter_gene_cor = ic, vif = vif, zref = T1[cls == cl, zref], zref_adj = T1[cls == cl, zref] / sqrt(max(1, vif))) }))
fwrite(MT, file.path(hsd, "p2_detectability_bypairing.tsv"), sep = "\t")
fwrite(Q1, file.path(hsd, "p2_detectability_presence.tsv"), sep = "\t")
fwrite(Q2, file.path(hsd, "p2_detectability_ratio_bypairing.tsv"), sep = "\t")
fwrite(Q2s, file.path(hsd, "p2_detectability_ratio.tsv"), sep = "\t")
fwrite(Q3, file.path(hsd, "p2_detectability_readability.tsv"), sep = "\t")
fwrite(D2, file.path(hsd, "p2_intergene_correlation.tsv"), sep = "\t")
say("=== 78_p2_detectability | %s | %d pairings x %d resamples ===", format(Sys.time()), K, NB)
say("\nPRESENCE"); print(Q1[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
say("\nMAGNITUDE: 37/42 ratio"); print(Q2s[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
say("\nREADABILITY (h2 = criterion; h3 = replicate-adjusted diagnostic)"); print(Q3[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
say("\ninter-gene correlation and class z diagnostic"); print(D2[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
