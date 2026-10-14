#!/usr/bin/env Rscript
## 76_p2_i25_resolvability.R -- is the response-scale interaction I(c) - I(25) resolvable, and why does its
## amplitude at 37 °C/10 min exceed both parts? Read-only on the voom matrix of script 18; no new model.
##
## (1) Per replicate pairing (3 pairings), the noise-corrected amplitude rms_tru, rms_null and SNR of I(c), I(25) and
##     I(c) - I(25), and the number of pairings at which rms_tru collapses to 0: the two quantities P2 reports.
## (2) Decomposition. With S = signal contrast (a+b)/2 and N = null contrast (b-a)/2 of one pairing,
##       A_d^2 = A_c^2 + A_25^2 - 2 (X_S - X_N),  X_S = mean(S_c S_25), X_N = mean(N_c N_25),
##     so X_S - X_N is the pairing-based estimate of the true cross moment, and cor_true = (X_S - X_N)/(A_c A_25).
##     This estimate relies on replicate r at condition c being paired with replicate r at 25 °C. The deposited records
##     do not say whether replicate r at 25 °C and at 37 °C come from one culture, so (3) and (4) check the pairing.
## (3) A pairing-free estimate: the 3-replicate group means of c and of 25 °C have independent noise if the cultures
##     are independent, so mean(Ibar_c Ibar_25) estimates the true cross moment without any pairing; A_d is rebuilt
##     from it with the same A_c, A_25.
## (4) Relabelling: all 6 assignments of the three 25 °C replicate labels to r1..r3, and the resulting A_d.
## (5) Gene-level correlation of the model coefficients I(c) and I(25) (script 72's fit, published normalisation),
##     which includes noise and is reported as such.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
say <- function(...) cat(sprintf(...), "\n", sep = "")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + factor(SM$grp)); colnames(design) <- levels(factor(SM$grp))
E <- voom(dge, design, plot = FALSE)$E
G4 <- c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"); SG <- c(1, -1, -1, 1)
col_of <- function(g, c, r) which(SM$knockout == g & SM$cond == c & SM$replicate == r)
I_rep <- function(c, r, relab = 1:3) { rr <- if (c == "25_0") relab[r] else r
  i <- vapply(G4, function(g) col_of(g, c, rr), 0L); as.vector(E[, i] %*% SG) }
rms_true <- function(s, n) sqrt(max(0, mean(s^2) - mean(n^2)))
HEAT <- c("37_10", "37_30", "42_10", "42_30"); PR <- combn(1:3, 2)
LA <- fread(file.path(hsd, "p2_i25_amplitude.tsv"))
per <- list(); dec <- list()
for (c_ in HEAT) for (j in 1:3) {
  a <- PR[1, j]; b <- PR[2, j]
  Sc <- (I_rep(c_, a) + I_rep(c_, b)) / 2; Nc <- (I_rep(c_, b) - I_rep(c_, a)) / 2
  S0 <- (I_rep("25_0", a) + I_rep("25_0", b)) / 2; N0 <- (I_rep("25_0", b) - I_rep("25_0", a)) / 2
  Sd <- Sc - S0; Nd <- Nc - N0
  per[[length(per) + 1]] <- data.table(condition = c_, pairing = sprintf("r%d-r%d", a, b),
    A_c = rms_true(Sc, Nc), N_c = sqrt(mean(Nc^2)), A_25 = rms_true(S0, N0), N_25 = sqrt(mean(N0^2)),
    A_d = rms_true(Sd, Nd), N_d = sqrt(mean(Nd^2)), X_S = mean(Sc * S0), X_N = mean(Nc * N0))
}
P <- rbindlist(per)
P[, `:=`(snr_d = A_d / N_d, cross_true = X_S - X_N)]
## fidelity against script 72's stored means
fid <- P[, .(A_d = mean(A_d)), by = condition][, q := sprintf("I(%s) - I(25)", condition)]
stopifnot(max(abs(fid$A_d - LA[match(fid$q, quantity), rms_tru])) < 1e-10)
say("fidelity: the four saved baseline-subtracted amplitudes are reproduced")
## pairing-free cross moment from group means
Ibar <- function(c_) rowMeans(sapply(1:3, function(r) I_rep(c_, r)))
PF <- rbindlist(lapply(HEAT, function(c_) {
  X <- mean(Ibar(c_) * Ibar("25_0")); Ac <- LA[quantity == sprintf("I(%s)", c_), rms_tru]; A0 <- LA[quantity == "I(25)", rms_tru]
  data.table(condition = c_, cross_groupmeans = X, A_c = Ac, A_25 = A0,
             A_d_pairing_free = sqrt(max(0, Ac^2 + A0^2 - 2 * X)), cor_true_pairing_free = X / (Ac * A0))
}))
## relabelling of the 25 °C replicates
perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
RL <- rbindlist(lapply(HEAT, function(c_) rbindlist(lapply(seq_along(perms), function(k) {
  pk <- perms[[k]]
  v <- vapply(1:3, function(j) { a <- PR[1, j]; b <- PR[2, j]
    Sd <- (I_rep(c_, a) + I_rep(c_, b)) / 2 - (I_rep("25_0", a, pk) + I_rep("25_0", b, pk)) / 2
    Nd <- (I_rep(c_, b) - I_rep(c_, a)) / 2 - (I_rep("25_0", b, pk) - I_rep("25_0", a, pk)) / 2
    c(rms_true(Sd, Nd), sqrt(mean(Nd^2))) }, numeric(2))
  data.table(condition = c_, labels_25 = paste(pk, collapse = ""), A_d = mean(v[1, ]), N_d = mean(v[2, ]),
             zero_pairings = sum(v[1, ] == 0))
}))))
RL[, snr_d := A_d / N_d]
## model coefficients, published normalisation (script 72's fit)
fit <- lmFit(voom(dge, design, plot = FALSE), design)
mk <- function(c_) { ct <- setNames(rep(0, ncol(design)), colnames(design)); ct[paste(G4, c_, sep = "__")] <- SG; ct }
cf <- contrasts.fit(fit, sapply(c("25_0", HEAT), mk))$coefficients
CORM <- data.table(condition = HEAT, cor_model_coef = vapply(HEAT, function(c_) cor(cf[, c_], cf[, "25_0"]), 0))
SUM <- merge(merge(P[, .(A_c = mean(A_c), A_25 = mean(A_25), A_d = mean(A_d), N_d = mean(N_d), snr_d = mean(A_d) / mean(N_d),
                         zero_pairings = sum(A_d == 0), cross_pairing = mean(cross_true),
                         cor_true_pairing = mean(cross_true) / (mean(A_c) * mean(A_25))), by = condition],
                   PF[, .(condition, cross_groupmeans, A_d_pairing_free, cor_true_pairing_free)], by = "condition"),
             merge(CORM, RL[, .(A_d_relabel_min = min(A_d), A_d_relabel_max = max(A_d), snr_relabel_min = min(snr_d),
                                snr_relabel_max = max(snr_d), zero_pairings_relabel_max = max(zero_pairings)), by = condition],
                   by = "condition"), by = "condition")
fwrite(P, file.path(hsd, "p2_i25_resolvability_bypair.tsv"), sep = "\t")
fwrite(RL, file.path(hsd, "p2_i25_relabel.tsv"), sep = "\t")
fwrite(SUM, file.path(hsd, "p2_i25_resolvability.tsv"), sep = "\t")
say("=== 76_p2_i25_resolvability | %s ===", format(Sys.time()))
print(P[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
print(t(SUM[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)]))
print(RL[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
