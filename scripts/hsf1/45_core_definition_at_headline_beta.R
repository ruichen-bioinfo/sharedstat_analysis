#!/usr/bin/env Rscript
## 45_core_definition_at_headline_beta.R -- evaluate the six core definitions at the beta the headline is
## reported at, instead of at beta = 0 only.
##
## Background. The abstract says the negative "holds across six ways of drawing the core boundary, including one
## from Gene Ontology alone, under which the core is resolvably less divergent". That claim rests on
## a5_core_definition_strict.tsv, which has no beta column. Reading script 04, cell_ratio() forms the ratio
## from S_M and N_M with no amplitude term, so it is the beta = 0 statistic; and the strict table's
## reported_core row, 1.1002348 / 1.1612972 / 1.0795665, is digit-for-digit the ratio_beta0 column of
## a5_cis_ratio_calibrated_37_30.tsv. The headline numbers, 1.162 / 1.235 / 1.025, are the beta = 0.383
## column. So the robustness check and the result it defends are computed at different betas, and the
## GO-only claim is not otherwise evaluated at the headline beta.
##
## Method. Transplants script 04's estimator exactly as script 34 does, rather than reimplementing it (the
## reimplementation in script 33 did not match script 04), then generalises the ratio to
##     ratio_beta = rms_tru(DoD | numerator) / rms_tru(DoD | denominator), each divided by amp^beta
## and reruns the same six definitions at beta = 0, at the fitted beta = 0.383444, and at beta = 1.
##
## Two fidelity checks:
##   1. at beta = 0 the generalised function must reproduce a5_core_definition_strict.tsv exactly, so the
##      generalisation is proved to be a generalisation and not a different statistic;
##   2. reported_core at beta = 0.383444 must reproduce the published cal_ratio 1.162484 / 1.234922 / 1.024726.
##
## Decision rule, written before running:
##   If the GO definition still gives a ratio below 1 with a CI excluding 1 at beta = 0.383, the abstract's
##   statement holds with an explicit beta. If it does not, the statement does not hold. Any other
##   pattern is reported as it comes.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1"); csd <- file.path(root, "data/02_processed/crossspecies")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

BETAS_REPORT <- c(0, 0.383444, 1)          # 0.383444 is the fitted beta in a5_cis_ratio_calibrated
NB <- 2000L

## ------------------------------------------------------------------ transplant, as in script 34 ----
SRC <- file.path(root, "scripts/hsf1/04_divergence_by_dependency.R")
EX  <- parse(SRC)
need <- c("cell_ratio", "cell_tru", "cell_amp", "idx_of", "HOGS", "DENOM", "TP", "tp_of", "CELLS",
          "S_M", "N_M", "R_M", "RN_M")
have <- function() all(vapply(need, exists, TRUE, envir = globalenv()))
n_eval <- 0L
for (e in EX) {
  s <- paste(deparse(e), collapse = " ")
  if (grepl("^sink\\(|^fwrite\\(|quit\\(", s)) next
  if (grepl("^say\\(", s)) next
  try(eval(e, envir = globalenv()), silent = TRUE)
  n_eval <- n_eval + 1L
  if (have()) break
}
say("transplant: evaluated %d of %d top-level expressions from %s", n_eval, length(EX), basename(SRC))
stopifnot(have(), is.function(cell_ratio), is.function(cell_tru), is.function(cell_amp))
say("the reported sets in scope: Hsf1_core n = %d, denominator '%s' n = %d, %d cells, %d timepoints",
    length(HOGS[["Hsf1_core"]]), DENOM[["Hsf1_core"]], length(HOGS[[DENOM[["Hsf1_core"]]]]),
    length(CELLS), length(TP))

## ------------------------------------------------------------------ the generalised ratio ----
## At beta = 0 this must be identical to script 04's cell_ratio. The amplitude is measured on the same gene
## set as the divergence it scales, per cell, using script 04's own cell_amp -- so the only thing added here
## is the exponent.
cell_ratio_beta <- function(i_num, i_den, beta) {
  if (beta == 0) return(cell_ratio(i_num, i_den))
  an <- cell_amp(i_num); ad <- cell_amp(i_den)
  num <- cell_tru(i_num) / pmax(an, 1e-9)^beta
  den <- cell_tru(i_den) / pmax(ad, 1e-9)^beta
  ifelse(den > 0 & ad > 0 & an > 0, num / den, NA_real_)
}
boot_ratio_beta <- function(i_num, i_den, beta, nb = NB) {
  out <- matrix(NA_real_, nb, length(TP))
  for (k in seq_len(nb)) {
    r <- cell_ratio_beta(sample(i_num, length(i_num), TRUE), sample(i_den, length(i_den), TRUE), beta)
    for (j in seq_along(TP)) out[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  out
}
summar <- function(i_num, i_den, beta, seed) {
  o <- cell_ratio_beta(i_num, i_den, beta)
  set.seed(seed); b <- boot_ratio_beta(i_num, i_den, beta)
  rbindlist(lapply(seq_along(TP), function(j) {
    sel <- tp_of == TP[j]; ci <- quantile(b[, j], c(.025, .975), na.rm = TRUE)
    data.table(time = TP[j], ratio = median(o[sel], na.rm = TRUE), lo = ci[1], hi = ci[2],
               n_cells = sum(sel), cells_gt1 = sum(o[sel] > 1, na.rm = TRUE))
  }))
}

## ------------------------------------------------------------------ the six definitions ----
CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
HB <- fread(file.path(csd, "hog_to_scer.tsv")); h2s <- unique(HB[, .(hog = gene, sys)])
sys2hog <- function(s) unique(h2s[sys %in% unique(as.character(s)), hog])
UN <- CL[, .(nd = uniqueN(class), cls = class[1]), by = gene][nd == 1L]
med_kd <- CL[class == "Hsf1_only", median(abs(lfc_kd), na.rm = TRUE)]
GOA <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "go_annotations_propagated_universe.tsv.gz"))))
PROTEO <- c("GO:0006457", "GO:0051082", "GO:0031072", "GO:0034620")
universe <- rownames(S_M)
to_idx <- function(h) { i <- match(intersect(h, universe), universe); i[!is.na(i)] }

DEFS <- list(
  reported_core            = HOGS[["Hsf1_core"]],
  D1_unanimous_exclusive = sys2hog(UN[cls == "Hsf1_only", gene]),
  D2_plus_both_required  = sys2hog(UN[cls %in% c("Hsf1_only", "both_required"), gene]),
  D3_any_method          = sys2hog(unique(CL[class == "Hsf1_only", gene])),
  D4_strong_effect_only  = sys2hog(intersect(UN[cls == "Hsf1_only", gene],
                                             CL[class == "Hsf1_only" & abs(lfc_kd) >= med_kd, gene])),
  D5_GO_proteostasis     = sys2hog(unique(GOA[go %in% PROTEO, sys])))
i_den_base <- HOGS[[DENOM[["Hsf1_core"]]]]

## ------------------------------------------------------------------ check 1: beta = 0 reproduces strict ----
say("\n=== fidelity check 1: beta = 0 must reproduce a5_core_definition_strict.tsv ===")
STRICT <- fread(file.path(hsd, "a5_core_definition_strict.tsv"))
worst <- 0
for (k in names(DEFS)) {
  fo <- to_idx(DEFS[[k]]); de <- to_idx(setdiff(i_den_base, DEFS[[k]]))
  if (length(fo) < 20L || length(de) < 20L) next
  o <- cell_ratio_beta(fo, de, 0)
  for (j in seq_along(TP)) {
    ref <- STRICT[definition == k & time == TP[j], ratio]
    if (!length(ref)) next
    worst <- max(worst, abs(median(o[tp_of == TP[j]], na.rm = TRUE) - ref))
  }
}
say("  worst absolute difference across %d definitions x %d timepoints: %.3g", length(DEFS), length(TP), worst)
stopifnot(worst < 1e-9)
say("  reproduced. cell_ratio_beta at beta = 0 IS script 04's cell_ratio.")

## ------------------------------------------------------------------ check 2: headline beta reproduces cal_ratio ----
say("\n=== fidelity check 2: reported_core at beta = 0.383444 must reproduce the published cal_ratio ===")
CAL <- fread(file.path(hsd, "a5_cis_ratio_calibrated_37_30.tsv"))[stratum == "Hsf1_core"]
fo <- to_idx(DEFS$reported_core); de <- to_idx(setdiff(i_den_base, DEFS$reported_core))
o <- cell_ratio_beta(fo, de, 0.383444)
w2 <- 0
for (j in seq_along(TP)) {
  ref <- CAL[time == TP[j], cal_ratio]
  got <- median(o[tp_of == TP[j]], na.rm = TRUE)
  say("  %2d min  published %.6f   reproduced %.6f   diff %.3g", TP[j], ref, got, abs(ref - got))
  w2 <- max(w2, abs(ref - got))
}
say("  worst: %.3g", w2)
stopifnot(w2 < 1e-6)
say("  reproduced. The generalisation lands on the headline numbers at the fitted beta.")

## ------------------------------------------------------------------ the answer ----
say("\n=== six definitions at beta = 0, 0.383444, 1 ===")
RES <- list()
for (k in names(DEFS)) {
  fo <- to_idx(DEFS[[k]]); de <- to_idx(setdiff(i_den_base, DEFS[[k]]))
  if (length(fo) < 20L || length(de) < 20L) { say("  %s skipped (%d focal / %d denom)", k, length(fo), length(de)); next }
  for (be in BETAS_REPORT) {
    z <- summar(fo, de, be, seed = 45L)
    z[, `:=`(definition = k, beta = be, n_focal = length(fo), n_denom = length(de))]
    RES[[length(RES) + 1L]] <- z
  }
}
R <- rbindlist(RES)
R[, `:=`(resolvably_less = ratio < 1 & hi < 1, resolvably_more = ratio > 1 & lo > 1)]
fwrite(R, file.path(hsd, "a5_core_definition_by_beta.tsv"), sep = "\t")

for (be in BETAS_REPORT) {
  say("\n  beta = %.6f", be)
  say("    %-24s %6s %s", "definition", "n", paste(sprintf("%22s", sprintf("%d min", TP)), collapse = ""))
  for (k in unique(R$definition)) {
    z <- R[definition == k & abs(beta - be) < 1e-9][order(time)]
    if (!nrow(z)) next
    cells <- sprintf("%22s", sprintf("%.3f[%.2f,%.2f]%s", z$ratio, z$lo, z$hi,
                     ifelse(z$resolvably_less, " <", ifelse(z$resolvably_more, " >", "  "))))
    say("    %-24s %6d %s", k, z$n_focal[1], paste(cells, collapse = ""))
  }
}
say("\n  '<' marks a CI entirely below 1, '>' entirely above 1.")

say("\n=== the pre-registered read-out ===")
for (be in BETAS_REPORT) {
  g <- R[definition == "D5_GO_proteostasis" & abs(beta - be) < 1e-9]
  say("  GO-only definition, beta %.3f: ratio %s | resolvably less in %d of %d timepoints",
      be, paste(sprintf("%.3f", g$ratio), collapse = " / "), sum(g$resolvably_less), nrow(g))
}
nz <- R[abs(beta - 0.383444) < 1e-9 & resolvably_more == TRUE, .N]
say("\n  definitions x timepoints at the headline beta with a CI ENTIRELY ABOVE 1: %d", nz)
if (nz > 0) print(R[abs(beta - 0.383444) < 1e-9 & resolvably_more == TRUE,
                    .(definition, time, ratio = round(ratio, 3), lo = round(lo, 2), hi = round(hi, 2))],
                 row.names = FALSE)
say("")
