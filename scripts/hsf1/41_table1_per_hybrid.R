#!/usr/bin/env Rscript
## 41_table1_per_hybrid.R -- build P1 Table 1 from the saved table, and write it as markdown.
##
## Purpose. This script supplies the Table 1 body in a standalone file
## that can be checked during the document build.
##
## Calibre. Model A, median over the three differential-expression frameworks, the convention locked for
## the primary analysis (script 37), so the table reproduces the four reported values.
##
## SOURCE: data/02_processed/crossspecies/gain_iv.tsv

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
say  <- function(...) cat(sprintf(...), "\n")

d <- fread(file.path(root, "data/02_processed/crossspecies/gain_iv.tsv"))
a <- d[model == "A"]
stopifnot(nrow(a) == 18L, length(unique(a$method)) == 3L, length(unique(a$group)) == 6L)

tab <- a[, .(slope = median(slope_ols), gain_iv = median(gain_iv), artefact = median(artefact)),
         by = group][order(group)]
## Report both columns as gains. A naive slope of -0.292 beside a corrected gain of 0.928 invites the
## reader to compare two different quantities; 1 + slope is the gain the naive reading implies, so the two
## columns then say the same kind of thing and the table speaks for itself.
tab[, naive := 1 + slope]
tab[, artefact_pct := 100 * abs(artefact) / abs(slope)]
## the two orderings are by the strength of the inferred divergence: |naive slope| under the naive
## reading, |gain - 1| under the instrumented one
tab[, rank_naive := rank(-abs(naive - 1))]
tab[, rank_iv    := rank(-abs(gain_iv - 1))]

say("=== fidelity gate against the four reported primary-analysis values ===")
gate <- function(what, got, want, tol) {
  ok <- is.finite(got) && abs(got - want) <= tol
  say("  %-52s %9.4f vs %8.4f  %s", what, got, want, if (ok) "OK" else "MISMATCH")
  if (!ok) stop(sprintf("fidelity gate failed: %s", what))
}
gate("smallest artefact share (%)",  min(tab$artefact_pct), 30.369, 5e-3)   # reading (i) of 11_gain_iv.R: no significance filter
gate("largest artefact share (%)",   max(tab$artefact_pct), 117.067, 5e-3)
gate("rank correlation of the two orderings",
     cor(tab$rank_naive, tab$rank_iv, method = "spearman"), 0.3143, 5e-4)
gate("ScxSk naive rank",  tab[group == "ScxSk"]$rank_naive, 1, 1e-9)
gate("ScxSk naive gain",  tab[group == "ScxSk"]$naive, 0.6487, 5e-4)
gate("ScxSp corrected gain (most divergent)", tab[group == "ScxSp"]$gain_iv, 0.8016, 5e-4)   # gain = 1 + slope_IV, reading (i)
gate("SpxSu corrected gain (crosses one)",    tab[group == "SpxSu"]$gain_iv, 1.0195, 5e-4)   # gain = 1 + slope_IV, reading (i)
gate("ScxSk instrumented rank", tab[group == "ScxSk"]$rank_iv, 5, 1e-9)
## the one hybrid whose inferred direction changes: the naive slope is negative, the corrected gain is
## above one, so a reported decrease becomes an increase
## With both columns as gains, the reversal is "a reported attenuation that is really none":
## naive gain below one, corrected gain above it.
sign_flip <- tab[naive < 1 & gain_iv > 1]
say("  %-52s %s", "hybrid whose direction changes", paste(sign_flip$group, collapse = ", "))
stopifnot(nrow(sign_flip) == 1L, sign_flip$group == "SpxSu")

## species abbreviations spelled out, because a table should be readable without the Methods
NAME <- c(ScxSk = "S. cerevisiae × S. kudriavzevii", ScxSp = "S. cerevisiae × S. paradoxus",
          ScxSu = "S. cerevisiae × S. uvarum",       SkxSu = "S. kudriavzevii × S. uvarum",
          SpxSk = "S. paradoxus × S. kudriavzevii",  SpxSu = "S. paradoxus × S. uvarum")

## the instrument-contamination breakdown (12_iv_breakdown.R): the cor(S, u) at which the corrected direction
## (gain above or below one) would be undone. Arithmetic on the saved IV slopes; see that script for the algebra.
BK <- fread(file.path(root, "data/02_processed/crossspecies/iv_breakdown.tsv"))
tab <- merge(tab, BK[, .(group, rho_star)], by = "group", all.x = TRUE)
stopifnot(!anyNA(tab$rho_star))
fmt_rho <- function(r) sub("^-", "\u2212", sprintf("%+.2f", r))
lines <- c(
 "| Hybrid | Naive gain | Corrected gain | Artefact, % | Rank, naive | Rank, corrected | Breakdown cor(*S*, *u*) |",
 "|---|---|---|---|---|---|---|",
 sprintf("| *%s* | %.3f | %.3f | %.0f | %d | %d | %s |",
         NAME[tab$group], tab$naive, tab$gain_iv, tab$artefact_pct,
         as.integer(tab$rank_naive), as.integer(tab$rank_iv), fmt_rho(tab$rho_star)))
## the third dataset's row from 54_p1_hovhannisyan_dataset.R. No instrument is
## ranked for it.
hv <- file.path(root, "data/02_processed/hsf1/p1x_hov2020_table1_row.tsv")
if (file.exists(hv)) { h3 <- fread(hv)
  lines <- c(lines, sprintf("| *S. cerevisiae × S. uvarum* (Hovhannisyan, 3 replicates) | %.3f | %.3f | %.0f | — | — | — |",
                            h3$naive, h3$gain_iv, h3$artefact_pct)) }

out <- file.path(root, "doc/04_figures/final/TABLE1_per_hybrid.md")
writeLines(lines, out)
say("\n=== written ===")
say("  %s  (%d rows)", basename(out), nrow(tab))
cat("\n"); cat(paste(lines, collapse = "\n"), "\n")
