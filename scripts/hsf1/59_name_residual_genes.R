#!/usr/bin/env Rscript
## 59_name_residual_genes.R -- name representative genes from the 104, taken from the table.
##
## Why names. The P3 metabolism result reports the surviving excess as a ratio over a set of 104 orthologous groups. A reader cannot
## tell from a ratio what kind of gene is involved, and "carbohydrate metabolism outside the Msn2/4 dependence
## class" is a description of a filter, not of biology. Three to five names, chosen by a stated rule and read off
## the saved tables, show what the set contains alongside the full set.
##
## The rule, fixed before the names were seen: rank the 104 by their noise-corrected divergence at 60 min, the
## timepoint at which the excess is resolvable at the measured exponent and at β = 0.75; list the top five with
## their S. cerevisiae gene names, their GO module, whether they are induced in S. cerevisiae, and their
## divergence at each timepoint. Gene recognition does not affect selection.
##
## The sets are rebuilt by re-running script 44's own construction (see script 56 for the same pattern), so the
## 23 / 7 / 104 split remains consistent across analyses.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
outd <- file.path(root, "data/02_processed/crossspecies")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

## ------------------------------------------------------------------ the three sets, from script 44 ----
src44 <- readLines(file.path(root, "scripts/hsf1/44_p3_module_stratum_crosstab.R"))
istop <- grep("^G3 <- ", src44)[1]
e44 <- new.env(parent = globalenv())
invisible(capture.output(suppressWarnings(suppressMessages(
  eval(parse(text = paste(src44[1:istop], collapse = "\n")), envir = e44)))))
msn <- get("msn", envir = e44); hsf <- get("hsf", envir = e44); nei <- get("nei", envir = e44)
DV  <- get("DV", envir = e44); J <- get("J", envir = e44)
stopifnot(length(msn) == 23L, length(hsf) == 7L, length(nei) == 124L)
say("=== sets rebuilt from script 44: Msn2/4 %d, Hsf1 %d, neither %d ===", length(msn), length(hsf), length(nei))
say("  of the 'neither' set, %d have divergence data (the manuscript's 104)", DV[hog %in% nei, uniqueN(hog)])
stopifnot(DV[hog %in% nei, uniqueN(hog)] == 104L)

## ------------------------------------------------------------------ rank by divergence at 60 min ----
say("\n=== the ranking quantity, and what it cannot be ===")
say("  DV columns: %s", paste(names(DV), collapse = ", "))
## The manuscript's divergence is rms_true(d_sig, d_nul), which subtracts a mean square of the null from a mean
## square of the signal. That is a set-level quantity: for a single gene the subtraction has no mean to take and
## can go negative. So a gene cannot be given "its noise-corrected divergence", and ranking on one would be
## inventing a per-gene version of a set estimator. The ranking here is on |d_sig| averaged over the six pairs,
## which is the uncorrected divergence, and the ranking is therefore explicitly a display choice: the noise floor
## is shown beside it so genes near the noise floor are identifiable.
## 
say("  ranking on mean |d_sig| over the six species pairs at 60 min (UNCORRECTED; see the note in this script)")
R60 <- DV[hog %in% nei & time == 60,
          .(div = mean(abs(d_sig), na.rm = TRUE), noise = mean(abs(d_nul), na.rm = TRUE),
            snr = mean(abs(d_sig), na.rm = TRUE) / pmax(mean(abs(d_nul), na.rm = TRUE), 1e-9),
            n_pair = .N), by = hog][order(-div)]
say("  %d genes ranked; |d_sig| from %.3f to %.3f; each has %d pairs",
    nrow(R60), max(R60$div), min(R60$div), R60$n_pair[1])
stopifnot(all(R60$n_pair == 6L))

## names, module and induction status, all from saved tables
HOGS <- fread(file.path(outd, "hog_to_scer.tsv"))
nmcol <- intersect(c("Name", "name", "gene_name"), names(HOGS))[1]
MOD  <- fread(file.path(hsd, "s4_module_membership.tsv"))
STR  <- J[, .SD, .SDcols = intersect(c("hog", "Hsf1_core", "Hsf1_joint", "Msn24_only", "TF_independent",
                                       "not_induced", "other_induced"), names(J))]
TOP <- merge(R60, HOGS[, .(hog = gene, sys, name = get(nmcol))], by = "hog", all.x = TRUE)
TOP <- merge(TOP, MOD[, .(hog, module)], by = "hog", all.x = TRUE)
TOP <- merge(TOP, STR, by = "hog", all.x = TRUE)
TOP[, induced := !(not_induced %in% TRUE)]
setorder(TOP, -div)

## per-timepoint divergence for the ones being named
per <- dcast(DV[hog %in% TOP$hog[1:5], .(hog, time, v = abs(d_sig))], hog ~ time,
             value.var = "v", fun.aggregate = mean)
setnames(per, as.character(c(15, 30, 60)), c("d15", "d30", "d60"), skip_absent = TRUE)
TOP5 <- merge(TOP[1:5], per, by = "hog")
setorder(TOP5, -div)
say("\n=== the five largest |d_sig| among the 104 at 60 min ===")
print(TOP5[, .(name, sys, module, induced, d15 = round(d15, 3), d30 = round(d30, 3), d60 = round(d60, 3),
               noise60 = round(noise, 3), snr60 = round(snr, 2))], row.names = FALSE)

fwrite(TOP[, .(rank = seq_len(.N), hog, sys, name, module, induced,
               mean_abs_d_sig_60 = div, mean_abs_d_nul_60 = noise, snr_60 = snr)],
       file.path(hsd, "residual_104_ranked.tsv"), sep = "\t")
say("\n  landed: residual_104_ranked.tsv (all %d, ranked)", nrow(TOP))

## ------------------------------------------------------------------ interpretation ----
say("\n=== read-out ===")
say("  names, in rank order: %s", paste(TOP5$name, collapse = ", "))
say("  all five carry a gene name: %s", all(!is.na(TOP5$name) & TOP5$name != ""))
say("  all five are in the metabolism module: %s", all(TOP5$module == "metabolism"))
say("  induced in S. cerevisiae: %d of 5", sum(TOP5$induced))
say("  Not shown: that these five drive the excess. The ratio is over all 104; these are the")
say("  most divergent members of it, which is a different statement.")
