#!/usr/bin/env Rscript
## Writes the seven-column 52-substrate resolution bound plotted in Figure 8 panel c,
## the bound-conversion analysis reads as `c_star`, and P3 Supplementary Figure S1c reports.
##
## Note: `10_promoter_hse_architecture_v2.R` writes a different, five-column summary over its own 10 substrates. The
## two tables share no column and either parses, so they are kept under different file names;
## The bound-conversion analysis needs the `p90` column of this table.
##
## The 52-substrate pooling (script 10's 10 substrates plus the three shards) is done here and not in `19_figures.R`,
## because the statistic is cited by three other files and the statistic must be available without running figure scripts.
##
## This script does the pooling, writes the seven columns, and does not overwrite the reference table if the values
## differ: a discrepancy is reported, not replaced.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(paste0("  ", fmt, "\n"), ...))

CANON <- file.path(hsd, "b3v2_power_resolution.tsv")
SCHEMA <- c("time", "n", "n_never", "min", "median", "p90", "max")
TP <- c(15L, 30L, 60L)

## ---- pool the substrates -----------------------------------------------------------------------------------------
## Script 10's calibration covers its 10 substrates; the three shards carry the rest. The pooled count is asserted,
## because 52 is cited in the Figure 8 caption and in P3's provenance table.
cal <- file.path(hsd, "b3v2_power_calibration.tsv")
if (!file.exists(cal)) stop(sprintf("missing %s: run 10_promoter_hse_architecture_v2.R first", basename(cal)))
P0 <- fread(cal)
dcol <- intersect(c("detected_infer_alpha", "detected"), names(P0))[1]
if (is.na(dcol))
  stop(sprintf("%s has no detection column; it has: %s", basename(cal), paste(names(P0), collapse = ", ")))
P0 <- P0[, .(substrate = paste0("s0_", substrate), c_inj, time, detected = get(dcol))]

shards <- lapply(1:3, function(s) {
  f <- file.path(hsd, sprintf("b3v2_resolution_shard%d.tsv", s))
  if (!file.exists(f)) { say("shard %d absent", s); return(NULL) }
  d <- fread(f)
  need <- c("substrate", "c_inj", "time", "detected")
  if (!all(need %in% names(d)))
    stop(sprintf("shard %d lacks %s", s, paste(setdiff(need, names(d)), collapse = ", ")))
  d[, ..need]
})
PWR <- rbindlist(c(list(P0), shards), use.names = TRUE)
n_sub <- uniqueN(PWR$substrate)
say("pooled %d substrates from the calibration set and %d shard(s)", n_sub, sum(!vapply(shards, is.null, TRUE)))
if (n_sub != 52L)
  stop(sprintf("pooled %d substrates rather than the required 52; the resolution bound is not computed",
               n_sub))

## ---- minimum detectable coupling per substrate and timepoint ----------------------------------------------------
MINC <- PWR[c_inj > 0 & detected == TRUE, .(min_c = min(c_inj)), by = .(substrate, time)]
ALL  <- CJ(substrate = unique(PWR$substrate), time = TP)
MINC <- merge(ALL, MINC, by = c("substrate", "time"), all.x = TRUE)

RES <- MINC[, .(n       = sum(!is.na(min_c)),
                n_never = sum(is.na(min_c)),
                min     = suppressWarnings(min(min_c, na.rm = TRUE)),
                median  = median(min_c, na.rm = TRUE),
                p90     = suppressWarnings(as.numeric(quantile(min_c, 0.90, na.rm = TRUE, names = FALSE))),
                max     = suppressWarnings(max(min_c, na.rm = TRUE))), by = time]
setorder(RES, time)
setcolorder(RES, SCHEMA)
stopifnot(identical(names(RES), SCHEMA), nrow(RES) == length(TP))
stopifnot(RES[, all(n + n_never == n_sub)])
say("n + n_never equals the pooled substrate count in every row")
print(RES, row.names = FALSE)

## ---- compare with the reference table --------------------------------------------------
if (!file.exists(CANON)) {
  fwrite(RES, CANON, sep = "\t")
  say("wrote reference table %s", basename(CANON))
  quit(save = "no", status = 0)
}

OLD <- fread(CANON)
if (!identical(sort(names(OLD)), sort(SCHEMA))) {
  say("Schema mismatch: the file on disk does not have the required columns.")
  say("  on disk : %s", paste(names(OLD), collapse = ", "))
  say("  expected: %s", paste(SCHEMA, collapse = ", "))
  say("The file contains the distinct five-column script-10 product, expected at")
  say("b3v2_power_resolution_script10_10substrate.tsv. The reference table is unchanged.")
  quit(save = "no", status = 1)
}

setkey(OLD, time); setkey(RES, time)
diffs <- rbindlist(lapply(SCHEMA[-1], function(cl) {
  m <- merge(OLD[, .(time, old = get(cl))], RES[, .(time, new = get(cl))], by = "time")
  m[, field := cl][!is.finite(old) | !is.finite(new) | abs(old - new) > 1e-9]
}), use.names = TRUE, fill = TRUE)

if (nrow(diffs) == 0L) {
  say("field by field, the computed table equals the reference table in all %d values", 3 * (length(SCHEMA) - 1))
  fwrite(RES, CANON, sep = "\t")
  say("the reference table has the required values and a declared producer")
} else {
  say("The computed table differs from the reference table; the reference table is unchanged.")
  print(diffs[, .(field, time, old, new, delta = new - old)], row.names = FALSE)
  say("Supplementary Figure S1c and the bound-conversion analysis's c_star use the canonical values.")
  say("This discrepancy is a provenance gap; the computed output is not equivalent to the cited values.")
  quit(save = "no", status = 1)
}
