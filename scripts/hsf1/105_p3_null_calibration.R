#!/usr/bin/env Rscript
## 105_p3_null_calibration.R -- null calibration of the P3 ratio intervals. The criterion (exclusion rate <= 0.075
## calibrated, > 0.10 anti-conservative, otherwise borderline) was fixed before this analysis was run.
## Random gene sets drawn from the orthologous groups of the module universe (the pool the module family tests use),
## compared with the remaining groups by the paper's own ratio_unc (script 85), beta = 0. True ratio 1 up to sampling.
suppressPackageStartupMessages({ library(data.table); library(parallel) })
root <- Sys.getenv("PBL_ROOT"); EX <- parse(file.path(root, "scripts/hsf1/85_p3_parental_axis_tests.R"))
for (e in EX) { s <- paste(deparse(e), collapse = " "); if (grepl("^FAM <- ", s)) break; if (grepl("^sink\\(|^fwrite\\(|quit\\(", s)) next; eval(e, envir = globalenv()) }
NDRAW <- 200L; SIZES <- c(134L, 83L); set.seed(20261003L)
jobs <- CJ(time = TP, size = SIZES, draw = seq_len(NDRAW)); jobs[, seed := sample.int(1e8, .N)]
run <- function(i) { j <- jobs[i]; D <- DV_heat[time == j$time]; hs <- unique(D$hog); set.seed(j$seed)
  pick <- sample(hs, j$size); inA <- D$hog %in% pick
  x <- ratio_unc(D, inA, !inA, 0, "heat", seed = j$seed)
  data.table(time = j$time, size = j$size, draw = j$draw, ratio = x$ratio, lo = x$lo, hi = x$hi, lo_rho2 = x$lo_rho2,
             hi_rho2 = x$hi_rho2, lo_boot = x$lo_boot, hi_boot = x$hi_boot, se_rho = x$se_rho, se_jk = x$se_jk) }
R <- rbindlist(mclapply(seq_len(nrow(jobs)), run, mc.cores = max(1L, detectCores() - 1L)))
fwrite(R, file.path(hsd, "p3_null_calibration_draws.tsv"), sep = "\t")
S <- R[, .(n = .N, excl_adjusted = mean(lo > 1 | hi < 1), excl_rho2 = mean(lo_rho2 > 1 | hi_rho2 < 1),
           excl_boot = mean(lo_boot > 1 | hi_boot < 1), median_ratio = median(ratio), jk_binding = mean(se_jk >= se_rho)), by = .(time, size)]
fwrite(S, file.path(hsd, "p3_null_calibration.tsv"), sep = "\t"); print(S)
v <- if (all(S$excl_adjusted <= 0.075)) "calibrated for random sets" else if (any(S$excl_adjusted > 0.10)) "anti-conservative" else "borderline"
cat("DECISION:", v, "\n")
