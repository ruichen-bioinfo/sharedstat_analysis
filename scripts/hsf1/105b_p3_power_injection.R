#!/usr/bin/env Rscript
## 105b_p3_power_injection.R -- follow-up to script 105, fixed before this analysis was run:
## (i) cache-key collisions among its random sets, (ii) power of the adjusted and bootstrap intervals after injecting a
## known divergence excess (true ratio 1.25, 1.5, 2.0; 100 draws per cell).
suppressPackageStartupMessages({ library(data.table); library(parallel) })
root <- Sys.getenv("PBL_ROOT"); EX <- parse(file.path(root, "scripts/hsf1/85_p3_parental_axis_tests.R"))
for (e in EX) { s <- paste(deparse(e), collapse = " "); if (grepl("^FAM <- ", s)) break; if (grepl("^sink\\(|^fwrite\\(|quit\\(", s)) next; eval(e, envir = globalenv()) }
## (i) cache-key collisions among the random sets of script 105 (same seeds)
NDRAW <- 200L; SIZES <- c(134L, 83L); set.seed(20261003L)
jobs <- CJ(time = TP, size = SIZES, draw = seq_len(NDRAW)); jobs[, seed := sample.int(1e8, .N)]
keys <- vapply(seq_len(nrow(jobs)), function(i) { j <- jobs[i]; D <- DV_heat[time == j$time]; set.seed(j$seed)
  paste(sort(sample(unique(D$hog), j$size)), collapse = ",") }, "")
dk <- vapply(strsplit(keys, ","), digest_key, "")
cat(sprintf("Gene sets: %d distinct gene sets, %d distinct cache keys (collisions %d)\n", uniqueN(keys), uniqueN(dk), uniqueN(keys) - uniqueN(dk)))
## (ii) power by injection
R_TRUE <- c(1.25, 1.5, 2.0); NP <- 100L; set.seed(20261004L)
pj <- CJ(time = TP, size = SIZES, r = R_TRUE, draw = seq_len(NP)); pj[, seed := sample.int(1e8, .N)]
run <- function(i) { j <- pj[i]; D <- copy(DV_heat[time == j$time]); set.seed(j$seed)
  pick <- sample(unique(D$hog), j$size); inA <- D$hog %in% pick
  D[inA, d_sig := d_sig * sqrt(j$r)]   # raises the set's signal mean square; noise unchanged
  RHO_CACHE <<- list()
  x <- ratio_unc(D, inA, !inA, 0, "heat", seed = j$seed)
  data.table(time = j$time, size = j$size, r = j$r, draw = j$draw, ratio = x$ratio, lo = x$lo, hi = x$hi, lo_boot = x$lo_boot, hi_boot = x$hi_boot) }
P <- rbindlist(mclapply(seq_len(nrow(pj)), run, mc.cores = 3L))
fwrite(P, file.path(hsd, "p3_power_injection_draws.tsv"), sep = "\t")
S <- P[, .(n = .N, median_ratio = median(ratio), power_adjusted = mean(lo > 1), power_boot = mean(lo_boot > 1)), by = .(size, r, time)]
fwrite(S, file.path(hsd, "p3_power_injection.tsv"), sep = "\t"); print(S[order(size, r, time)])
