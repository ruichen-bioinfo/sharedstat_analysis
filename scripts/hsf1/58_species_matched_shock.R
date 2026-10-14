#!/usr/bin/env Rscript
## 58_species_matched_shock.R -- module divergence when each species is shocked at its own temperature.
##
## Question. P3's surviving excess is measured under a common 37 C, which is not the same intensity of stress
## for a thermophile and a cryophile. If part of the excess is that mismatch rather than regulatory rewiring, it
## should shrink when each species is shocked at its own temperature. Keyport Kik et al. did exactly that:
## S. cerevisiae 33 to 40 C and S. kudriavzevii 24 to 31 C, eight minutes, two biological replicates each.
##
## Decision rule, fixed before running:
##   Quantity : module-level divergence of the S. cerevisiae and S. kudriavzevii responses under species-specific
##              shock, using the same five module definitions projected onto S. cerevisiae orthologs; metabolism
##              against the remaining modules; beta grid including the anchor 0.657. Two replicates give the
##              signal and null contrasts.
##   Branch (a): the excess shrinks, interval covering 1: then part of the excess under a common 37 C is
##              ecological distance.
##   Branch (b): the excess persists: then the rewiring reading stands.
##   Power statement: eight minutes corresponds to this project's 15 min, the timepoint at which the heat arm is
##              not resolvable. So a negative in branch (a) is informative only if the interval is narrow; a wide
##              interval is reported as no information and is not evidence in either direction.
##   Not done  : changing the module definitions, changing beta, or changing MINCPM to obtain agreement.
##
## K. marxianus is excluded: it sits outside Saccharomyces, carries a single Msn2/4 paralogue, and including it
## would change the phylogenetic scale fixed for this project.
##
## Note: the ratios in this table are not comparable with P3's. The magnitude used for the beta normalisation here is
## the S. cerevisiae response alone, because that is the only reference arm a two-species design has. P3
## normalises by the mean of the two allele responses inside a hybrid. The two are different denominators, so a
## value of 1.15 here and a value of 1.146 in the P3 metabolism result are not the same quantity and cannot validly be
## placed side by side. This result does not provide evidence in either direction, because the result below carries no information
## either way, but the table is saved and the caveat applies to anyone reading it.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
ext  <- file.path(root, "data/03_external/keyportkik2024")
outd <- file.path(root, "data/02_processed/crossspecies")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
MINCPM <- 1; PRIOR <- 2; NBOOT <- 2000L; SEED <- 20260913L
## The heat-arm table s4_module_divergence_beta.tsv this is compared against uses the grid 0 / 0.25 / 0.5 /
## 0.75 / 1 and does not contain 0.657 (taking the anchor from it returns an empty row), so the grid here
## matches the table's own grid and the comparison is made at values that exist on both sides.
BETAS <- c(0, 0.25, 0.5, 0.75, 1)
ANCHOR <- 0.75

## ------------------------------------------------------------------ design, asserted ----
say("=== design, from the per-sample GEO records ===")
gsm <- readLines(file.path(ext, "GSE234499_samples.txt"), warn = FALSE)
G <- data.table(
  gsm  = sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", gsm, value = TRUE)),
  title = sub("^!Sample_title = ", "", grep("^!Sample_title", gsm, value = TRUE)),
  temp = as.numeric(sub("^!Sample_characteristics_ch1 = temperature: ", "",
                        grep("^!Sample_characteristics_ch1 = temperature:", gsm, value = TRUE))),
  br   = as.integer(sub("^!Sample_characteristics_ch1 = biological replicate: ", "",
                        grep("^!Sample_characteristics_ch1 = biological replicate:", gsm, value = TRUE))))
G[, species := sub("_.*$", "", title)]
G[, `:=`(t_lo = min(temp), t_hi = max(temp)), by = species]
G[, cond := ifelse(temp == t_lo, "ctrl", "shock")]
say("  shifts: %s", paste(unique(sprintf("%s %d->%d C", G$species, G$t_lo, G$t_hi)), collapse = "; "))
stopifnot(G[, .N, by = .(species, cond, br)][, all(N == 1L)], uniqueN(G$species) == 3L)
say("  two replicates in each of ctrl and shock, for each species: asserted")
say("  K. marxianus is present in the data and is EXCLUDED by pre-registration")

## ------------------------------------------------------------------ counts, mapped ----
MAP <- fread(file.path(outd, "skdz_to_scer_map.tsv"))
MAP[module == "" , module := NA_character_]     # fread turns the empty field into "", not NA
say("\n=== the ortholog map (script 57) ===")
say("  %d pairs, %d with a module", nrow(MAP), MAP[!is.na(module), .N])
stopifnot(nrow(MAP) > 4000L, MAP[!is.na(module), .N] > 1000L, MAP[!is.na(module), .N] < nrow(MAP))

X <- fread(cmd = paste("zcat", shQuote(file.path(ext, "GSE234499_counts.tsv.gz"))))
X <- X[!grepl("^__", ORF) & species %in% c("Scer", "Skud")]
SM <- data.table(sample = paste0("SK", 1:12),
                 cond = rep(c("ctrl", "ctrl", "shock", "shock"), 3),
                 br = rep(c(1L, 2L, 1L, 2L), 3))
X <- merge(X, SM, by = "sample")
X[, cpm := 1e6 * counts / sum(counts), by = sample]
X[, lcpm := log2(cpm + PRIOR)]
## put both species on the S. cerevisiae identifier
X[species == "Skud", key := MAP[match(ORF, skdz), sys]]
X[species == "Scer", key := ORF]
say("\n=== mapping the counts onto S. cerevisiae identifiers ===")
say("  S. kudriavzevii rows mapped: %.1f%%", 100 * mean(!is.na(X[species == "Skud", key])))
X <- X[!is.na(key)]
ok <- X[, .(k = sum(cpm >= MINCPM)), by = .(key, species)][k >= 2L]
keep <- intersect(ok[species == "Scer", key], ok[species == "Skud", key])
keep <- intersect(keep, MAP[!is.na(module), sys])
say("  genes expressed in both species and carrying a module: %d", length(keep))
X <- X[key %in% keep]
stopifnot(length(keep) > 700L)

## ------------------------------------------------------------------ responses and contrasts ----
W <- dcast(X, key + species ~ cond + br, value.var = "lcpm")
R <- W[, .(key, species,
           d1 = shock_1 - ctrl_1,
           d2 = shock_2 - ctrl_2)]
R <- R[is.finite(d1) & is.finite(d2)]
## Divergence between the two species' responses, per replicate block, then the signal and null contrasts of
## that difference. This is the same construction P3 uses on the Fay design: signal is the mean over blocks,
## null is the half-difference, and the noise-corrected amplitude subtracts the null mean square.
DV <- merge(R[species == "Scer", .(key, sc1 = d1, sc2 = d2)],
            R[species == "Skud", .(key, sk1 = d1, sk2 = d2)], by = "key")
DV[, `:=`(b1 = sk1 - sc1, b2 = sk2 - sc2)]
DV[, `:=`(sig = (b1 + b2) / 2, nul = (b2 - b1) / 2)]
DV <- merge(DV, MAP[!is.na(module), .(key = sys, module)], by = "key")
say("\n=== divergence of the two responses ===")
say("  genes: %d", nrow(DV))
print(DV[, .(genes = .N, rms_sig = round(sqrt(mean(sig^2)), 3), rms_nul = round(sqrt(mean(nul^2)), 3)),
         by = module][order(-genes)], row.names = FALSE)
h_all <- 1 - mean(DV$nul^2) / mean(DV$sig^2)
say("  reliability of the divergence contrast, h = %.3f", h_all)
if (h_all <= 0) say("  WARNING: the null exceeds the signal; the divergence is not resolved from noise at all")

## ------------------------------------------------------------------ the ratio, on the beta grid ----
## beta is the exponent on the response magnitude in the normalisation, exactly as in P3: a value of 0 is no
## normalisation and 1 is full normalisation by the S. cerevisiae response.
rms_true <- function(s, n) sqrt(max(0, mean(s^2) - mean(n^2)))
ratio_beta <- function(D, be, idx = NULL) {
  d <- if (is.null(idx)) D else D[idx]
  d <- copy(d)
  amp <- pmax(abs((d$sc1 + d$sc2) / 2), 0.05)      # the S. cerevisiae response magnitude
  d[, `:=`(s = sig / amp^be, n = nul / amp^be)]
  m <- d[module == "metabolism", rms_true(s, n)]
  o <- d[module != "metabolism", rms_true(s, n)]
  if (!is.finite(m) || !is.finite(o) || o == 0) return(NA_real_)
  m / o
}
set.seed(SEED)
gene <- unique(DV$key); gi <- split(seq_len(nrow(DV)), DV$key)
OUT <- rbindlist(lapply(BETAS, function(be) {
  obs <- ratio_beta(DV, be)
  bt <- replicate(NBOOT, {
    idx <- unlist(gi[sample(gene, length(gene), TRUE)], use.names = FALSE)
    ratio_beta(DV, be, idx)
  })
  data.table(beta = be, ratio = obs,
             lo = quantile(bt, .025, na.rm = TRUE), hi = quantile(bt, .975, na.rm = TRUE),
             n_na = sum(!is.finite(bt)))
}))
OUT[, `:=`(width = hi - lo, covers_1 = lo <= 1 & hi >= 1)]
say("\n=== metabolism against the other four modules, species-specific shock, 8 min ===")
print(OUT[, .(beta, ratio = round(ratio, 3), ci = sprintf("[%.3f, %.3f]", lo, hi),
              width = round(width, 3), covers_1, degenerate_boots = n_na)], row.names = FALSE)
fwrite(OUT, file.path(hsd, "species_matched_module_ratio.tsv"), sep = "\t")

## ------------------------------------------------------------------ comparison with the common-temperature analysis ----
P3 <- fread(file.path(hsd, "s4_module_divergence_beta.tsv"))
say("\n=== the same quantity under a COMMON 37 C (this project's heat arm, script s4) ===")
cmp <- P3[module == "metabolism" & abs(beta - ANCHOR) < 1e-9,
          .(time, ratio = round(ratio, 3), ci = sprintf("[%.3f, %.3f]", lo, hi), degenerate)]
print(cmp, row.names = FALSE)
a <- OUT[abs(beta - ANCHOR) < 1e-9]
p15 <- P3[module == "metabolism" & abs(beta - ANCHOR) < 1e-9 & time == 15]
stopifnot(nrow(a) == 1L, nrow(p15) == 1L)
say("\n  side by side at every beta on the shared grid, against the 15 min heat arm:")
SBS <- merge(OUT[, .(beta, matched = ratio, m_lo = lo, m_hi = hi, m_w = width)],
             P3[module == "metabolism" & time == 15, .(beta, common = ratio, c_lo = lo, c_hi = hi)],
             by = "beta")
SBS[, c_w := c_hi - c_lo]
print(SBS[, .(beta, matched = round(matched, 3), matched_ci = sprintf("[%.2f, %.2f]", m_lo, m_hi),
              w_matched = round(m_w, 3), common15 = round(common, 3),
              common_ci = sprintf("[%.2f, %.2f]", c_lo, c_hi), w_common = round(c_w, 3),
              matched_wider = m_w > c_w)], row.names = FALSE)

say("\n=== read-out against the pre-registration ===")
say("  branch (a) shrinks with an interval covering 1 : %s", a$covers_1)
say("  branch (b) persists, interval excluding 1      : %s", !a$covers_1 && a$lo > 1)
say("  POWER: the pre-registration says a negative here is informative only if the interval is narrow.")
say("  interval width at the anchor: %.3f; the common-37 C 15 min width is %.3f",
    a$width, p15$hi - p15$lo)
if (a$covers_1 && a$width > (p15$hi - p15$lo)) {
  say("  => the interval is WIDER than the comparison it would have to beat, so this is NO INFORMATION")
  say("     and does not provide evidence that the excess shrinks.")
} else if (a$covers_1) {
  say("  => interval covers 1 and is not wider than the comparison: branch (a), reportable")
} else {
  say("  => interval excludes 1: branch (b), the excess persists under species-specific shock")
}
say("\n  n = %d genes; %d in metabolism", nrow(DV), DV[module == "metabolism", .N])
say("  two replicates only, so the divergence reliability (h = %.3f) is measured from a single pair of blocks",
    h_all)
say("")
