# 03_fay_parental_dimensions.R -- does the heat-stress response differ among the four
# Saccharomyces species, on the Hsf1-associated core specifically?
#
# Scope. Fay et al. 2023 (PMID 37972247) already published that the
# cryophilic species respond more strongly. That is not what this script claims to add.
# The question here follows from 02_core_function.R: the Hsf1-associated core is a narrow
# proteostasis module, so does that module diverge, and does it diverge like the rest of
# the heat-induced response or differently? 04_divergence_by_dependency.R uses the
# answer before it can partition divergence by regulatory dependency.
#
# Projection limitation, which governs every number below. Hsf1
# dependence can only be measured in S. cerevisiae (only the Mühlhofer dataset perturbs Hsf1).
# So what is compared across species is the orthologs of an Sc-defined gene set, not each
# species' own Hsf1 regulon. Regulon membership may itself have changed; no available data
# can exclude that. Nothing here is evidence about any species' Hsf1 protein or activity.
#
# Design, read off design.csv. Parents are 4 species x 2 conditions x
# 4 timepoints (0/15/30/60) x 2 replicate blocks = 32 rows each. There is exactly
# one sample per (species, condition, time, block): the two "replicates" are two complete
# time-course blocks, not two samples at one timepoint. That has one useful consequence
# and one hard limit:
#   useful -- each block yields an independent complete estimate of the whole trajectory, so
#             the between-block spread is a direct, assumption-free measure of how much of an
#             apparent species difference is just replicate noise. This is the calibration the
#             Species differences are therefore shown side by side with within-replicate
#             variation.
#   limit  -- with 2 blocks there is no useful per-species significance test. n=4 species is a
#             structural ceiling. Therefore species-level statements are descriptive only.
#             No species-level p-value is computed because the design does not support one.
#
# Amplitude is measured threshold-free: counts are only an amplitude ordering
# when the noise floor is shared, and it is not shared here (4 species, different mapping
# rates). q90|lfc| and rms(lfc), never "number of DE genes".
#
# Usage: Rscript scripts/hsf1/03_fay_parental_dimensions.R

suppressPackageStartupMessages({ library(data.table); library(edgeR) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
ext  <- file.path(root, "data/03_external/fay2023")
csd  <- file.path(root, "data/02_processed/crossspecies")
hsd  <- file.path(root, "data/02_processed/hsf1")
outd <- hsd
repd <- file.path(root, "doc/05_reports")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

ANCHOR  <- "37_30"
SPECIES <- c("Sc", "Sp", "Sk", "Su")
# Thermal labels are Fay's, and the 2x2 has one strain per cell, so the thermophily x
# domestication interaction is not interpretable as biology. Labels are for reading the
# table, not for a factorial test.
THERMAL <- c(Sc = "thermophilic", Sp = "thermophilic", Sk = "cryophilic", Su = "cryophilic")
ORIGIN  <- c(Sc = "domesticated", Sp = "wild",         Sk = "wild",       Su = "domesticated")
TP      <- c(15L, 30L, 60L)
PRIOR   <- 2       # prior.count for logCPM; damps low-count lfc without thresholding
MINCPM  <- 1       # gene must clear this in every species (see filtering rationale below)

say <- function(...) cat(sprintf(...), "\n", sep = "")
q90 <- function(x) as.numeric(quantile(abs(x), 0.90, na.rm = TRUE))
rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))

logf <- file.path(repd, "logs/hsf1_fay_parental.log")
dir.create(dirname(logf), showWarnings = FALSE, recursive = TRUE)
sink(logf, split = TRUE)
say("=== 03_fay_parental_dimensions | %s | anchor %s C ===", format(Sys.time()), ANCHOR)
say("Fay 2023 PARENTAL level, heat arm. Descriptive at species level.")
say("Projection: orthologs of an Sc-defined Hsf1-associated set, not each species' regulon.")

## --------------------------------------------------------------- inputs ---
# read.csv, not fread: 172 rows of htseq.ortho.csv are ragged (267-269 fields instead of
# 270). All 172 carry one==0 & all4==0, so every ragged row is outside the 5121
# one-to-one ortholog set. read.csv fills the missing trailing fields with NA; those rows
# are dropped immediately.
dat <- read.csv(file.path(ext, "htseq.ortho.csv"))
design <- read.csv(file.path(ext, "design.csv"))
meta_cols <- c("HOG","OG","one","all4","multi","Name","c1","c2","c3","c4",
               "Len.Sc","Len.Sp","Len.Sk","Len.Su")
stopifnot(identical(colnames(dat)[1:14], meta_cols))
stopifnot(nrow(design) == ncol(dat) - 14L)

one <- dat$one == 1
say("ortholog table: %d HOGs; one-to-one in all four species: %d", nrow(dat), sum(one))
stopifnot(sum(one) == 5121L)
cts <- as.matrix(dat[one, -(1:14)])
rownames(cts) <- dat$HOG[one]
stopifnot(!anyNA(cts))          # ragged rows were all one==0; this must hold
say("count matrix after restriction: %d x %d, no NA", nrow(cts), ncol(cts))

D <- as.data.table(design)
D[, col := Name]
D[, block := fifelse(as.integer(sub("^[A-Za-z]+", "", Name)) <= 32L, "block1", "block2")]
P <- D[Strain %in% SPECIES & Allele == Strain & Condition == "heat"]
P[, Time := as.integer(Time)]
say("\nparental heat-arm samples: %d", nrow(P))
print(dcast(P, Strain + block ~ Time, value.var = "col"))
# one sample per cell is the premise of the whole design section above -- assert it
stopifnot(nrow(P) == 32L, !anyDuplicated(P[, .(Strain, block, Time)]))
stopifnot(all(P$col %in% colnames(cts)))

## ------------------------------------------------- gene set assignment ---
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
CS <- fread(file.path(hsd, sprintf("hsf1_core_set_%s.tsv", ANCHOR)))
CL <- fread(file.path(hsd, sprintf("hsf1_epistasis_classes_%s.tsv", ANCHOR)))

# 01_epistasis_decomposition.R assigns classes per method; use the unanimous assignment
# only, so that a stratum label never depends on which DE method happened to be run.
UC <- dcast(CL[, .(gene, method, class)], gene ~ method, value.var = "class")
UC[, unanimous := fifelse(DESeq2 == edgeR & edgeR == limma, DESeq2, NA_character_)]
core  <- CS[set == "Hsf1_core_unanimous",  gene]
joint <- CS[set == "Hsf1_joint_candidate", gene]

# Bridge Sc gene -> HOG. The class table's `gene` column is the standard name where SGD
# has one and the systematic name otherwise, so match on both columns of the bridge.
map_hog <- function(g) {
  h <- HB[Name %in% g | sys %in% g, unique(gene)]
  h[h %in% rownames(cts)]
}
STRATA <- list(
  Hsf1_core        = core,
  Hsf1_joint       = joint,
  Msn24_only       = UC[unanimous == "Msn24_only",     gene],
  TF_independent   = UC[unanimous == "TF_independent",  gene],
  other_induced    = UC[unanimous %in% c("Msn24_only", "TF_independent"), gene],
  not_induced      = UC[unanimous == "not_induced",     gene]
)
HOGS <- lapply(STRATA, map_hog)

say("\n=== projection of each Sc stratum onto one-to-one orthologs ===")
say("%-16s %7s %7s %7s", "stratum", "Sc", "mapped", "frac")
for (s in names(STRATA))
  say("%-16s %7d %7d %6.1f%%", s, length(STRATA[[s]]), length(HOGS[[s]]),
      100 * length(HOGS[[s]]) / length(STRATA[[s]]))

lost <- setdiff(core, unique(c(HB[gene %in% HOGS$Hsf1_core, Name], HB[gene %in% HOGS$Hsf1_core, sys])))
say("\ncore genes with NO one-to-one ortholog (%d of %d) -- these are DROPPED from every",
    length(lost), length(core))
say("cross-species statement below:")
say("  %s", paste(sort(lost), collapse = " "))
say("This loss is NOT random. A one-to-one requirement preferentially removes recently")
say("duplicated genes, which is exactly the paralog-rich part of the core")
say("(HSP and anCWMP families). So the projected core is depleted of duplicated")
say("chaperones. This loss limits any inference that the core is conserved.")

## ------------------------------------------------- per-species lfc ---
# Normalisation is done within species (8 samples each). Rationale: every quantity used
# below is a within-species difference across time, so a species' own library composition
# cancels; pooling species into one TMM fit would instead let one species' composition
# move another's baseline. Cross-species comparability is bought by the gene filter, not
# by a shared normalisation.
say("\n=== per-species logCPM and time-course lfc (heat arm) ===")
keep_by_sp <- list()
LCPM <- list()
LFC  <- list()
NUL  <- list()
for (sp in SPECIES) {
  ps <- P[Strain == sp][order(block, Time)]
  m  <- cts[, ps$col, drop = FALSE]
  y  <- DGEList(counts = m)
  y  <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = PRIOR)
  LCPM[[sp]] <- lcpm
  keep_by_sp[[sp]] <- rowSums(cpm(y) >= MINCPM) >= 4L
  say("  %-3s libs %s | TMM factors %s | genes >= %g CPM in >=4 libs: %d",
      sp, paste(sprintf("%.1fM", colSums(m) / 1e6), collapse = " "),
      paste(sprintf("%.2f", y$samples$norm.factors), collapse = " "),
      MINCPM, sum(keep_by_sp[[sp]]))
  for (bk in c("block1", "block2")) {
    c0 <- ps[block == bk & Time == 0L, col]
    for (t in TP) {
      ct <- ps[block == bk & Time == t, col]
      LFC[[length(LFC) + 1L]] <- data.table(
        hog = rownames(lcpm), species = sp, block = bk, time = t,
        lfc = lcpm[, ct] - lcpm[, c0])
    }
  }
  # Null contrast: same species, same timepoint, the two replicate blocks against each
  # other. By construction it contains no heat signal, so its dispersion is this species'
  # own noise floor -- sequencing depth, library prep and block effects all included.
  for (t in c(0L, TP)) {
    NUL[[length(NUL) + 1L]] <- data.table(
      hog = rownames(lcpm), species = sp, time = t,
      nul = lcpm[, ps[block == "block2" & Time == t, col]] -
            lcpm[, ps[block == "block1" & Time == t, col]])
  }
}
LFC <- rbindlist(LFC)
NUL <- rbindlist(NUL)

# Filter that makes species comparable. A gene detected well in Sc but marginally in Sk
# would show a larger |lfc| in Sk purely from counting noise, which would masquerade as a
# stronger response. Requiring adequate expression in all four species puts every species on
# the same gene set. It costs genes and it is the conservative choice.
keep_all <- Reduce(`&`, keep_by_sp)
say("\ngene filter: expressed in all four species: %d of %d one-to-one orthologs (%.1f%%)",
    sum(keep_all), length(keep_all), 100 * mean(keep_all))
for (sp in SPECIES)
  say("  lost by requiring %s too: %d", sp, sum(Reduce(`&`, keep_by_sp[setdiff(SPECIES, sp)])) - sum(keep_all))
KEEP <- rownames(cts)[keep_all]
LFC <- LFC[hog %in% KEEP]
NUL <- NUL[hog %in% KEEP]

say("\nstratum sizes AFTER the all-species expression filter:")
for (s in names(HOGS)) {
  HOGS[[s]] <- intersect(HOGS[[s]], KEEP)
  say("  %-16s %d", s, length(HOGS[[s]]))
}

## ------------------------------------------------- the four dimensions ---
SETS <- c("Hsf1_core", "Hsf1_joint", "Msn24_only", "TF_independent", "other_induced", "not_induced")
ST <- rbindlist(lapply(SETS, function(s) data.table(hog = HOGS[[s]], stratum = s)))
X  <- merge(LFC, ST, by = "hog", allow.cartesian = TRUE)

# amplitude, per (stratum, species, block, time): threshold-free only
AMP <- X[, .(n = .N, q90 = q90(lfc), rms = rms(lfc), med = median(lfc)),
         by = .(stratum, species, block, time)]

## Noise-floor control; dimension 1 is only interpretable after it.
## The four species were not sequenced to comparable depth: the libraries differ by well
## over an order of magnitude (printed above). A shallow library has a higher noise floor,
## and |lfc| noise adds in quadrature to real signal, so a shallow species can look like it
## responds more strongly for no biological reason at all. Using a threshold-free statistic
## instead of a DE count does not fix this: it is the shared-noise-floor problem in a new
## coordinate.
## The null contrast (block2 - block1 at the same timepoint) measures each species' own floor.
## Both the signal lfc and the null carry the same 2*sigma^2 term, so
##       rms_true = sqrt( max(0, rms_signal^2 - rms_null^2) )
## removes it. The null also absorbs any block effect, which makes the correction
## conservative (it can only shrink an amplitude, never inflate one).
NZ <- merge(NUL, ST, by = "hog", allow.cartesian = TRUE)
FLOOR <- NZ[, .(rms_null = rms(nul), q90_null = q90(nul)), by = .(stratum, species, time)]
FLOOR0 <- NZ[time == 0L, .(rms_null0 = rms(nul)), by = .(stratum, species)]

say("\n================ NOISE FLOOR (prerequisite for DIMENSION 1) ================")
say("Sequencing depth is NOT matched across species, so amplitudes are not comparable until")
say("each species' own floor is measured. Null contrast = block2 - block1 at the SAME")
say("timepoint: no heat signal by construction.")
say("%-16s %-4s %10s %10s %10s %10s", "stratum", "sp", "rms_null_0", "rms_null_15", "rms_null_30", "rms_null_60")
for (s in c("Hsf1_core", "other_induced")) for (sp in SPECIES) {
  v <- sapply(c(0L, TP), function(t) FLOOR[stratum == s & species == sp & time == t, rms_null])
  say("%-16s %-4s %10.3f %10.3f %10.3f %10.3f", s, sp, v[1], v[2], v[3], v[4])
}
say("\nfloor at 0 min (pre-shift, the cleanest read on depth alone), Hsf1_core:")
o0 <- FLOOR0[stratum == "Hsf1_core"][order(-rms_null0)]
say("  %s", paste(sprintf("%s %.3f", o0$species, o0$rms_null0), collapse = " > "))
say("  ratio highest/lowest floor = %.2fx  => a raw amplitude comparison across species is",
    max(o0$rms_null0) / min(o0$rms_null0))
say("  %s", if (max(o0$rms_null0) / min(o0$rms_null0) > 1.5)
  "NOT admissible without the correction below." else "only mildly affected.")

say("\n================ DIMENSION 1: AMPLITUDE ================")
say("q90|lfc| per species, per replicate block. The two blocks are complete independent")
say("time courses, so the block-to-block gap is replicate noise on the SAME quantity.")
for (s in c("Hsf1_core", "other_induced")) {
  say("\n-- %s (n = %d genes) --", s, length(HOGS[[s]]))
  say("%-4s %-7s %8s %8s %8s", "sp", "block", "15min", "30min", "60min")
  for (sp in SPECIES) for (bk in c("block1", "block2")) {
    v <- sapply(TP, function(t) AMP[stratum == s & species == sp & block == bk & time == t, q90])
    say("%-4s %-7s %8.2f %8.2f %8.2f", sp, bk, v[1], v[2], v[3])
  }
}

# The calibration that decides whether any of it is interpretable.
say("\n---- is the between-species spread larger than replicate noise? ----")
say("For each stratum and timepoint: spread across species (block-averaged) vs the median")
say("|block1 - block2| gap within a species. A ratio near or below 1 means NOT RESOLVED.")
CAL <- AMP[, .(mean_q90 = mean(q90)), by = .(stratum, species, time)
           ][, .(sp_spread = max(mean_q90) - min(mean_q90)), by = .(stratum, time)]
GAP <- dcast(AMP, stratum + species + time ~ block, value.var = "q90"
             )[, .(block_gap = median(abs(block1 - block2))), by = .(stratum, time)]
CAL <- merge(CAL, GAP, by = c("stratum", "time"))
CAL[, ratio := sp_spread / block_gap]
CAL[, resolved := fifelse(ratio >= 2, "yes", fifelse(ratio >= 1, "marginal", "NO"))]
say("%-16s %5s %10s %10s %7s %9s", "stratum", "time", "sp_spread", "block_gap", "ratio", "resolved")
for (i in seq_len(nrow(CAL)))
  say("%-16s %5d %10.3f %10.3f %7.2f %9s", CAL$stratum[i], CAL$time[i],
      CAL$sp_spread[i], CAL$block_gap[i], CAL$ratio[i], CAL$resolved[i])

say("\n---- species ordering by q90|lfc|, averaged over blocks ----")
say("Reported to compare with Fay 2023's published result (cryophilic respond more strongly).")
say("This reproduces their finding on these gene sets rather than establishing a new result.")
for (s in c("Hsf1_core", "other_induced", "not_induced")) {
  say("\n-- %s --", s)
  say("%-4s %-13s %-13s %8s %8s %8s", "sp", "thermal", "origin", "15min", "30min", "60min")
  o <- AMP[stratum == s, .(q = mean(q90)), by = .(species, time)]
  for (sp in SPECIES) {
    v <- sapply(TP, function(t) o[species == sp & time == t, q])
    say("%-4s %-13s %-13s %8.2f %8.2f %8.2f", sp, THERMAL[sp], ORIGIN[sp], v[1], v[2], v[3])
  }
  pk <- o[, .(q = mean(q)), by = species][order(-q)]
  say("  ordering (all timepoints pooled): %s", paste(sprintf("%s %.2f", pk$species, pk$q), collapse = " > "))
  th <- o[, .(q = mean(q)), by = .(species)][, .(m = mean(q)), by = .(t = THERMAL[species])]
  say("  cryophilic mean %.2f vs thermophilic mean %.2f (n=2 strains each; DESCRIPTIVE, no test)",
      th[t == "cryophilic", m], th[t == "thermophilic", m])
}

say("\n---- the same ordering AFTER removing each species' own noise floor ----")
say("rms_true = sqrt(rms_signal^2 - rms_null^2). This is the only amplitude comparison that")
say("is admissible across species given the depth imbalance. If an ordering appears here that")
say("is absent above, or vice versa, the raw ordering was a depth artefact.")
SIG <- AMP[, .(rms_sig = mean(rms)), by = .(stratum, species, time)]
COR <- merge(SIG, FLOOR, by = c("stratum", "species", "time"))
COR[, rms_true := sqrt(pmax(0, rms_sig^2 - rms_null^2))]
COR[, snr := rms_sig / rms_null]
for (s in c("Hsf1_core", "other_induced", "not_induced")) {
  say("\n-- %s --", s)
  say("%-4s %-13s %8s %8s %8s %8s", "sp", "thermal", "true15", "true30", "true60", "min_SNR")
  for (sp in SPECIES) {
    v <- sapply(TP, function(t) COR[stratum == s & species == sp & time == t, rms_true])
    sn <- COR[stratum == s & species == sp & time %in% TP, min(snr)]
    say("%-4s %-13s %8.2f %8.2f %8.2f %8.2f", sp, THERMAL[sp], v[1], v[2], v[3], sn)
  }
  pk <- COR[stratum == s & time %in% TP, .(q = mean(rms_true)), by = species][order(-q)]
  say("  corrected ordering: %s", paste(sprintf("%s %.2f", pk$species, pk$q), collapse = " > "))
  raw <- AMP[stratum == s, .(q = mean(q90)), by = species][order(-q)]
  say("  raw ordering was  : %s", paste(raw$species, collapse = " > "))
  say("  %s", if (identical(pk$species, raw$species))
    "SAME ordering before and after correction -> the depth imbalance does not drive it."
    else "ORDERING CHANGES under correction -> the raw ordering is not interpretable.")
}

## Specificity. The corrected table above shows the same species ordering on `not_induced`
## as on the induced strata. That matters more than the ordering itself: if the genes that
## are not heat-induced in S. cerevisiae move just as much, then a large part of what looks
## like "species X responds more strongly" is not specific to the heat programme at all --
## it is how much that species' transcriptome moves in general, plus whatever residual
## depth/normalisation structure survived the correction.
## The fix is an internal control: divide each stratum's reproducible amplitude by the
## not_induced amplitude in the same species. Being a within-species ratio it is immune to
## depth, to library composition and to any global scale factor.
## It also has a second reading. `not_induced` is defined in S. cerevisiae. If those genes
## move in another species, the Sc-defined partition does not transfer -- which is a
## statement about regulon conservation, not an artefact.
say("\n---- specificity: stratum amplitude relative to the not_induced baseline ----")
say("Within-species ratio, so depth and global scale cancel. High = the Sc-defined induced")
say("set moves much more than the Sc-defined non-induced set in that species, i.e. the Sc")
say("partition transfers. Near 1 = it does not transfer.")
BASE <- COR[stratum == "not_induced", .(species, time, base = rms_true)]
SPEC <- merge(COR[stratum != "not_induced"], BASE, by = c("species", "time"))
SPEC[, spec := fifelse(base > 0, rms_true / base, NA_real_)]
for (s in c("Hsf1_core", "Hsf1_joint", "other_induced")) {
  say("\n-- %s --", s)
  say("%-4s %-13s %8s %8s %8s", "sp", "thermal", "15min", "30min", "60min")
  for (sp in SPECIES) {
    v <- sapply(TP, function(t) SPEC[stratum == s & species == sp & time == t, spec])
    say("%-4s %-13s %8s %8s %8s", sp, THERMAL[sp],
        ifelse(is.na(v[1]), "  n/a", sprintf("%8.2f", v[1])),
        ifelse(is.na(v[2]), "  n/a", sprintf("%8.2f", v[2])),
        ifelse(is.na(v[3]), "  n/a", sprintf("%8.2f", v[3])))
  }
}
say("\n'n/a' means the not_induced baseline floored to zero in that cell (all of its apparent")
say("amplitude was explained by that species' own noise), so no ratio is defined. That is the")
say("expected behaviour for a negative control and is reported rather than filled in.")
say("INDEPENDENCE OF THE GENE-SET DEFINITION: the strata were defined in the Mühlhofer dataset -- a")
say("different laboratory, a different S. cerevisiae background and a different heat protocol")
say("from Fay's T73. So the Sc column here is a CROSS-DATASET transfer check, not a")
say("self-fulfilling one. Its value is the benchmark the other three are read against.")

say("\n================ DIMENSION 2: TIMING ================")
say("Peak timepoint of the stratum-level amplitude, per species per block. With only three")
say("post-shift timepoints the resolution is coarse: a peak can only be 15, 30 or 60 min.")
PEAK <- AMP[, .(peak = time[which.max(rms)], rms_peak = max(rms)), by = .(stratum, species, block)]
for (s in c("Hsf1_core", "other_induced")) {
  say("\n-- %s --", s)
  say("%-4s %10s %10s %12s", "sp", "block1", "block2", "agree?")
  for (sp in SPECIES) {
    p1 <- PEAK[stratum == s & species == sp & block == "block1", peak]
    p2 <- PEAK[stratum == s & species == sp & block == "block2", peak]
    say("%-4s %10d %10d %12s", sp, p1, p2, ifelse(p1 == p2, "yes", "NO -> unresolved"))
  }
}
# per-gene peak time is better powered than the stratum summary (n = genes, not n = 4)
PG <- X[, .(peak = time[which.max(abs(lfc))]), by = .(stratum, species, block, hog)]
say("\nper-gene peak-time distribution (pooled over blocks), %% of genes peaking at each time:")
say("%-16s %-4s %7s %7s %7s", "stratum", "sp", "15min", "30min", "60min")
for (s in c("Hsf1_core", "other_induced")) for (sp in SPECIES) {
  d <- PG[stratum == s & species == sp]
  f <- sapply(TP, function(t) 100 * mean(d$peak == t))
  say("%-16s %-4s %6.1f%% %6.1f%% %6.1f%%", s, sp, f[1], f[2], f[3])
}

say("\n================ DIMENSION 3: PERSISTENCE ================")
say("rms(60 min) / rms(peak). CENSORED STATISTIC -- read the caveat before the table.")
say("This ratio is bounded above by 1 and equals EXACTLY 1 whenever the peak falls on the")
say("last observed timepoint. For such a species the response had not yet turned over inside")
say("the 60 min window, so the value is not a measurement of persistence at all: timing and")
say("persistence are not separable there. Those cells are marked 'censored' and must not be")
say("compared with the others, nor averaged into a species contrast.")
PER <- merge(AMP[time == 60L, .(stratum, species, block, rms60 = rms)],
             PEAK, by = c("stratum", "species", "block"))
PER[, persist := rms60 / rms_peak]
PER[, censored := peak == 60L]
for (s in c("Hsf1_core", "other_induced")) {
  say("\n-- %s --", s)
  say("%-4s %10s %10s %10s   %s", "sp", "block1", "block2", "mean", "status")
  for (sp in SPECIES) {
    d <- PER[stratum == s & species == sp][order(block)]
    say("%-4s %10.3f %10.3f %10.3f   %s", sp, d$persist[1], d$persist[2], mean(d$persist),
        if (all(d$censored)) "CENSORED (peak at 60 min in both blocks)"
        else if (any(d$censored)) "partly censored" else "measurable")
  }
  U <- PER[stratum == s & !censored]
  say("  species with a measurable turnover: %d of 4", U[, uniqueN(species)])
  if (U[, uniqueN(species)] < 3L) {
    say("  => too few uncensored species to compare persistence. not estimable as a separate dimension.")
    say("     The censoring pattern itself is informative about")
    say("     timing: a peak at 60 min means the window ended too early.")
  } else {
    pr <- U[, .(p = mean(persist)), by = species]
    gp <- U[, .(g = abs(diff(persist))), by = species][, median(g)]
    say("  among uncensored species only: range %.3f, block gap %.3f -> %s",
        diff(range(pr$p)), gp, ifelse(diff(range(pr$p)) > 2 * gp, "resolved", "NOT resolved"))
  }
}

say("\n================ DIMENSION 4: COMPOSITION ================")
say("Does the SAME core gene respond in every species, or do different members carry the")
say("response? Measured as the between-species correlation of per-gene lfc within a stratum,")
say("calibrated against the between-BLOCK correlation of the same species (the ceiling that")
say("replicate noise alone permits). A species pair cannot be more discordant than noise.")
## The statistic is the ratio between-species r / within-species r. Using the ratio rather
## than the raw between-species r is what makes it comparable across strata of different
## size and different noise: both terms are computed on the same genes with the same
## replicate structure, so the shared attenuation cancels. Value 1 = species are as
## concordant as replicates of one species, i.e. no detectable compositional divergence.
## Inference is at the gene level (n = 85 and 769), which is the only level with power here
## -- never at the species level (n = 4).
## Attenuation must be matched on both sides of the ratio. Averaging the two blocks before
## correlating across species halves the noise variance in the between-species term, while
## the within-species ceiling is necessarily computed from two single blocks. Since
## observed r = true r / (1 + sigma^2/S^2) per noisy input, the averaged version is
## attenuated less, so the ratio comes out too high -- species look more conserved than
## they are, which is the wrong direction to be wrong in when the finding under test is a
## conservation claim.
## Fix: compute the between-species correlation block-wise too (block1 vs block1, block2 vs
## block2, then average the two r values). Every correlation entering the ratio is then
## between two single-block measurements carrying the same noise, and the attenuation cancels.
NBOOT <- 2000L
set.seed(20260905)
ratio_stat <- function(w, idx) {
  g <- function(sp, bk) w[[paste0(sp, "_", bk)]][idx]
  # within-species: the replicate ceiling
  wc <- sapply(SPECIES, function(sp) cor(g(sp, "block1"), g(sp, "block2"), use = "pairwise.complete.obs"))
  # between-species: block-matched, so both inputs are single blocks like the ceiling
  pairs <- combn(SPECIES, 2, simplify = FALSE)
  bc <- sapply(pairs, function(p)
    mean(c(cor(g(p[1], "block1"), g(p[2], "block1"), use = "pairwise.complete.obs"),
           cor(g(p[1], "block2"), g(p[2], "block2"), use = "pairwise.complete.obs"))))
  c(within = median(wc), between = median(bc), ratio = median(bc) / median(wc))
}
CMP <- list()
for (s in c("Hsf1_core", "other_induced")) {
  say("\n-- %s (n = %d) --", s, length(HOGS[[s]]))
  W <- dcast(X[stratum == s], hog + time ~ species + block, value.var = "lfc")
  for (t in TP) {
    w <- W[time == t]
    n <- nrow(w)
    pt <- ratio_stat(w, seq_len(n))
    bs <- replicate(NBOOT, ratio_stat(w, sample.int(n, n, replace = TRUE))["ratio"])
    ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
    say("  %d min: within (replicate ceiling) %.3f | between %.3f | ratio %.3f [95%% CI %.3f-%.3f]",
        t, pt["within"], pt["between"], pt["ratio"], ci[1], ci[2])
    CMP[[length(CMP) + 1L]] <- data.table(stratum = s, time = t, n = n,
      within = pt["within"], between = pt["between"], ratio = pt["ratio"],
      lo = ci[1], hi = ci[2])
  }
}
CMP <- rbindlist(CMP)
say("\nCI is a gene-level bootstrap (%d resamples). It does NOT cover species sampling:", NBOOT)
say("with 4 species and one strain each, the species axis has no inferential width at all.")

say("\n---- does the Hsf1-associated core diverge MORE than other heat-induced genes? ----")
say("Compared as ratios, which removes the shared replicate attenuation. The core is the")
say("smaller set (%d vs %d), so a difference could be a sample-size effect; the n-matched",
    CMP[stratum == "Hsf1_core", n[1]], CMP[stratum == "other_induced", n[1]])
say("control below subsamples other_induced to the core's size and asks where the core sits")
say("in that null distribution. This control CAN fail.")
Wc <- dcast(X[stratum == "Hsf1_core"],     hog + time ~ species + block, value.var = "lfc")
Wo <- dcast(X[stratum == "other_induced"], hog + time ~ species + block, value.var = "lfc")
say("%5s %10s %10s %12s %12s %s", "time", "core", "other", "core_CI", "n-matched", "verdict")
MATCH <- list()
for (t in TP) {
  wc_ <- Wc[time == t]; wo_ <- Wo[time == t]
  rc <- ratio_stat(wc_, seq_len(nrow(wc_)))["ratio"]
  ro <- ratio_stat(wo_, seq_len(nrow(wo_)))["ratio"]
  nsub <- nrow(wc_)
  sub <- replicate(NBOOT, ratio_stat(wo_, sample.int(nrow(wo_), nsub, replace = FALSE))["ratio"])
  pemp <- mean(sub <= rc, na.rm = TRUE)          # how often does a same-size random
                                                  # heat-induced set look as divergent?
  ci <- quantile(replicate(NBOOT, ratio_stat(wc_, sample.int(nrow(wc_), nrow(wc_), TRUE))["ratio"]),
                 c(0.025, 0.975), na.rm = TRUE)
  say("%5d %10.3f %10.3f %12s %12s %s", t, rc, ro,
      sprintf("%.2f-%.2f", ci[1], ci[2]),
      sprintf("%.3f (%.0f%%)", median(sub, na.rm = TRUE), 100 * pemp),
      if (pemp < 0.05) "core MORE divergent than size-matched controls"
      else if (pemp > 0.95) "core LESS divergent (more conserved)"
      else "NOT distinguishable from a random same-size heat-induced set")
  MATCH[[length(MATCH) + 1L]] <- data.table(time = t, core_ratio = rc, other_ratio = ro,
    core_lo = ci[1], core_hi = ci[2], matched_median = median(sub, na.rm = TRUE), p_emp = pemp)
}
MATCH <- rbindlist(MATCH)
say("\n'%%' is the fraction of size-matched random heat-induced sets that are at least as")
say("divergent as the core. It is a one-sided empirical p on the CONSERVATION side when large")
say("and on the divergence side when small. This is the quantity used to compare the core with other heat-induced genes.")

## ------------------------------------------------- outputs ---
fwrite(AMP,  file.path(outd, sprintf("fay_parental_amplitude_%s.tsv", ANCHOR)), sep = "\t")
fwrite(CAL,  file.path(outd, sprintf("fay_parental_resolution_%s.tsv", ANCHOR)), sep = "\t")
fwrite(PER,  file.path(outd, sprintf("fay_parental_persistence_%s.tsv", ANCHOR)), sep = "\t")
fwrite(PEAK, file.path(outd, sprintf("fay_parental_peak_%s.tsv", ANCHOR)), sep = "\t")
fwrite(COR,  file.path(outd, sprintf("fay_parental_noise_corrected_%s.tsv", ANCHOR)), sep = "\t")
fwrite(SPEC, file.path(outd, sprintf("fay_parental_specificity_%s.tsv", ANCHOR)), sep = "\t")
fwrite(CMP,  file.path(outd, sprintf("fay_parental_composition_%s.tsv", ANCHOR)), sep = "\t")
fwrite(MATCH,file.path(outd, sprintf("fay_parental_core_vs_matched_%s.tsv", ANCHOR)), sep = "\t")
fwrite(X,    file.path(outd, sprintf("fay_parental_lfc_%s.tsv.gz", ANCHOR)), sep = "\t")
say("\nwrote 5 tables to %s", outd)
say("Species-level results in these files are DESCRIPTIVE (n=4, one strain per cell).")
writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_fay_parental.txt"))
say("=== done | %s ===", format(Sys.time()))
sink()
