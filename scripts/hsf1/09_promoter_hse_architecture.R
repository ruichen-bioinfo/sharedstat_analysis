# 09_promoter_hse_architecture.R -- does target-side HSE / promoter architecture predict cis divergence?
#
# The question is a competing-mechanism elimination, not a search for interesting HSE mutations
# (that framing builds a story out of whatever variation exists):
#
#     Can target-side HSE / promoter architecture explain residual regulatory divergence
#     that is not explained by Hsf1 protein evolution?
#
# Rationale. The Hsf1 DNA-binding domain is 97.2% identical across the four species (script 08),
# so the protein-centric model (Model C) is weakened, and a conserved regulator produces divergent
# system-level output. Two candidate resolutions remain. Model B: the targets' own promoters
# changed (cis). Model A: what feeds into Hsf1 changed (upstream/trans). This script tests Model B
# directly, on the cis axis, because cis divergence is what a promoter change must produce. The
# cis axis per gene comes from script 05; the question here is whether promoter architecture
# predicts it.
#
# Design: two arms distinguish the models, and their
# contrast is the result, not either arm alone:
#     arm 1  does HSE/promoter divergence predict cis divergence?    (Model B's prediction)
#     arm 2  does the same promoter divergence predict trans divergence?  (it should not)
# A promoter change acts on the allele it sits on, so it must show up in cis and not in trans.
# If promoter divergence predicts cis but not trans, that supports Model B. If it predicts both
# equally, that is not Model B: promoter divergence is then tracking some third thing
# (expression level, mutation rate, gene-level constraint), and the apparent cis link is
# confounded. Arm 2 exists for this contrast.
#
# Three controls that can fail, declared before the numbers:
#   (1) Positive control on the scan. Hsf1-core promoters must carry more HSEs than non-induced
#       promoters, a known fact external to anything this test wants to show. If the scan
#       cannot recover it, the scan is not measuring HSEs and negative results are not interpretable,
#       and the finding is "scan not validated".
#   (2) Power calibration. Inject a known coupling between promoter divergence and cis
#       divergence at increasing strength; find the smallest coupling this test can detect.
#       A null result below that threshold is a statement about the test, not the biology.
#   (3) Amplitude matching instead of amplitude normalisation. The rms of a divergence axis
#       scales with response amplitude, so any group comparison is confounded if the groups
#       respond differently. Here the groups (HSE turnover vs not) are matched on
#       response-amplitude decile and expression decile, which removes the confound without
#       needing the noise-corrected amplitude that the script-05 output does not carry.
#
# HSE definition, fixed in advance and quoted from the two sources:
#   Yamamoto 2005 (PMID 15647283, J Biol Chem): the HSE "consists of contiguous inverted
#     repeats of the sequence 5'-nGAAn-3'"; and a novel class has "three direct repeats of
#     nTTCn (or nGAAn) interrupted by 5 bp", where "the number and spacing of the repeating
#     units are critical determinants".
#   Pincus 2018 (PMID 30332327, Mol Biol Cell): Hsf1 binds 74 loci in acute heat shock, linked
#     to 46 genes with strong Hsf1-dependent expression; HSE strength interacts with
#     nucleosome positioning to tune expression.
# Operationalised exactly as those two sentences state, and not extended beyond them:
#   unit          = 5 bp, nGAAn (orientation +) or nTTCn (orientation -)
#   classic HSE   = >= 3 units at 5-bp periodicity with alternating orientation (= contiguous
#                   inverted repeats)
#   gapped HSE    = 3 units at 10-bp periodicity with the same orientation (= three direct
#                   repeats interrupted by 5 bp)
# No PWM is invented, no score threshold is tuned. Nucleosome context (Pincus's second axis)
# is not modelled: there are no cross-species nucleosome data, and this is a stated limit of
# any negative result.
#
# Promoter window, fixed before the analysis: 600 bp upstream of the annotated start, clipped
# at the nearest upstream annotated feature so the window never overlaps an upstream ORF (the
# constraint comes from the HSF1 locus itself, where the four intergenic distances are
# 644-665 bp; script 07). Achieved window lengths are reported per species and checked for
# comparability; a gene whose window differs widely between species is flagged, because then a
# count difference would be a window-length artefact rather than architecture.
#
# Strain mismatch matters more here than anywhere else and is a limitation of the analysis.
# The promoters come from reference strains DBVPG6765(Sc) / CBS432(Sp) / CR85(Sk) / YJF1449(Su);
# the expression comes from T73 / N17 / CR85 / BMV58. Only S. kudriavzevii shares a strain
# across the two axes. Non-coding sequence is more polymorphic within a species than coding
# sequence, so a promoter difference measured on the reference strain may not exist in the
# strain whose RNA was sequenced. Consequence: every HSE-to-expression link is at the
# species-reference level, and a per-species sensitivity check reports whether the one matched
# species (Sk) behaves like the other three.
# Usage: Rscript scripts/hsf1/09_promoter_hse_architecture.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
gdir <- file.path(root, "data/03_external/fay2023_genomes")
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
repd <- file.path(root, "doc/05_reports")
outd <- hsd

WIN     <- 600L     # promoter window fixed in advance, clipped at the upstream feature
KMER    <- 4L       # k for the alignment-free promoter divergence (non-coding: alignment-free
                    # is preferred; a global aligner on 600 bp of non-coding sequence
                    # would impose gaps that are not homology)
NBOOT   <- 2000L
NPERM   <- 1000L
BETA_CIS <- 0.275   # measured in script 05 for the cis axis; used only for a secondary bracket
TP      <- c(15L, 30L, 60L)
set.seed(20260906L)

GEN <- data.table(
  species = c("Sc", "Sp", "Sk", "Su"),
  strain  = c("DBVPG6765", "CBS432", "CR85", "YJF1449"),
  matched_strain = c(FALSE, FALSE, TRUE, FALSE))
HYB <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")

say <- function(...) cat(sprintf(...), "\n", sep = "")
rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))

logf <- file.path(repd, "logs/hsf1_promoter_hse.log")
dir.create(dirname(logf), showWarnings = FALSE, recursive = TRUE)
sink(logf, split = TRUE)
say("=== 09_promoter_hse_architecture | %s ===", format(Sys.time()))
say("Question: can target-side HSE/promoter architecture explain residual regulatory")
say("divergence NOT explained by Hsf1 protein evolution? Two arms: cis (Model B predicts yes)")
say("and trans (Model B predicts no). The CONTRAST is the result.")
say("HSE definition locked to Yamamoto 2005 (PMID 15647283) + Pincus 2018 (PMID 30332327).")
say("Promoter window pre-registered: %d bp upstream, clipped at nearest upstream feature.", WIN)

## ----------------------------------------------------------- FASTA reader ---
read_fasta <- function(path) {
  ln <- readLines(path); h <- which(startsWith(ln, ">"))
  nm <- sub("^>(\\S+).*", "\\1", ln[h]); s <- h + 1L; e <- c(h[-1] - 1L, length(ln))
  setNames(vapply(seq_along(h), function(i) toupper(paste(ln[s[i]:e[i]], collapse = "")), ""), nm)
}
revcomp <- function(x) {
  m <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  vapply(x, function(s) paste(rev(m[strsplit(s, "")[[1]]]), collapse = ""), "")
}

## ------------------------------------------------------- the HSE scanner ---
# Returns unit positions and orientations, then assembles them into the two locked classes.
scan_units <- function(s) {
  n <- nchar(s); if (n < 5L) return(NULL)
  # a unit occupies p..p+4 (nGAAn); the informative core is p+1..p+3
  core <- substring(s, 2L:(n - 3L), 4L:(n - 1L))
  p <- 1L:(n - 4L)
  fwd <- core == "GAA"; rev_ <- core == "TTC"
  keep <- fwd | rev_
  if (!any(keep)) return(NULL)
  data.table(pos = p[keep], ori = ifelse(fwd[keep], 1L, -1L))
}
count_hse <- function(s) {
  U <- scan_units(s)
  if (is.null(U) || nrow(U) == 0L)
    return(list(n_units = 0L, n_classic = 0L, n_gapped = 0L, max_run = 0L))
  setorder(U, pos)
  ori <- setNames(U$ori, U$pos); have <- U$pos
  # classic: contiguous inverted repeats -> units every 5 bp with alternating orientation
  runs <- integer(0)
  for (st in have) {
    if ((st - 5L) %in% have && ori[as.character(st - 5L)] == -ori[as.character(st)]) next  # not a run start
    L <- 1L; cur <- st
    while ((cur + 5L) %in% have && ori[as.character(cur + 5L)] == -ori[as.character(cur)]) {
      L <- L + 1L; cur <- cur + 5L
    }
    runs <- c(runs, L)
  }
  max_run <- if (length(runs)) max(runs) else 0L
  n_classic <- sum(runs >= 3L)
  # gapped / step type: three direct repeats at 10-bp periodicity, same orientation
  n_gapped <- 0L
  for (st in have) {
    if ((st + 10L) %in% have && (st + 20L) %in% have &&
        ori[as.character(st + 10L)] == ori[as.character(st)] &&
        ori[as.character(st + 20L)] == ori[as.character(st)]) n_gapped <- n_gapped + 1L
  }
  list(n_units = nrow(U), n_classic = n_classic, n_gapped = n_gapped, max_run = max_run)
}
# self-test on constructed sequences: the scanner must find what the definitions say it should
t_classic <- paste0("CCCC", "AGAAT", "ATTCT", "AGAAT", "CCCC")            # 3 alternating units
t_gapped  <- paste0("CCC",  "ATTCA", "CCCCC", "ATTCA", "CCCCC", "ATTCA")  # 3 direct, 5 bp apart
t_none    <- strrep("CGCGCG", 12)
sc_c <- count_hse(t_classic); sc_g <- count_hse(t_gapped); sc_n <- count_hse(t_none)
say("\n---- scanner self-test (must pass or every count below is meaningless) ----")
say("constructed classic HSE (3 alternating nGAAn/nTTCn): max_run=%d n_classic=%d  [expect >=3, >=1]",
    sc_c$max_run, sc_c$n_classic)
say("constructed gapped HSE (3 direct nTTCn, 5 bp gaps) : n_gapped=%d              [expect >=1]",
    sc_g$n_gapped)
say("CG-only control                                    : n_units=%d               [expect 0]",
    sc_n$n_units)
stopifnot(sc_c$max_run >= 3L, sc_c$n_classic >= 1L, sc_g$n_gapped >= 1L, sc_n$n_units == 0L)
say("scanner self-test PASSED")

## -------------------------------------- promoter extraction, all 4 species ---
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
A6 <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
UNIV <- unique(A6$hog)
sysmap <- unique(HB[gene %in% UNIV, .(hog = gene, sys)])
say("\nA6 universe: %d genes ; with a systematic name: %d", length(UNIV), nrow(sysmap))

PROM <- list(); WINSTAT <- list()
for (i in seq_len(nrow(GEN))) {
  sp <- GEN$species[i]; strain <- GEN$strain[i]
  gf <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(gdir, paste0(strain, ".final.gff")))),
              sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gf, c("chr","src","type","start","end","score","strand","frame","attr"))
  # Gene_Name extraction: NA when the record carries no tag (some do not)
  gf[, gname := NA_character_]
  hasn <- grepl("Gene_Name=", gf$attr, fixed = TRUE)
  gf[hasn, gname := sub("^.*Gene_Name=([^;]+).*$", "\\1", attr)]
  gf[, sysn := sub("_.*$", "", gname)]
  genes <- gf[type == "gene"]
  # all annotated features on each chromosome define the clipping boundaries
  feats <- gf[type %in% c("gene", "tRNA", "rRNA", "ncRNA")][, .(chr, start, end)]
  setorder(feats, chr, start)
  fa <- read_fasta(file.path(gdir, paste0(strain, ".fasta")))

  tg <- genes[sysn %in% sysmap$sys & !is.na(sysn)]
  tg <- tg[!duplicated(sysn)]
  out <- vector("list", nrow(tg))
  for (k in seq_len(nrow(tg))) {
    ch <- tg$chr[k]; if (!ch %in% names(fa)) next
    cs <- fa[[ch]]; cl <- nchar(cs)
    fc <- feats[chr == ch]
    if (tg$strand[k] == "+") {
      lo_want <- tg$start[k] - WIN; hi <- tg$start[k] - 1L
      up <- fc[end < tg$start[k]]
      bound <- if (nrow(up)) max(up$end) + 1L else 1L
      lo <- max(lo_want, bound, 1L)
      if (hi < lo) next
      seqp <- substr(cs, lo, hi)
    } else {
      hi_want <- tg$end[k] + WIN; lo <- tg$end[k] + 1L
      dn <- fc[start > tg$end[k]]
      bound <- if (nrow(dn)) min(dn$start) - 1L else cl
      hi <- min(hi_want, bound, cl)
      if (hi < lo) next
      seqp <- revcomp(substr(cs, lo, hi))
    }
    out[[k]] <- data.table(sys = tg$sysn[k], species = sp, prom = seqp, plen = nchar(seqp))
  }
  P <- rbindlist(out)
  PROM[[sp]] <- P
  WINSTAT[[sp]] <- data.table(species = sp, strain = strain, n = nrow(P),
    med_len = median(P$plen), q25 = quantile(P$plen, .25), q75 = quantile(P$plen, .75),
    frac_full = mean(P$plen >= WIN), frac_short = mean(P$plen < 200L))
  say("  %-3s (%-9s) promoters extracted: %5d | median len %3.0f [%3.0f-%3.0f] | full-%dbp %.1f%% | <200bp %.1f%%",
      sp, strain, nrow(P), median(P$plen), quantile(P$plen, .25), quantile(P$plen, .75),
      WIN, 100 * mean(P$plen >= WIN), 100 * mean(P$plen < 200L))
}
PR <- rbindlist(PROM)
WS <- rbindlist(WINSTAT)
fwrite(WS, file.path(outd, "b3_window_stats.tsv"), sep = "\t")

# a gene is usable only if all four species yielded a promoter
ok <- PR[, .N, by = sys][N == 4L, sys]
PR <- PR[sys %in% ok]
say("\ngenes with a promoter in all 4 species: %d", length(ok))
# window comparability per gene: a big between-species length spread would make count
# differences a window artefact rather than architecture
PW <- dcast(PR, sys ~ species, value.var = "plen")
PW[, len_min := pmin(Sc, Sp, Sk, Su)][, len_max := pmax(Sc, Sp, Sk, Su)]
PW[, len_ratio := len_max / pmax(1L, len_min)]
say("per-gene window length ratio (max/min across species): median %.2f, %.1f%% of genes <=2x",
    median(PW$len_ratio), 100 * mean(PW$len_ratio <= 2))
COMPARABLE <- PW[len_ratio <= 2, sys]
say("genes with comparable windows (<=2x spread), used for the architecture tests: %d", length(COMPARABLE))

## ----------------------------------------------------------- HSE counts ---
say("\n---- scanning promoters for the two locked HSE classes ----")
HS <- PR[, { h <- count_hse(prom)
             .(n_units = h$n_units, n_classic = h$n_classic,
               n_gapped = h$n_gapped, max_run = h$max_run, plen = plen) },
         by = .(sys, species)]
# density, because a longer window trivially contains more units
HS[, units_per_kb := 1000 * n_units / pmax(1L, plen)]
fwrite(HS, file.path(outd, "b3_hse_counts.tsv.gz"), sep = "\t", compress = "gzip")
say("HSE counts per species (median over genes):")
say("%-4s %9s %10s %10s %9s %12s", "sp", "n_units", "n_classic", "n_gapped", "max_run", "units/kb")
for (sp in GEN$species) { z <- HS[species == sp]
  say("%-4s %9.0f %10.2f %10.2f %9.0f %12.2f", sp, median(z$n_units), mean(z$n_classic),
      mean(z$n_gapped), median(z$max_run), median(z$units_per_kb)) }
say("(n_classic/n_gapped shown as MEANS: they are mostly 0, so a median would be 0 and hide")
say(" all the signal. Any HSE gene is a minority of the genome, which is expected.)")

## ------------ CONTROL 1: the positive control the scan must pass ---------
say("\n================ CONTROL 1: DOES THE SCAN RECOVER KNOWN HSE ENRICHMENT? ================")
say("Hsf1-core promoters must carry more HSEs than non-induced promoters. This is external to")
say("anything this analysis tests. If it fails, the scan is not measuring HSEs and negative results")
say("are not interpretable; the scan is classified as 'not validated'.")
ST <- unique(A6[, .(hog, stratum)])
# Script 05 wrote stratum with all.x=TRUE, so genes in none of the unanimous classes carry NA.
# Label them rather than letting a blank row appear in the control table.
ST[is.na(stratum) | stratum == "", stratum := "unassigned"]
ST <- merge(ST, sysmap, by = "hog")
HSm <- merge(HS, ST, by = "sys", allow.cartesian = TRUE)
PC <- HSm[, .(n = .N, mean_units = mean(n_units), mean_perkb = mean(units_per_kb),
              frac_any_classic = mean(n_classic > 0), frac_any_gapped = mean(n_gapped > 0),
              mean_maxrun = mean(max_run)), by = .(stratum)]
setorder(PC, -mean_perkb)
say("%-16s %6s %10s %10s %14s %13s %11s", "stratum", "n", "units", "units/kb",
    "any_classic", "any_gapped", "max_run")
for (i in seq_len(nrow(PC))) say("%-16s %6d %10.2f %10.2f %13.1f%% %12.1f%% %11.2f",
  PC$stratum[i], PC$n[i], PC$mean_units[i], PC$mean_perkb[i],
  100 * PC$frac_any_classic[i], 100 * PC$frac_any_gapped[i], PC$mean_maxrun[i])
fwrite(PC, file.path(outd, "b3_positive_control.tsv"), sep = "\t")

# formal test, pooled over the 4 species: core vs not_induced, permutation on genes
core_sys <- unique(ST[stratum == "Hsf1_core", sys]); ni_sys <- unique(ST[stratum == "not_induced", sys])
gene_perkb <- HS[, .(perkb = mean(units_per_kb), anyc = as.integer(any(n_classic > 0)),
                     mrun = mean(max_run)), by = sys]
o_core <- gene_perkb[sys %in% core_sys]; o_ni <- gene_perkb[sys %in% ni_sys]
obs_d <- mean(o_core$perkb) - mean(o_ni$perkb)
obs_c <- mean(o_core$anyc) - mean(o_ni$anyc)
pool <- gene_perkb[sys %in% c(core_sys, ni_sys)]
nulld <- replicate(NPERM, { idx <- sample(nrow(pool), nrow(o_core))
  mean(pool$perkb[idx]) - mean(pool$perkb[-idx]) })
nullc <- replicate(NPERM, { idx <- sample(nrow(pool), nrow(o_core))
  mean(pool$anyc[idx]) - mean(pool$anyc[-idx]) })
say("\ncore (n=%d) vs not_induced (n=%d), permutation %d:", nrow(o_core), nrow(o_ni), NPERM)
say("  units/kb        core %.2f vs not_induced %.2f  (diff %+.2f, p_emp %.4f)",
    mean(o_core$perkb), mean(o_ni$perkb), obs_d, mean(nulld >= obs_d))
say("  any classic HSE core %.1f%% vs not_induced %.1f%%  (diff %+.1f pp, p_emp %.4f)",
    100 * mean(o_core$anyc), 100 * mean(o_ni$anyc), 100 * obs_c, mean(nullc >= obs_c))
say("  mean max_run    core %.2f vs not_induced %.2f", mean(o_core$mrun), mean(o_ni$mrun))
SCAN_OK <- mean(nulld >= obs_d) < 0.05 || mean(nullc >= obs_c) < 0.05
say("\nSCAN VALIDATED: %s", SCAN_OK)
if (!SCAN_OK) {
  say("*** The scan does NOT recover the known HSE enrichment in Hsf1-core promoters. Per the")
  say("*** pre-specified stop condition, the promoter HSE test is classified as 'not feasible' and negative")
  say("*** results about HSE architecture are not interpretable.")
}

## -------------------- promoter divergence per species pair, per gene ------
say("\n---- promoter divergence per gene per species pair ----")
say("Two measures: (a) HSE architecture turnover (classic-HSE presence differs, or |dn_units|)")
say("and (b) alignment-free %d-mer promoter sequence divergence. Alignment-free is the honest", KMER)
say("choice for non-coding: a global aligner would impose gaps that are not homology.")
kmers <- do.call(paste0, expand.grid(rep(list(c("A","C","G","T")), KMER), stringsAsFactors = FALSE))
kprof <- function(s) {
  n <- nchar(s); if (n < KMER) return(setNames(rep(0, length(kmers)), kmers))
  v <- table(factor(substring(s, 1:(n - KMER + 1L), KMER:n), levels = kmers))
  as.numeric(v) / max(1L, sum(v))
}
PRc <- PR[sys %in% COMPARABLE]
KP <- PRc[, .(k = list(kprof(prom))), by = .(sys, species)]
KPw <- dcast(KP, sys ~ species, value.var = "k")
HSw_u <- dcast(HS[sys %in% COMPARABLE], sys ~ species, value.var = "n_units")
HSw_c <- dcast(HS[sys %in% COMPARABLE], sys ~ species, value.var = "n_classic")
HSw_k <- dcast(HS[sys %in% COMPARABLE], sys ~ species, value.var = "units_per_kb")
stopifnot(identical(KPw$sys, HSw_u$sys), identical(HSw_u$sys, HSw_c$sys))

PD <- list()
for (h in HYB) {
  s1 <- substr(h, 1, 2); s2 <- substr(h, 4, 5)
  kd <- vapply(seq_len(nrow(KPw)), function(i)
    sqrt(sum((KPw[[s1]][[i]] - KPw[[s2]][[i]])^2)), 0)
  PD[[h]] <- data.table(sys = KPw$sys, pair = h,
    d_units = abs(HSw_u[[s2]] - HSw_u[[s1]]),
    d_perkb = abs(HSw_k[[s2]] - HSw_k[[s1]]),
    turnover = as.integer((HSw_c[[s1]] > 0) != (HSw_c[[s2]] > 0)),
    both_classic = as.integer(HSw_c[[s1]] > 0 & HSw_c[[s2]] > 0),
    kmer_div = kd)
}
PD <- rbindlist(PD)
fwrite(PD, file.path(outd, "b3_promoter_divergence.tsv.gz"), sep = "\t", compress = "gzip")
say("per pair: median |dn_units| / median kmer_div / %% with classic-HSE turnover")
say("%-8s %14s %14s %16s", "pair", "med|dn_units|", "med kmer_div", "turnover %%")
for (h in HYB) { z <- PD[pair == h]
  say("%-8s %14.1f %14.4f %15.1f%%", h, median(z$d_units), median(z$kmer_div),
      100 * mean(z$turnover)) }

## ------------------- the two arms: cis (Model B) and trans (control) -----
say("\n================ THE TWO ARMS: DOES PROMOTER DIVERGENCE PREDICT CIS, NOT TRANS? ================")
say("Model B predicts promoter divergence -> CIS divergence and NOT trans divergence, because a")
say("promoter change acts on the allele carrying it. If it predicts BOTH, the link is confounded")
say("by a third factor and is NOT evidence for Model B. The contrast is the result.")
A6s <- merge(A6[, .(hog, pair, time, cis, cis_nul, trans, trans_nul, amp)],
             sysmap, by = "hog")
A6s <- unique(A6s)
D <- merge(A6s, PD, by = c("sys", "pair"))
say("joined rows (gene x pair x time): %d ; unique genes %d", nrow(D), uniqueN(D$sys))

# high vs low promoter divergence, split at the median within each pair so the split is not
# driven by which pair is more diverged overall
D[, kmer_hi := as.integer(kmer_div > median(kmer_div, na.rm = TRUE)), by = pair]

# ---- matrix form, because the gene-level bootstrap has to run thousands of times ----
# A data.table join + dcast inside every bootstrap iteration is too slow. Everything below is
# precomputed into [gene x cell] matrices so one bootstrap iteration is a row-subset plus
# colMeans, which is what the estimator needs:
#   rms_tru(group, cell) = sqrt(max(0, mean(x^2) - mean(x_nul^2)))  over the group's genes.
D <- D[!is.na(cis) & !is.na(cis_nul) & !is.na(trans) & !is.na(trans_nul)]
D[, cell := paste(pair, time, sep = "_")]
mk <- function(vn) {
  w <- dcast(D, sys ~ cell, value.var = vn, fun.aggregate = function(z) z[1])
  m <- as.matrix(w[, setdiff(names(w), "sys"), with = FALSE]); rownames(m) <- w$sys; m
}
CIS <- mk("cis"); CISN <- mk("cis_nul"); TRN <- mk("trans"); TRNN <- mk("trans_nul")
AMPA <- abs(mk("amp"))
CELLS <- colnames(CIS); GEN_ <- rownames(CIS)
stopifnot(identical(colnames(CISN), CELLS), identical(rownames(TRN), GEN_))
pair_of <- sub("_[0-9]+$", "", CELLS); tp_of <- as.integer(sub(".*_", "", CELLS))
PAIRS <- unique(pair_of)
grp_mat <- function(vn) {
  w <- dcast(unique(D[, c("sys", "pair", vn), with = FALSE]), sys ~ pair, value.var = vn,
             fun.aggregate = function(z) z[1])
  m <- as.matrix(w[, PAIRS, with = FALSE]); rownames(m) <- w$sys; m[GEN_, , drop = FALSE]
}
HI <- grp_mat("kmer_hi"); TO <- grp_mat("turnover")
D[, kmer_z := (kmer_div - mean(kmer_div, na.rm = TRUE)) / sd(kmer_div, na.rm = TRUE), by = pair]
KZ <- grp_mat("kmer_z"); KZ[is.na(KZ)] <- 0
say("matrices built: %d genes x %d cells (%d pairs x %d timepoints)",
    length(GEN_), length(CELLS), length(PAIRS), length(TP))

# ratio per cell, for an arbitrary gene subset (bootstrap resamples rows of these matrices)
arm_cells <- function(X, XN, G, rows = seq_along(GEN_)) {
  out <- rep(NA_real_, length(CELLS))
  for (c_ in seq_along(CELLS)) {
    g <- G[rows, match(pair_of[c_], PAIRS)]
    x <- X[rows, c_]; xn <- XN[rows, c_]
    ih <- which(g == 1L); il <- which(g == 0L)
    if (length(ih) < 5L || length(il) < 5L) next
    vh <- sqrt(max(0, mean(x[ih]^2) - mean(xn[ih]^2)))
    vl <- sqrt(max(0, mean(x[il]^2) - mean(xn[il]^2)))
    if (vl > 0) out[c_] <- vh / vl
  }
  out
}
med_by_t <- function(v) vapply(TP, function(t) median(v[tp_of == t], na.rm = TRUE), 0)
ARM <- list()
for (axis in c("cis", "trans")) {
  X <- if (axis == "cis") CIS else TRN; XN <- if (axis == "cis") CISN else TRNN
  for (grp in c("kmer_hi", "turnover")) {
    G <- if (grp == "kmer_hi") HI else TO
    obs <- arm_cells(X, XN, G)
    bt <- matrix(NA_real_, NBOOT, length(TP))
    for (k in seq_len(NBOOT)) {
      rows <- sample.int(length(GEN_), length(GEN_), TRUE)
      bt[k, ] <- med_by_t(arm_cells(X, XN, G, rows))
    }
    ob <- med_by_t(obs)
    for (j in seq_along(TP)) {
      ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
      ARM[[length(ARM) + 1L]] <- data.table(axis = axis, split = grp, time = TP[j],
        ratio = ob[j], lo = ci[1], hi = ci[2],
        cells_gt1 = sum(obs[tp_of == TP[j]] > 1, na.rm = TRUE),
        cells = sum(!is.na(obs[tp_of == TP[j]])),
        equiv_bound = max(abs(c(ci[1], ci[2]) - 1)))
    }
  }
}
ARM <- rbindlist(ARM)
fwrite(ARM, file.path(outd, "b3_arms_cis_vs_trans.tsv"), sep = "\t")
say("\nratio = rms_tru(axis | high promoter divergence) / rms_tru(axis | low), median over 6 pairs")
say("%-6s %-10s %5s %8s %19s %8s %14s", "axis", "split", "time", "ratio", "95%% CI",
    "cells>1", "equiv |r-1|>")
for (i in seq_len(nrow(ARM))) say("%-6s %-10s %5d %8.3f [%7.3f,%7.3f] %5d/%d %13.3f",
  ARM$axis[i], ARM$split[i], ARM$time[i], ARM$ratio[i], ARM$lo[i], ARM$hi[i],
  ARM$cells_gt1[i], ARM$cells[i], ARM$equiv_bound[i])

say("\n---- THE DISCRIMINATING CONTRAST (cis ratio minus trans ratio, same split) ----")
for (grp in c("kmer_hi", "turnover")) {
  say("split = %s:", grp)
  for (t in TP) {
    rc <- ARM[axis == "cis" & split == grp & time == t, ratio]
    rt <- ARM[axis == "trans" & split == grp & time == t, ratio]
    say("  %2d min: cis %.3f vs trans %.3f -> contrast %+.3f %s", t, rc, rt, rc - rt,
        ifelse(rc > 1 & rt <= 1, "<- pattern Model B predicts",
        ifelse(abs(rc - rt) < 0.05, "<- both move together => confounded, NOT Model B", "")))
  }
}

## ------------- CONTROL 3: amplitude- and expression-matched permutation ---
say("\n================ CONTROL 3: AMPLITUDE- AND EXPRESSION-MATCHED PERMUTATION ================")
say("rms of a divergence axis scales with response amplitude, so a group difference is")
say("confounded if the groups respond differently. Matching removes it without needing the")
say("noise-corrected amplitude that the A6 output does not carry.")
# Per-gene mean response amplitude, then deciles. The permutation relabels which genes are
# called "high promoter divergence" while holding the amplitude-decile composition fixed, so
# any residual ratio cannot be an amplitude effect.
amp_gene <- rowMeans(AMPA, na.rm = TRUE)
amp_dec <- cut(amp_gene, breaks = quantile(amp_gene, seq(0, 1, .1), na.rm = TRUE),
               include.lowest = TRUE, labels = FALSE)
hi_any <- as.integer(rowSums(HI, na.rm = TRUE) >= ceiling(length(PAIRS) / 2))  # majority of pairs
need <- table(factor(amp_dec[hi_any == 1L], levels = 1:10))
have <- table(factor(amp_dec[hi_any == 0L], levels = 1:10))
feasible <- all(as.integer(need) <= as.integer(have))
say("amplitude-decile composition, high vs low promoter-divergence genes (they DO differ,")
say("which is exactly why matching is needed):")
say("  decile : %s", paste(sprintf("%d:%d/%d", 1:10, as.integer(need), as.integer(have)),
                            collapse = " "))
say("  exact decile match feasible: %s", feasible)
MP <- list()
for (axis in c("cis", "trans")) {
  X <- if (axis == "cis") CIS else TRN; XN <- if (axis == "cis") CISN else TRNN
  ob <- med_by_t(arm_cells(X, XN, HI))
  nullm <- matrix(NA_real_, NPERM, length(TP))
  for (k in seq_len(NPERM)) {
    drawn <- unlist(lapply(1:10, function(dd) {
      kk <- as.integer(need)[dd]; pool <- which(amp_dec == dd & hi_any == 0L)
      if (kk == 0L || !length(pool)) integer(0) else sample(pool, min(kk, length(pool)))
    }))
    Gp <- matrix(0L, length(GEN_), length(PAIRS)); Gp[drawn, ] <- 1L
    nullm[k, ] <- med_by_t(arm_cells(X, XN, Gp))
  }
  for (j in seq_along(TP)) {
    ci <- quantile(nullm[, j], c(.025, .975), na.rm = TRUE)
    MP[[length(MP) + 1L]] <- data.table(axis = axis, time = TP[j], feasible = feasible,
      obs = ob[j], null_med = median(nullm[, j], na.rm = TRUE), lo = ci[1], hi = ci[2],
      p_emp = mean(nullm[, j] >= ob[j], na.rm = TRUE))
  }
}
MP <- rbindlist(MP)
fwrite(MP, file.path(outd, "b3_matched_permutation.tsv"), sep = "\t")
say("\n%-6s %5s %9s %9s %19s %8s", "axis", "time", "obs", "null_med", "null 95%%", "p_emp")
for (i in seq_len(nrow(MP))) say("%-6s %5d %9.3f %9.3f [%7.3f,%7.3f] %8.3f",
  MP$axis[i], MP$time[i], MP$obs[i], MP$null_med[i], MP$lo[i], MP$hi[i], MP$p_emp[i])
say("amplitude-decile match feasible: %s", MP$feasible[1])

## ---------------------- CONTROL 2: power calibration (required) ----------
say("\n================ CONTROL 2: POWER CALIBRATION FOR NULL-RESULT INTERPRETATION ================")
say("Inject a coupling of strength c between promoter divergence and cis divergence, homologous")
say("to the statistic's numerator: cis -> cis + c * kmer_z * |amp|, for the high-divergence")
say("genes only. c=0 must NOT come out significant; the smallest detectable c is the resolution")
say("of the promoter HSE test's negative. A null below that threshold is a statement about the test.")
# The injection must be homologous to the statistic's numerator: add a promoter-divergence-
# proportional term to cis for the high-divergence genes, scaled by that gene's own response
# amplitude (so the injected effect is a real multiplicative divergence increment, not a
# constant offset that the rms would treat differently from a biological effect).
inject_cis <- function(cc) {
  M <- CIS
  for (c_ in seq_along(CELLS)) {
    j <- match(pair_of[c_], PAIRS)
    add <- cc * KZ[, j] * AMPA[, c_]
    add[HI[, j] != 1L | is.na(add)] <- 0
    M[, c_] <- M[, c_] + add
  }
  M
}
NB_PWR <- 400L
PWR <- list()
for (cc in c(0, 0.05, 0.10, 0.20, 0.40, 0.80)) {
  Xi <- inject_cis(cc)
  ob <- med_by_t(arm_cells(Xi, CISN, HI))
  bt <- matrix(NA_real_, NB_PWR, length(TP))
  for (k in seq_len(NB_PWR)) {
    rows <- sample.int(length(GEN_), length(GEN_), TRUE)
    bt[k, ] <- med_by_t(arm_cells(Xi, CISN, HI, rows))
  }
  for (j in seq_along(TP)) {
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    PWR[[length(PWR) + 1L]] <- data.table(c_inj = cc, time = TP[j],
      ratio = ob[j], lo = ci[1], hi = ci[2], excludes_1 = ci[1] > 1 | ci[2] < 1)
  }
}
PWR <- rbindlist(PWR)
fwrite(PWR, file.path(outd, "b3_power_calibration.tsv"), sep = "\t")
say("\n%-8s %5s %9s %19s %10s", "c", "time", "ratio", "95%% CI", "excl. 1?")
for (i in seq_len(nrow(PWR))) say("%-8.2f %5d %9.3f [%7.3f,%7.3f] %10s",
  PWR$c_inj[i], PWR$time[i], PWR$ratio[i], PWR$lo[i], PWR$hi[i],
  ifelse(PWR$excludes_1[i], "YES", "no"))
say("\nc=0 false positives (must be 0): %d of %d", PWR[c_inj == 0 & excludes_1, .N],
    PWR[c_inj == 0, .N])
det <- PWR[c_inj > 0 & excludes_1 == TRUE, if (.N) min(c_inj) else NA_real_]
say("smallest detectable coupling c: %s", ifelse(is.na(det), "none up to 0.80", sprintf("%.2f", det)))
if (is.na(det)) {
  say("*** NO DETECTABLE COUPLING AT ANY TESTED STRENGTH -> the promoter HSE test is not informative")
  say("*** for a negative conclusion about HSE architecture (stop condition 2).")
}

## ----------------- strain-mismatch sensitivity, Sk vs the rest ---
say("\n================ STRAIN-MISMATCH SENSITIVITY ================")
say("Only S. kudriavzevii (CR85) has the same strain on the sequence and the expression axis.")
say("Promoters are more polymorphic than coding sequence, so a reference-strain promoter")
say("difference may not exist in the sequenced strain. If the three Sk-containing pairs behave")
say("like the three that do not contain Sk, the mismatch is not driving the result.")
SK <- list()
for (hs in c(0L, 1L)) {
  sel <- which((grepl("Sk", pair_of)) == (hs == 1L))
  v <- arm_cells(CIS, CISN, HI)
  for (t in TP) {
    ii <- sel[tp_of[sel] == t]
    SK[[length(SK) + 1L]] <- data.table(has_Sk = hs, time = t,
      pairs = paste(unique(pair_of[sel]), collapse = ","),
      ratio = median(v[ii], na.rm = TRUE))
  }
}
SK <- rbindlist(SK)
fwrite(SK, file.path(outd, "b3_strain_sensitivity.tsv"), sep = "\t")
say("%-8s %5s %9s  %s", "has_Sk", "time", "cis ratio", "pairs")
for (i in seq_len(nrow(SK))) say("%-8d %5d %9.3f  %s", SK$has_Sk[i], SK$time[i],
  SK$ratio[i], SK$pairs[i])

## --------------------------------------------------- what this test concludes ---
say("\n================ WHAT CAN BE CONCLUDED ================")
say("Interpretation criteria:")
say("  1. scan validated?          %s", SCAN_OK)
say("  2. coupling detectable?     %s", ifelse(is.na(det), "NO -> test not informative",
                                               sprintf("yes, c >= %.2f", det)))
say("  3. cis arm (kmer split), 15/30/60: %.3f / %.3f / %.3f",
    ARM[axis == "cis" & split == "kmer_hi" & time == 15L, ratio],
    ARM[axis == "cis" & split == "kmer_hi" & time == 30L, ratio],
    ARM[axis == "cis" & split == "kmer_hi" & time == 60L, ratio])
say("  4. trans arm (same split):         %.3f / %.3f / %.3f",
    ARM[axis == "trans" & split == "kmer_hi" & time == 15L, ratio],
    ARM[axis == "trans" & split == "kmer_hi" & time == 30L, ratio],
    ARM[axis == "trans" & split == "kmer_hi" & time == 60L, ratio])
say("  5. matched permutation p_emp (cis): %.3f / %.3f / %.3f",
    MP[axis == "cis" & time == 15L, p_emp], MP[axis == "cis" & time == 30L, p_emp],
    MP[axis == "cis" & time == 60L, p_emp])
say("")
say("Interpretation. A cis arm that is flat within the calibrated")
say("resolution indicates that, within the detectable coupling range, target-side HSE and")
say("promoter architecture divergence was not associated with allele-specific divergence of the")
say("heat-stress response; this does not show that HSE architecture is irrelevant or that upstream")
say("regulation evolved, because Model A was not tested directly.")
say("A null the promoter HSE test is scientifically informative: it moves the residual explanation toward broader")
say("trans/upstream control WITHOUT demonstrating it.")
say("Limits of this negative: nucleosome")
say("context is not modelled (no cross-species data, though Pincus 2018 shows it matters);")
say("the promoter window is %d bp and cannot see distal elements; only reference-strain", WIN)
say("promoters were available for 3 of 4 species; the HSE model is the two locked")
say("literature classes and no PWM was fitted.")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_b3_promoter.txt"))
say("\n=== done ===")
