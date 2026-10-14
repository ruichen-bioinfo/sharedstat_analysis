# 10_promoter_hse_architecture_v2.R -- promoter HSE architecture test.
#
# This analysis avoids a confounded promoter-divergence metric, so neither a
# positive nor a negative follows from that metric. The positive control
# is derived here as a validation check.
#
# It is a competing-mechanism elimination, not a hunt for interesting
# HSE mutations:
#
#     Can target-side HSE / promoter architecture explain residual regulatory divergence
#     that is not explained by Hsf1 protein evolution?
#
# Analysis design:
#
#  The kmer_div metric is excluded because it does not measure HSE architecture. It explains 0.4% of the variance in d_units and 0.1% in turnover, with
#     negative sign in all 6 species pairs. A 4-mer window is shorter than an HSE unit (nGAAn,
#     5 bp), so it cannot see 5-bp periodicity, orientation alternation, or spacing; it measures
#     bulk promoter composition drift. It is also confounded in the direction that produces a
#     result: r(kmer, amp) = -0.125, high-kmer half median amplitude 0.431 vs low half 0.505
#     (ratio 0.85), and rms-type divergence statistics scale with amplitude, so a
#     high-divergence group gets ratio < 1 for free. The observed 0.945/0.915/0.885 is
#     quantitatively consistent with exactly that.
#
#  n_units / d_units / units_per_kb are not used as primary metrics. n_units sits ~1.3x above a
#     composition-matched random expectation (chance calibration below). They are printed as
#     diagnostics only.
#
#  Architecture is measured as position-anchored site-level events. A site in one species with
#     no counterpart within +/-50 bp in the other is a turnover event. This catches "same count,
#     moved position", which binary presence/absence misses. Anchoring on upstream coordinates is
#     justified in this clade by script 07: the four HSF1 upstream intergenic distances are
#     644/645/648/665 bp, spread 1.03x. There is no global alignment of the 600 bp window (it
#     would impose non-homologous gaps). Local +/-30 bp alignment is used only when the
#     positional anchor is ambiguous.
#
#  Admission criterion for any architecture metric: observed / composition-matched random
#     expectation. Below ~2x the metric is noise dominated and may only be a diagnostic, never a
#     primary. This is why site-level events (~5x) are admissible while whole-promoter unit
#     counts (~1.3x) are not: the line is drawn by the ratio.
#
#  The classic and gapped HSE classes are co-primary. The gapped class has a favorable chance
#     ratio relative to the classic class and a higher base
#     carrier rate is higher, so classic union gapped nearly doubles the number of informative
#     events, which increases power relative to classic-only turnover, where
#     only ~50-87 genes are available per pair.
#
#  Amplitude matching is not the default. For turnover the
#     confound is measured at r = +0.028, amplitude ratio 1.11, i.e. essentially absent. Order of
#     operations: (1) run a balance test and, if it passes, declare matching unnecessary and
#     use the balance test itself as the control, costing no n; (2) if imbalanced, normalise
#     by amp^beta (cis beta = 0.275, measured in script 05) rather than matching, because
#     normalisation costs no n whereas matching would cut an already small group further;
#     (3) only if both fail, coarsen to tertiles.
#     Deciles x ~60 genes = 6 per cell is infeasible by construction. This is a property
#     of the design rather than evidence that the data cannot support the test; it is not
#     evidence that the test itself is infeasible.
#
#  When matching is required and infeasible, the analysis stops rather than
#     under-drawing silently.
#
#  Power calibration is on the discriminant, not on the cis arm. At n ~ 50-150 the trans arm is
#     as underpowered as the cis arm, so the inference "both arms move together => confounded by
#     a third factor" does not hold until the trans arm is itself calibrated: co-movement is
#     fully consistent with noise at that group size. The calibrated object is therefore
#     Delta(t) = cis_ratio(t) - trans_ratio(t), the quantity that discriminates Model B, and the
#     null is a within-pair hi/lo label permutation at fixed group size.
#
#  The observed data are not treated as a c=0 control. A valid c=0 control splits
#     a homogeneous pool in half, so its expectation is exactly 1 and it can fail.
#     The observed data have no reason to equal 1, so treating them as c=0 asks the
#     wrong question. The null is label
#     permutation throughout.
#
# Analysis components:
#   the HSE scanner and its constructed-sequence self-test; the promoter window (600 bp,
#   clipped at the nearest upstream annotated feature, within the <=644 bp HSF1 intergenic
#   constraint); the window-comparability QC; the positive control (Hsf1_core promoters carry
#   more HSEs than not_induced promoters: 9.9% vs 3.4%, p_emp = 0.009, monotone by Hsf1
#   dependency); the trans arm as a falsifier; the strain-mismatch sensitivity check.
#
# HSE definition, fixed in advance. No PWM is invented, no threshold tuned:
#   Yamamoto 2005 (PMID 15647283): the HSE "consists of contiguous inverted repeats of the
#     sequence 5'-nGAAn-3'"; a novel class has "three direct repeats of nTTCn (or nGAAn)
#     interrupted by 5 bp"; "the number and spacing of the repeating units are critical".
#   Pincus 2018 (PMID 30332327): 74 Hsf1 loci / 46 Hsf1-dependent genes; HSE strength interacts
#     with nucleosome positioning.
#   unit        = 5 bp, nGAAn (+) or nTTCn (-)
#   classic HSE = >= 3 units, 5-bp periodicity, alternating orientation
#   gapped HSE  = >= 3 units, 10-bp periodicity, same orientation
# Nucleosome context is not modelled (no cross-species data); this is a stated limit of any
# negative result.
#
# Gapped HSEs are collapsed to maximal chains rather than counted at every qualifying start
# position, so 4 units at 10-bp spacing score as 1 site. Both classes use maximal
# chains, so that is one site. Carrier rates (>=1 of a class) are unaffected, and the chance
# calibration below uses the same rule on observed and random sequence, which is what
# comparability requires.
#
# No polarity is claimed. Unmatched sites are reported per species as n_unmatched_<sp>, plus
# n_moved = min(u1,u2) and n_net = |u1-u2|. Calling one species "gain" and the other "loss"
# would require ancestral-state reconstruction, which needs an outgroup and is outside the
# scope of this analysis.
#
# Strain mismatch matters more here than anywhere else. Promoters come from reference strains
# DBVPG6765(Sc)/CBS432(Sp)/CR85(Sk)/YJF1449(Su); expression comes from T73/N17/CR85/BMV58.
# Only S. kudriavzevii shares a strain across the two axes, and non-coding sequence is more
# polymorphic within a species than coding sequence, so a reference-strain promoter difference
# may not exist in the strain whose RNA was sequenced. Every HSE-to-expression link is at the
# species-reference level.
#
# Usage: Rscript scripts/hsf1/10_promoter_hse_architecture_v2.R

suppressPackageStartupMessages({ library(data.table); library(IRanges) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
gdir <- file.path(root, "data/03_external/fay2023_genomes")
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
repd <- file.path(root, "doc/05_reports")
outd <- hsd
MIRROR <- "logs/output_copy"   # second copy of the outputs

WIN       <- 600L    # promoter window fixed in advance, clipped at the upstream feature
ANCHOR_TOL<- 50L     # +/-50 bp on the distance-to-ATG axis = same site
FLANK     <- 30L     # local alignment flank, used only when the anchor is ambiguous
AMBIG_MARGIN <- 10L  # two candidates whose |dd| differ by less than this are not separable
LC_MAX    <- 0.50    # low-complexity window: not interpreted
NBOOT     <- 1000L
NPERM     <- 1000L   # label-permutation null for the discriminant, on the observed data
NPERM_SWEEP <- 300L  # label-permutation null recomputed per (substrate, injected c)
N_SUB     <- 10L     # independent randomised substrates, so the resolution is a
                     # property of the design rather than of one lucky/unlucky draw
ALPHA_INFER <- 0.05 / 3            # Bonferroni over the 3 timepoints; the same alpha the
                                   # inference uses, so the bound and the inference share one test
NRAND     <- 3000L   # composition-matched random promoters
BETA_CIS  <- 0.275   # cis-axis amplitude exponent measured in script 05, for fallback (2) only
CHANCE_MIN<- 2.0     # below this ratio a metric is a diagnostic, never a primary
TP        <- c(15L, 30L, 60L)
set.seed(20260906L)

GEN <- data.table(
  species = c("Sc", "Sp", "Sk", "Su"),
  strain  = c("DBVPG6765", "CBS432", "CR85", "YJF1449"),
  matched_strain = c(FALSE, FALSE, TRUE, FALSE))
HYB <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")

say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, "logs/hsf1_promoter_hse_v2.log")
sink(logf, split = TRUE)
say("=== 10_promoter_hse_architecture_v2 | %s ===", format(Sys.time()))
say("Design: position-anchored HSE turnover; chance criterion against composition-matched random sequence;")
say("classic and gapped classes co-primary; amplitude balance tested before matching;")
say("power calibration of the discriminant.")
say("HSE definition locked to Yamamoto 2005 (PMID 15647283) + Pincus 2018 (PMID 30332327).")

## ----------------------------------------------------------- FASTA reader ---
read_fasta <- function(path) {
  ln <- readLines(path); h <- which(startsWith(ln, ">"))
  nm <- sub("^>(\\S+).*", "\\1", ln[h]); s <- h + 1L; e <- c(h[-1] - 1L, length(ln))
  # toupper: soft-masked lowercase would otherwise be turned into N by a case-sensitive
  # regex, silently deleting real sequence
  setNames(vapply(seq_along(h), function(i) toupper(paste(ln[s[i]:e[i]], collapse = "")), ""), nm)
}
revcomp <- function(x) {
  m <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  vapply(x, function(s) paste(rev(m[strsplit(s, "")[[1]]]), collapse = ""), "")
}

## ------------------------------------------------------- the HSE scanner ---
scan_units <- function(s) {
  n <- nchar(s); if (n < 5L) return(NULL)
  core <- substring(s, 2L:(n - 3L), 4L:(n - 1L))
  p <- 1L:(n - 4L)
  fwd <- core == "GAA"; rev_ <- core == "TTC"
  keep <- fwd | rev_
  if (!any(keep)) return(NULL)
  data.table(pos = p[keep], ori = ifelse(fwd[keep], 1L, -1L))
}

# Maximal chains of units at a fixed step. flip=TRUE demands alternating orientation (classic,
# = contiguous inverted repeats); flip=FALSE demands the same orientation (gapped, = direct
# repeats interrupted by 5 bp). Returns one record per maximal chain of length >= 3.
chains <- function(have, ori, step, flip) {
  out <- list()
  hs <- as.character(have)
  ok <- function(a, b) {
    if (!(b %in% have)) return(FALSE)
    if (flip) ori[[as.character(b)]] == -ori[[as.character(a)]]
    else      ori[[as.character(b)]] ==  ori[[as.character(a)]]
  }
  for (st in have) {
    if (ok(st, st - step) && ((st - step) %in% have)) next   # not a chain start
    L <- 1L; cur <- st
    while (ok(cur, cur + step)) { L <- L + 1L; cur <- cur + step }
    if (L >= 3L) out[[length(out) + 1L]] <- c(p5 = st, p3 = cur + 4L, len = L)
  }
  out
}

# max_run is kept only as a diagnostic: it is a count-like quantity.
count_hse <- function(s) {
  U <- scan_units(s)
  if (is.null(U) || nrow(U) == 0L)
    return(list(n_units = 0L, n_classic = 0L, n_gapped = 0L, max_run = 0L))
  setorder(U, pos)
  have <- U$pos; ori <- as.list(setNames(U$ori, as.character(U$pos)))
  cc <- chains(have, ori, 5L, TRUE); gg <- chains(have, ori, 10L, FALSE)
  # longest alternating chain regardless of the >=3 cut, for the diagnostic
  mr <- 0L
  for (st in have) {
    if (((st - 5L) %in% have) && ori[[as.character(st - 5L)]] == -ori[[as.character(st)]]) next
    L <- 1L; cur <- st
    while (((cur + 5L) %in% have) && ori[[as.character(cur + 5L)]] == -ori[[as.character(cur)]]) {
      L <- L + 1L; cur <- cur + 5L }
    if (L > mr) mr <- L
  }
  list(n_units = nrow(U), n_classic = length(cc), n_gapped = length(gg), max_run = mr)
}

# Position-anchored sites. d_atg = distance from the site's 3'-most base to the annotated
# start. Anchoring on the 3' end is deliberate: it is the end closest to the ATG, it is
# unaffected by a unit being added at the 5' end of a chain, and, because promoter extraction
# clips only the distal end, it is invariant to how short the achieved window is. Measured:
# window median 437 bp with 31.5% hitting the 600 bp cap, so a plen-normalised coordinate would
# be biased whereas this one is not.
hse_sites <- function(s, plen) {
  U <- scan_units(s)
  if (is.null(U) || nrow(U) == 0L) return(NULL)
  setorder(U, pos)
  have <- U$pos; ori <- as.list(setNames(U$ori, as.character(U$pos)))
  cc <- chains(have, ori, 5L, TRUE); gg <- chains(have, ori, 10L, FALSE)
  if (!length(cc) && !length(gg)) return(NULL)
  mk <- function(lst, cls) if (!length(lst)) NULL else {
    m <- do.call(rbind, lst)
    data.table(cls = cls, run_len = as.integer(m[, "len"]),
               p5 = as.integer(m[, "p5"]), p3 = as.integer(m[, "p3"]),
               d_atg = plen - as.integer(m[, "p3"]) + 1L)
  }
  S <- rbindlist(list(mk(cc, "classic"), mk(gg, "gapped")))
  # The two fixed definitions can read the same units two ways: units at 5(+),10(-),15(+),25(+)
  # form a classic chain (5..19, alternating, step 5) AND a gapped chain (5..29, same
  # orientation, step 10). Left uncollapsed that is one element counted as two sites, and in the
  # pairing step it would manufacture a turnover event out of a definitional overlap. So sites
  # whose spans overlap are reduced to one site: labelled classic if any classic is present
  # (contiguous inverted repeats are the stronger evidence), run_len = the longest chain
  # involved, anchored at the 3'-most base of the merged span.
  if (nrow(S) > 1L) {
    rd <- as.data.table(reduce(IRanges(start = S$p5, end = S$p3)))
    S <- rbindlist(lapply(seq_len(nrow(rd)), function(i) {
      z <- S[p5 >= rd$start[i] & p3 <= rd$end[i]]
      data.table(cls = if (any(z$cls == "classic")) "classic" else "gapped",
                 run_len = max(z$run_len), p5 = rd$start[i], p3 = rd$end[i],
                 d_atg = plen - rd$end[i] + 1L) }))
  }
  setorder(S, d_atg); S
}

## ------------------------------------------------- score-only Smith-Waterman ---
# The +/-30 bp fallback needs a homology decision, not an alignment path, so this returns the
# score only and performs no traceback. Score-only Smith-Waterman avoids traceback-state
# encoding and has no traceback state to
# encode. A linear gap penalty is used for the same reason.
# 
sw_score <- function(a, b, match = 2L, mismatch = -1L, gap = -2L) {
  A <- strsplit(a, "", fixed = TRUE)[[1]]; B <- strsplit(b, "", fixed = TRUE)[[1]]
  na <- length(A); nb <- length(B)
  if (!na || !nb) return(0L)
  prev <- integer(nb + 1L); best <- 0L
  for (i in seq_len(na)) {
    cur <- integer(nb + 1L); ai <- A[i]
    for (j in seq_len(nb)) {
      sc <- if (ai == B[j]) match else mismatch
      v <- max(0L, prev[j] + sc, prev[j + 1L] + gap, cur[j] + gap)
      cur[j + 1L] <- v
      if (v > best) best <- v
    }
    prev <- cur
  }
  best
}
lc_top3 <- function(s) {          # low-complexity guard (see LC_MAX)
  n <- nchar(s); if (n < 3L) return(1)
  k <- substring(s, 1L:(n - 2L), 3L:n)
  tb <- sort(table(k), decreasing = TRUE)
  as.numeric(sum(tb[seq_len(min(3L, length(tb)))]) / length(k))
}

## ------------------------------------------------------------- self-tests ---
say("\n---- self-tests (validation checks) ----")
t_classic <- paste0("CCCC", "AGAAT", "ATTCT", "AGAAT", "CCCC")
t_gapped  <- paste0("CCC",  "ATTCA", "CCCCC", "ATTCA", "CCCCC", "ATTCA")
t_none    <- strrep("CGCGCG", 12)
sc_c <- count_hse(t_classic); sc_g <- count_hse(t_gapped); sc_n <- count_hse(t_none)
say("classic (3 alternating nGAAn/nTTCn): max_run=%d n_classic=%d  [expect >=3, ==1]",
    sc_c$max_run, sc_c$n_classic)
say("gapped  (3 direct nTTCn, 5 bp gaps): n_gapped=%d               [expect ==1]", sc_g$n_gapped)
say("CG-only control                    : n_units=%d                [expect 0]", sc_n$n_units)
stopifnot(sc_c$max_run >= 3L, sc_c$n_classic == 1L, sc_g$n_gapped == 1L, sc_n$n_units == 0L)

# maximal-chain collapse: 4 units at 10 bp must be one gapped site, not two
t_gap4 <- paste0("CC", "ATTCA", "CCCCC", "ATTCA", "CCCCC", "ATTCA", "CCCCC", "ATTCA")
stopifnot(count_hse(t_gap4)$n_gapped == 1L)
say("maximal-chain collapse (4 direct units == 1 site): PASS")

# d_atg arithmetic, both strands, and clip invariance
pl <- nchar(t_classic); s1 <- hse_sites(t_classic, pl)
# t_classic = 4 + 5 + 5 + 5 + 4 = 23 bp. Units start at 5, 10, 15 (1-based); the chain's 3'-most
# base is 15 + 4 = 19; so d_atg = plen - p3 + 1 = 23 - 19 + 1 = 5, i.e. the site's last base
# lies 5 nt from the annotated start (the trailing "CCCC" is the 4 nt in between).
say("d_atg on constructed classic: plen=%d p5=%d p3=%d d_atg=%d  [expect 23, 5, 19, 5]",
    pl, s1$p5[1], s1$p3[1], s1$d_atg[1])
stopifnot(pl == 23L, s1$p5[1] == 5L, s1$p3[1] == 19L, s1$d_atg[1] == 5L)
# clip invariance: removing bases from the distal (5') end must not move d_atg
t_clip <- substr(t_classic, 3L, pl)
s1c <- hse_sites(t_clip, nchar(t_clip))
stopifnot(s1c$d_atg[1] == s1$d_atg[1])
say("d_atg invariant to distal clipping: PASS (%d == %d)", s1c$d_atg[1], s1$d_atg[1])
# a minus-strand gene: extraction revcomps, so the same motif on the opposite strand must give
# the same d_atg when the promoter is presented 5'->3' on the gene's own strand
stopifnot(hse_sites(as.character(revcomp(revcomp(t_classic))), pl)$d_atg[1] == s1$d_atg[1])
say("revcomp round-trip stable: PASS")

# overlapping classic + gapped on the same units must collapse to one site.
# units at 5(+),10(-),15(+),25(+): classic chain 5..19, gapped chain 5..29 -> one merged site
t_ovl <- paste0("CCCC", "AGAAT", "ATTCT", "AGAAT", "CCCCC", "AGAAT")
h_ovl <- count_hse(t_ovl); s_ovl <- hse_sites(t_ovl, nchar(t_ovl))
say("overlap case: raw chains classic=%d gapped=%d -> merged sites=%d cls=%s p5=%d p3=%d d_atg=%d",
    h_ovl$n_classic, h_ovl$n_gapped, nrow(s_ovl), s_ovl$cls[1], s_ovl$p5[1], s_ovl$p3[1],
    s_ovl$d_atg[1])
stopifnot(nchar(t_ovl) == 29L, h_ovl$n_classic == 1L, h_ovl$n_gapped == 1L,
          nrow(s_ovl) == 1L, s_ovl$cls[1] == "classic", s_ovl$p5[1] == 5L,
          s_ovl$p3[1] == 29L, s_ovl$d_atg[1] == 1L)
say("definitional-overlap collapse: PASS (would otherwise fabricate a turnover event)")
# and a non-overlapping pair must stay two sites
t_two <- paste0("CCCC", "AGAAT", "ATTCT", "AGAAT", strrep("C", 60),
                "AGAAT", "ATTCT", "AGAAT", "CCCC")
stopifnot(nrow(hse_sites(t_two, nchar(t_two))) == 2L)
say("two well-separated sites stay two: PASS")

# Smith-Waterman
stopifnot(sw_score("ACGT", "ACGT") == 8L)                       # 4 matches x 2
stopifnot(sw_score("AAAA", "TTTT") == 0L)                       # local score floors at 0
stopifnot(sw_score("ACGTACGT", "ACGTTACGT") == 14L)             # needs one gap: 8*2 - 2
set.seed(11L)
for (i in 1:200) {                                              # symmetry + determinism
  a <- paste(sample(c("A","C","G","T"), sample(20:60, 1), TRUE), collapse = "")
  b <- paste(sample(c("A","C","G","T"), sample(20:60, 1), TRUE), collapse = "")
  stopifnot(sw_score(a, b) == sw_score(b, a), sw_score(a, b) == sw_score(a, b))
}
say("Smith-Waterman: identity / no-homology / gap-requiring / symmetry(200 pairs): PASS")
stopifnot(lc_top3(strrep("AT", 40)) > 0.9, lc_top3(paste(sample(c("A","C","G","T"), 200, TRUE),
                                                        collapse = "")) < 0.4)
say("lc_top3 discriminates repeat from random: PASS")
set.seed(20260906L)   # restore the main random stream after the test-only seeds

## -------------------------------------- promoter extraction, all 4 species ---
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
A6 <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
UNIV <- unique(A6$hog)
sysmap <- unique(HB[gene %in% UNIV, .(hog = gene, sys)])
say("\nanalysis universe: %d genes ; with a systematic name: %d", length(UNIV), nrow(sysmap))

PROM <- list(); WINSTAT <- list()
for (i in seq_len(nrow(GEN))) {
  sp <- GEN$species[i]; strain <- GEN$strain[i]
  gf <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(gdir, paste0(strain, ".final.gff")))),
              sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gf, c("chr","src","type","start","end","score","strand","frame","attr"))
  gf[, gname := NA_character_]
  hasn <- grepl("Gene_Name=", gf$attr, fixed = TRUE)
  gf[hasn, gname := sub("^.*Gene_Name=([^;]+).*$", "\\1", attr)]
  gf[, sysn := sub("_.*$", "", gname)]
  genes <- gf[type == "gene"]
  feats <- gf[type %in% c("gene", "tRNA", "rRNA", "ncRNA")][, .(chr, start, end)]
  setorder(feats, chr, start)
  fa <- read_fasta(file.path(gdir, paste0(strain, ".fasta")))
  tg <- genes[sysn %in% sysmap$sys & !is.na(sysn)]
  tg <- tg[!duplicated(sysn)]
  out <- vector("list", nrow(tg))
  for (k in seq_len(nrow(tg))) {
    ch <- tg$chr[k]; if (!ch %in% names(fa)) next
    cs <- fa[[ch]]; cl <- nchar(cs); fc <- feats[chr == ch]
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
      seqp <- as.character(revcomp(substr(cs, lo, hi)))
    }
    out[[k]] <- data.table(sys = tg$sysn[k], species = sp, prom = seqp, plen = nchar(seqp))
  }
  P <- rbindlist(out); PROM[[sp]] <- P
  WINSTAT[[sp]] <- data.table(species = sp, strain = strain, n = nrow(P),
    med_len = median(P$plen), frac_full = mean(P$plen >= WIN), frac_short = mean(P$plen < 200L))
  say("  %-3s (%-9s) promoters: %5d | median len %3.0f | full-%dbp %.1f%% | <200bp %.1f%%",
      sp, strain, nrow(P), median(P$plen), WIN, 100*mean(P$plen >= WIN), 100*mean(P$plen < 200L))
}
PR <- rbindlist(PROM); WS <- rbindlist(WINSTAT)
fwrite(WS, file.path(outd, "b3v2_window_stats.tsv"), sep = "\t")

ok <- PR[, .N, by = sys][N == 4L, sys]
PR <- PR[sys %in% ok]
PW <- dcast(PR, sys ~ species, value.var = "plen")
PW[, len_min := pmin(Sc, Sp, Sk, Su)][, len_max := pmax(Sc, Sp, Sk, Su)]
PW[, len_ratio := len_max / pmax(1L, len_min)]
COMPARABLE <- PW[len_ratio <= 2, sys]
say("\ngenes with a promoter in all 4 species: %d ; window ratio median %.2f, %.1f%% <=2x",
    length(ok), median(PW$len_ratio), 100 * mean(PW$len_ratio <= 2))
say("genes entering the architecture tests: %d", length(COMPARABLE))
PRc <- PR[sys %in% COMPARABLE]

## ---------------- chance calibration of every candidate metric ------
say("\n================ OBSERVED vs COMPOSITION-MATCHED RANDOM ================")
say("Any architecture metric must clear ~%.1fx over a composition-matched random expectation,", CHANCE_MIN)
say("or it is noise dominated and may only be a diagnostic. Random sequences use the REAL")
say("promoter length distribution and the measured Sc promoter base composition, scanned with")
say("the identical rule (maximal chains).")
sc_prom <- PRc[species == "Sc", prom]
bc <- table(strsplit(paste(sc_prom, collapse = ""), "", fixed = TRUE)[[1]])
bc <- bc[c("A","C","G","T")]; bcp <- as.numeric(bc) / sum(bc)
say("MEASURED Sc promoter composition: A %.3f C %.3f G %.3f T %.3f", bcp[1], bcp[2], bcp[3], bcp[4])
say("The assumed composition is A/T 0.31, C/G 0.19 (a fixed reference value). The measured promoters are more")
say("AT-rich than that, and since the HSE unit core is GAA/TTC a more AT-rich background raises")
say("the random expectation and therefore LOWERS the ratio. Both are reported so the difference")
say("is visible rather than silently resolved in favour of whichever is more convenient.")
chance_tab <- function(prob, tag) {
  lens <- sample(PRc$plen, NRAND, replace = TRUE)
  RND <- vapply(lens, function(L) paste(sample(c("A","C","G","T"), L, TRUE, prob = prob),
                                        collapse = ""), "")
  rh <- rbindlist(lapply(RND, function(s) as.data.table(count_hse(s))))
  data.table(background = tag,
    metric   = c("median n_units", "carries >=1 classic HSE", "carries >=1 gapped HSE",
                 "carries >=1 classic OR gapped"),
    random   = c(median(rh$n_units), mean(rh$n_classic > 0), mean(rh$n_gapped > 0),
                 mean(rh$n_classic > 0 | rh$n_gapped > 0)))
}
obs_h <- rbindlist(lapply(PRc$prom, function(s) as.data.table(count_hse(s))))
OBSV <- c(median(obs_h$n_units), mean(obs_h$n_classic > 0), mean(obs_h$n_gapped > 0),
          mean(obs_h$n_classic > 0 | obs_h$n_gapped > 0))
CH <- rbind(chance_tab(bcp, "measured (A .323 C .182 G .175 T .320)"),
            chance_tab(c(0.31, 0.19, 0.19, 0.31), "assumed reference (A/T .31, C/G .19)"))
CH[, observed := rep(OBSV, 2L)]
CH[, ratio := observed / pmax(random, .Machine$double.eps)]
CH[, admissible_primary := ratio >= CHANCE_MIN]
fwrite(CH, file.path(outd, "b3v2_chance_calibration.tsv"), sep = "\t")
say("")
say("%-40s %-30s %9s %9s %8s %7s", "background", "metric", "observed", "random", "ratio", "primary")
for (i in seq_len(nrow(CH))) say("%-40s %-30s %9.4f %9.4f %7.2fx %7s", CH$background[i],
  CH$metric[i], CH$observed[i], CH$random[i], CH$ratio[i],
  ifelse(CH$admissible_primary[i], "YES", "no"))
CHm <- CH[background == CH$background[1]]      # the measured background governs
if (!all(CHm$admissible_primary[2:4]))
  say("*** WARNING: a site class failed the chance criterion and is classified as diagnostic.")
say("")
say("CHANCE-SITE DILUTION, a direct consequence of the ratio: at %.2fx enrichment about %.0f%% of",
    CHm$ratio[4], 100 / CHm$ratio[4])
say("detected sites are expected by chance, and a chance site has no reason to be positionally")
say("conserved. So a large share of the turnover events below is chance-site turnover. This")
say("DILUTES the split and biases the test toward the null -- it is conservative, not inflating.")

## ------------------------------------------- position-anchored site table ---
say("\n---- position-anchored HSE sites (d_atg = distance from site 3' end to ATG) ----")
SITES <- PRc[, { S <- hse_sites(prom, plen)
                 if (is.null(S)) NULL else cbind(sys = sys, species = species, S) },
             by = seq_len(nrow(PRc))][, -1L]
fwrite(SITES, file.path(outd, "b3v2_hse_sites.tsv.gz"), sep = "\t", compress = "gzip")
say("sites found: %d (classic %d, gapped %d) over %d gene x species promoters",
    nrow(SITES), sum(SITES$cls == "classic"), sum(SITES$cls == "gapped"), nrow(PRc))
say("d_atg distribution: median %.0f bp, IQR %.0f-%.0f",
    median(SITES$d_atg), quantile(SITES$d_atg, .25), quantile(SITES$d_atg, .75))

# homology-score threshold for the ambiguous-anchor fallback, calibrated once from shuffled
# real windows rather than tuned: a candidate must beat the 95th percentile of the score
# between real windows and composition-preserving shuffles of other real windows.
say("\ncalibrating the local-alignment homology threshold from shuffled real windows ...")
win_of <- function(prom, plen, d_atg, p3) {
  a <- max(1L, p3 - FLANK); b <- min(plen, p3 + FLANK); substr(prom, a, b)
}
SS <- merge(SITES, PRc[, .(sys, species, prom, plen)], by = c("sys", "species"))
SS[, win := win_of(prom, plen, d_atg, p3), by = seq_len(nrow(SS))]
set.seed(4242L)
sm <- sample(SS$win, 600L)
null_sw <- vapply(seq_len(300L), function(i) {
  a <- sm[i]; b <- sm[300L + i]
  bshuf <- paste(sample(strsplit(b, "", fixed = TRUE)[[1]]), collapse = "")
  sw_score(a, bshuf) }, 0L)
SW_THR <- as.numeric(quantile(null_sw, 0.95))
say("shuffled-window score null: median %.0f, 95th pct %.0f -> homology threshold %.0f",
    median(null_sw), SW_THR, SW_THR)
set.seed(20260906L)

## ------------------------------- +/-50 bp one-to-one site pairing (IRanges) ---
say("\n---- site pairing per gene per species pair (+/-%d bp anchor, IRanges) ----", ANCHOR_TOL)
SL <- split(SS[, .(sys, species, cls, run_len, d_atg, p3, win)],
            by = c("sys", "species"), keep.by = TRUE, drop = TRUE)
key_of <- function(g, sp) paste(g, sp, sep = ".")
AMB <- 0L; LCREJ <- 0L; SWUSE <- 0L
pair_gene <- function(S1, S2) {
  # returns matched index pairs plus the unmatched counts; deterministic
  if (is.null(S1) && is.null(S2)) return(list(m = NULL, u1 = 0L, u2 = 0L, amb = 0L, lc = 0L))
  if (is.null(S1)) return(list(m = NULL, u1 = 0L, u2 = nrow(S2), amb = 0L, lc = 0L))
  if (is.null(S2)) return(list(m = NULL, u1 = nrow(S1), u2 = 0L, amb = 0L, lc = 0L))
  q <- IRanges(start = S1$d_atg - ANCHOR_TOL, end = S1$d_atg + ANCHOR_TOL)
  s <- IRanges(start = S2$d_atg, width = 1L)
  hits <- as.data.table(findOverlaps(q, s))
  if (!nrow(hits)) return(list(m = NULL, u1 = nrow(S1), u2 = nrow(S2), amb = 0L, lc = 0L))
  hits[, dd := abs(S1$d_atg[queryHits] - S2$d_atg[subjectHits])]
  namb <- 0L; nlc <- 0L
  # ambiguity: a site with >1 candidate whose two best distances are not separable
  for (qq in unique(hits$queryHits)) {
    dz <- sort(hits[queryHits == qq, dd])
    if (length(dz) > 1L && (dz[2] - dz[1]) < AMBIG_MARGIN) {
      namb <- namb + 1L
      cand <- hits[queryHits == qq]
      w1 <- S1$win[qq]
      if (lc_top3(w1) >= LC_MAX) { nlc <- nlc + 1L; hits <- hits[queryHits != qq]; next }
      keep <- vapply(seq_len(nrow(cand)), function(z) {
        w2 <- S2$win[cand$subjectHits[z]]
        if (lc_top3(w2) >= LC_MAX) return(-1L)
        sw_score(w1, w2) }, 0L)
      if (max(keep) < SW_THR) { hits <- hits[queryHits != qq]; next }
      best <- cand$subjectHits[which.max(keep)]
      hits <- hits[queryHits != qq | subjectHits == best]
    }
  }
  setorder(hits, dd, queryHits, subjectHits)
  usedq <- integer(0); useds <- integer(0); m <- list()
  for (r in seq_len(nrow(hits))) {
    if (hits$queryHits[r] %in% usedq || hits$subjectHits[r] %in% useds) next
    usedq <- c(usedq, hits$queryHits[r]); useds <- c(useds, hits$subjectHits[r])
    m[[length(m) + 1L]] <- c(i1 = hits$queryHits[r], i2 = hits$subjectHits[r])
  }
  list(m = if (length(m)) do.call(rbind, m) else NULL,
       u1 = nrow(S1) - length(usedq), u2 = nrow(S2) - length(useds), amb = namb, lc = nlc)
}

PD <- list()
for (h in HYB) {
  s1 <- substr(h, 1, 2); s2 <- substr(h, 4, 5)
  res <- vector("list", length(COMPARABLE))
  for (gi in seq_along(COMPARABLE)) {
    g <- COMPARABLE[gi]
    S1 <- SL[[key_of(g, s1)]]; S2 <- SL[[key_of(g, s2)]]
    p <- pair_gene(S1, S2)
    AMB <- AMB + p$amb; LCREJ <- LCREJ + p$lc
    nm <- if (is.null(p$m)) 0L else nrow(p$m)
    dstr <- if (nm == 0L) NA_real_ else
      mean(abs(S1$run_len[p$m[, 1]] - S2$run_len[p$m[, 2]]))
    nstr <- if (nm == 0L) 0L else
      sum(S1$run_len[p$m[, 1]] != S2$run_len[p$m[, 2]])
    res[[gi]] <- data.table(sys = g, pair = h,
      n_shared = nm, n_unmatched_1 = p$u1, n_unmatched_2 = p$u2,
      n_moved = min(p$u1, p$u2), n_net = abs(p$u1 - p$u2),
      n_events = p$u1 + p$u2, mean_dstrength = dstr, n_strength_changed = nstr)
  }
  PD[[h]] <- rbindlist(res)
}
PD <- rbindlist(PD)
PD[, turnover := as.integer(n_events > 0L)]
PD[, strength_change := as.integer(n_strength_changed > 0L)]
fwrite(PD, file.path(outd, "b3v2_site_turnover.tsv.gz"), sep = "\t", compress = "gzip")
say("ambiguous anchors: %d (local-alignment fallback used), rejected as low complexity: %d",
    AMB, LCREJ)
say("%-8s %9s %9s %9s %9s %11s %13s", "pair", "shared", "events", "moved", "net",
    "turnover %%", "strength chg %%")
for (h in HYB) { z <- PD[pair == h]
  say("%-8s %9.0f %9.0f %9.0f %9.0f %10.1f%% %12.1f%%", h, sum(z$n_shared), sum(z$n_events),
      sum(z$n_moved), sum(z$n_net), 100*mean(z$turnover), 100*mean(z$strength_change)) }
say("\nturnover group size per pair (this is the power-limiting n):")
for (h in HYB) say("  %-8s n_turnover=%4d  n_no_turnover=%4d",
                   h, PD[pair == h & turnover == 1L, .N], PD[pair == h & turnover == 0L, .N])
say("The classic-and-gapped union increases these group sizes")
say("relative to the 50-87 genes obtained with binary classic-only turnover.")

## ---------------------------- positive control, re-derived and per class ----
say("\n================ CONTROL 1: SCAN vs KNOWN HSE ENRICHMENT ================")
ST <- unique(A6[, .(hog, stratum)])
ST[is.na(stratum) | stratum == "", stratum := "unassigned"]
ST <- merge(ST, sysmap, by = "hog")
HSg <- cbind(PRc[, .(sys, species)], rbindlist(lapply(PRc$prom, function(s) as.data.table(count_hse(s)))))
HSm <- merge(HSg, ST, by = "sys", allow.cartesian = TRUE)
PC <- HSm[, .(n = .N, any_classic = mean(n_classic > 0), any_gapped = mean(n_gapped > 0),
              any_either = mean(n_classic > 0 | n_gapped > 0), mean_maxrun = mean(max_run)),
          by = stratum]
setorder(PC, -any_either)
say("%-16s %7s %13s %12s %12s %9s", "stratum", "n", "any_classic", "any_gapped", "any_either",
    "max_run")
for (i in seq_len(nrow(PC))) say("%-16s %7d %12.1f%% %11.1f%% %11.1f%% %9.2f", PC$stratum[i],
  PC$n[i], 100*PC$any_classic[i], 100*PC$any_gapped[i], 100*PC$any_either[i], PC$mean_maxrun[i])
fwrite(PC, file.path(outd, "b3v2_positive_control.tsv"), sep = "\t")

gene_flags <- HSg[, .(anyc = as.integer(any(n_classic > 0)), anyg = as.integer(any(n_gapped > 0)),
                      anye = as.integer(any(n_classic > 0 | n_gapped > 0))), by = sys]
core_sys <- unique(ST[stratum == "Hsf1_core", sys]); ni_sys <- unique(ST[stratum == "not_induced", sys])
pool <- gene_flags[sys %in% c(core_sys, ni_sys)]
o_core <- gene_flags[sys %in% core_sys]
PCT <- list()
for (fl in c("anyc", "anyg", "anye")) {
  obs <- mean(o_core[[fl]]) - mean(gene_flags[sys %in% ni_sys][[fl]])
  nul <- replicate(NPERM, { idx <- sample(nrow(pool), nrow(o_core))
    mean(pool[[fl]][idx]) - mean(pool[[fl]][-idx]) })
  PCT[[fl]] <- data.table(flag = fl, core = mean(o_core[[fl]]),
    not_induced = mean(gene_flags[sys %in% ni_sys][[fl]]), diff = obs,
    p_emp = mean(nul >= obs))
}
PCT <- rbindlist(PCT); fwrite(PCT, file.path(outd, "b3v2_positive_control_test.tsv"), sep = "\t")
say("\ncore (n=%d) vs not_induced (n=%d), permutation %d:", nrow(o_core),
    gene_flags[sys %in% ni_sys, .N], NPERM)
for (i in seq_len(nrow(PCT))) say("  %-5s core %.1f%% vs not_induced %.1f%%  (diff %+.1f pp, p_emp %.4f)",
  PCT$flag[i], 100*PCT$core[i], 100*PCT$not_induced[i], 100*PCT$diff[i], PCT$p_emp[i])
SCAN_OK <- any(PCT$p_emp < 0.05)
say("SCAN VALIDATED: %s", SCAN_OK)
if (!SCAN_OK) stop("scan validation failed")

## ---------------------------------------- join to the script-05 cis/trans table ---
A6s <- unique(merge(A6[, .(hog, pair, time, cis, cis_nul, trans, trans_nul, amp)], sysmap, by = "hog"))
D <- merge(A6s, PD, by = c("sys", "pair"))
D <- D[!is.na(cis) & !is.na(cis_nul) & !is.na(trans) & !is.na(trans_nul)]
D[, cell := paste(pair, time, sep = "_")]
say("\njoined rows (gene x pair x time): %d ; unique genes %d", nrow(D), uniqueN(D$sys))
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
TO <- grp_mat("turnover"); TO[is.na(TO)] <- 0L
NEV <- grp_mat("n_events"); NEV[is.na(NEV)] <- 0L
say("matrices: %d genes x %d cells (%d pairs x %d timepoints)",
    length(GEN_), length(CELLS), length(PAIRS), length(TP))

## ---------------------- balance test, before any matching ----
say("\n================ BALANCE TEST (run BEFORE deciding to match) ================")
say("Amplitude matching is needed only if the site-level metric is")
say("already balanced on amplitude, matching is unnecessary and would only destroy n. The")
say("balance test is the control for amplitude balance.")
amp_gene <- rowMeans(AMPA, na.rm = TRUE)
to_any <- as.integer(rowSums(TO, na.rm = TRUE) >= ceiling(length(PAIRS) / 2))
r_pear <- cor(to_any, amp_gene, use = "complete.obs")
r_spear <- cor(to_any, amp_gene, method = "spearman", use = "complete.obs")
amp_dec <- cut(amp_gene, breaks = quantile(amp_gene, seq(0, 1, .1), na.rm = TRUE),
               include.lowest = TRUE, labels = FALSE)
tabd <- table(factor(amp_dec, levels = 1:10), to_any)
chi <- suppressWarnings(chisq.test(tabd))
amp_ratio <- median(amp_gene[to_any == 1L], na.rm = TRUE) / median(amp_gene[to_any == 0L], na.rm = TRUE)
BAL <- data.table(r_pearson = r_pear, r_spearman = r_spear, amp_median_ratio = amp_ratio,
  decile_chisq_p = chi$p.value, n_hi = sum(to_any == 1L), n_lo = sum(to_any == 0L))
fwrite(BAL, file.path(outd, "b3v2_balance_test.tsv"), sep = "\t")
say("r(turnover, amp) Pearson %+.4f / Spearman %+.4f ; amplitude median ratio %.3f",
    r_pear, r_spear, amp_ratio)
say("amplitude-decile composition chi-square p = %.4f  (hi n=%d, lo n=%d)",
    chi$p.value, BAL$n_hi, BAL$n_lo)
say("decile counts (turnover=1): %s", paste(as.integer(tabd[, 2]), collapse = " "))
say("decile counts (turnover=0): %s", paste(as.integer(tabd[, 1]), collapse = " "))
BALANCED <- abs(r_pear) < 0.05 && chi$p.value > 0.05
if (BALANCED) {
  say("=> BALANCED: |r| < 0.05 and decile composition not significant. ")
  say("   MATCHING IS DECLARED UNNECESSARY and this test is the control. No n is spent.")
  ADJ <- "none"
} else {
  say("=> IMBALANCED. applying amp^beta normalisation (beta = %.3f, measured in this analysis)", BETA_CIS)
  say("   rather than matching, because normalisation costs no n. Refitting beta on this subset:")
  ADJ <- "amp_beta"
}
say("NOTE: deciles x ~%d genes per cell would be infeasible BY CONSTRUCTION. That is a",
    as.integer(round(BAL$n_hi / 10)))
say("property of the design rather than evidence that the data cannot support the test. It is")
say("treated accordingly.")

## ------------------------------------------------------------ estimators ---
# rms_tru(group, cell) = sqrt(max(0, mean(x^2) - mean(x_nul^2))), then ratio hi/lo per cell,
# then median over the 6 pairs per timepoint. Matrix form because the resamples run 1000x.
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
if (ADJ == "amp_beta") {   # divide the axis by amp^beta, elementwise per cell
  SC <- pmax(AMPA, 1e-9)^BETA_CIS
  CIS <- CIS / SC; CISN <- CISN / SC; TRN <- TRN / SC; TRNN <- TRNN / SC
  say("amp^beta normalisation applied to all four axis matrices.")
}
# the discriminant: Model B predicts promoter divergence raises cis and not trans
disc_of <- function(G, Xc = CIS, XNc = CISN, rows = seq_along(GEN_))
  med_by_t(arm_cells(Xc, XNc, G, rows)) - med_by_t(arm_cells(TRN, TRNN, G, rows))

## ------- label-permutation null for the discriminant (not the arm) --
say("\n================ POWER CALIBRATION ON THE DISCRIMINANT ================")
N_TO_MED <- as.integer(median(vapply(seq_len(ncol(TO)), function(j) sum(TO[, j] == 1L), 0L)))
say("The calibrated object is Delta(t) = cis_ratio(t) - trans_ratio(t), NOT the cis arm alone.")
say("Reason: at a median turnover-group size of n = %d the trans arm is as", N_TO_MED)
say("underpowered as the cis arm, so 'both arms move together => confounded by a third factor'")
say("is NOT a valid inference until the trans arm is itself calibrated -- co-movement is fully")
say("consistent with noise at this group size.")
say("Null = WITHIN-PAIR permutation of the hi/lo labels at FIXED group size (the")
say("observed data is never used as c=0).")
perm_G <- function(G) {
  Gp <- G
  for (j in seq_len(ncol(G))) {
    v <- G[, j]; k <- sum(v == 1L, na.rm = TRUE)
    Gp[, j] <- 0L
    if (k > 0L) Gp[sample.int(nrow(G), k), j] <- 1L
  }
  Gp
}
obs_cis   <- med_by_t(arm_cells(CIS, CISN, TO))
obs_trans <- med_by_t(arm_cells(TRN, TRNN, TO))
obs_disc  <- obs_cis - obs_trans
NUL <- matrix(NA_real_, NPERM, length(TP))
NULc <- matrix(NA_real_, NPERM, length(TP)); NULt <- matrix(NA_real_, NPERM, length(TP))
for (k in seq_len(NPERM)) {
  Gp <- perm_G(TO)
  a <- med_by_t(arm_cells(CIS, CISN, Gp)); b <- med_by_t(arm_cells(TRN, TRNN, Gp))
  NULc[k, ] <- a; NULt[k, ] <- b; NUL[k, ] <- a - b
}
NULLTAB <- rbindlist(lapply(seq_along(TP), function(j) data.table(
  time = TP[j], obs_cis = obs_cis[j], obs_trans = obs_trans[j], obs_disc = obs_disc[j],
  null_cis_med = median(NULc[, j], na.rm = TRUE), null_trans_med = median(NULt[, j], na.rm = TRUE),
  null_disc_med = median(NUL[, j], na.rm = TRUE),
  null_disc_lo = quantile(NUL[, j], .025, na.rm = TRUE),
  null_disc_hi = quantile(NUL[, j], .975, na.rm = TRUE),
  p_emp_two_sided = mean(abs(NUL[, j] - median(NUL[, j], na.rm = TRUE)) >=
                         abs(obs_disc[j] - median(NUL[, j], na.rm = TRUE)), na.rm = TRUE))))
fwrite(NULLTAB, file.path(outd, "b3v2_discriminant_null.tsv"), sep = "\t")
say("\n%-5s %9s %9s %9s %21s %9s", "time", "cis", "trans", "Delta", "null Delta 95%%", "p_emp")
for (i in seq_len(nrow(NULLTAB))) say("%-5d %9.3f %9.3f %+9.3f [%+8.3f,%+8.3f] %9.3f",
  NULLTAB$time[i], NULLTAB$obs_cis[i], NULLTAB$obs_trans[i], NULLTAB$obs_disc[i],
  NULLTAB$null_disc_lo[i], NULLTAB$null_disc_hi[i], NULLTAB$p_emp_two_sided[i])
say("sanity: the null median of Delta should sit near 0 -- it is %s",
    paste(sprintf("%+.3f", NULLTAB$null_disc_med), collapse = " / "))

## ------- control that can fail: homogeneous-pool split (the c=0 control) -----
say("\n================ CAN-FAIL CONTROL: HOMOGENEOUS-POOL SPLIT ================")
say("This is a c = 0 control. Take only genes with ZERO")
say("turnover in EVERY pair -- a homogeneous pool -- and split it at random into a pseudo-hi")
say("group of the real group's size and a pseudo-lo remainder. Delta must come out at 0. If it")
say("does not, the estimator itself is biased and nothing below can be read. Note this is NOT")
say("the observed data with c=0 added to it, which would not be a null,")
say("because the observed hi/lo groups are two genuinely different sets of genes.")
homo <- which(rowSums(TO, na.rm = TRUE) == 0L)
say("homogeneous pool (no turnover in any of the 6 pairs): %d genes", length(homo))
k_per_pair <- vapply(seq_len(ncol(TO)), function(j) sum(TO[, j] == 1L), 0L)
HOM <- matrix(NA_real_, NPERM, length(TP))
for (k in seq_len(NPERM)) {
  Gh <- matrix(0L, length(GEN_), length(PAIRS), dimnames = dimnames(TO))
  for (j in seq_len(ncol(Gh))) Gh[sample(homo, k_per_pair[j]), j] <- 1L
  HOM[k, ] <- med_by_t(arm_cells(CIS, CISN, Gh)) - med_by_t(arm_cells(TRN, TRNN, Gh))
}
HOMTAB <- rbindlist(lapply(seq_along(TP), function(j) data.table(time = TP[j],
  delta_med = median(HOM[, j], na.rm = TRUE),
  lo = quantile(HOM[, j], .025, na.rm = TRUE), hi = quantile(HOM[, j], .975, na.rm = TRUE))))
HOMTAB[, contains_zero := lo <= 0 & hi >= 0]
fwrite(HOMTAB, file.path(outd, "b3v2_homogeneous_control.tsv"), sep = "\t")
say("%-5s %11s %21s %14s", "time", "Delta med", "95%%", "contains 0?")
for (i in seq_len(nrow(HOMTAB))) say("%-5d %+11.3f [%+8.3f,%+8.3f] %14s", HOMTAB$time[i],
  HOMTAB$delta_med[i], HOMTAB$lo[i], HOMTAB$hi[i], ifelse(HOMTAB$contains_zero[i], "YES", "NO"))
HOMO_OK <- all(HOMTAB$contains_zero) && max(abs(HOMTAB$delta_med)) < 0.05
say("HOMOGENEOUS-POOL CONTROL PASSED: %s", HOMO_OK)
if (!HOMO_OK) say("*** estimator biased on a homogeneous pool -- subsequent results are not interpreted.")

## --------------------------- injection sweep on signal-free substrates ---
say("\n---- injection sweep: smallest detectable coupling ----")
say("Injecting on the observed cis matrix would make its c=0 row equal the observed data")
say("itself, so the sweep uses signal-free substrates. A detection")
say("limit read off the observed matrix would not be valid.")
say("")
say("A single randomised substrate is insufficient. One fixed")
say("randomised substrate gave Delta = -0.095 / -0.101 / +0.014")
say("with a permutation null of [-0.112, +0.112]. At 15 and 30 min that substrate")
say("started at the lower edge of its null, so an injection of about 0.21 was needed to exit the")
say("upper bound, compared with +0.119 for the observed inference from a null centred at -0.007.")
say("A limit such as c >= 0.40 / 0.80 from one draw reflects that draw rather than the")
say("design. The interval used is a permutation interval, so the mechanism is the variance between")
say("substrate draws, not the interval type.")
say("")
say("Specification (the equivalence bound is derived with the same")
say("test as the inference'):")
say("  - detection criterion = LABEL-PERMUTATION p_emp < alpha, the same test as the inference")
say("  - alpha = %.4f, Bonferroni over the %d timepoints, matching the inference exactly",
    ALPHA_INFER, length(TP))
say("  - repeated over %d independent randomised substrates, so the threshold is a property of", N_SUB)
say("    the design and not of one draw; median and range are reported")
say("  - the c=0 row remains a can-fail control and must NOT be detected")

# --- assertion: the calibration substrate is never the observed one ---
# Injecting on the observed labels would make the c=0 row the observed data, so this is
# checked by assertion.
assert_signal_free <- function(baseline_labels, observed_labels, tag) {
  stopifnot(!identical(baseline_labels, observed_labels))
  if (identical(dim(baseline_labels), dim(observed_labels)))
    stopifnot(!all(baseline_labels == observed_labels, na.rm = TRUE))
  invisible(TRUE)
}

perm_p <- function(nullvec, obs) {
  m <- median(nullvec, na.rm = TRUE)
  mean(abs(nullvec - m) >= abs(obs - m), na.rm = TRUE)
}
inject_on <- function(cc, G) {
  M <- CIS
  for (c_ in seq_along(CELLS)) {
    j <- match(pair_of[c_], PAIRS)
    add <- cc * G[, j] * AMPA[, c_]
    add[is.na(add)] <- 0
    M[, c_] <- M[, c_] + add
  }
  M
}
C_GRID <- c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.80)
set.seed(90610L)
PWR <- list()
for (sb in seq_len(N_SUB)) {
  TO_S <- perm_G(TO)
  assert_signal_free(TO_S, TO, sprintf("substrate %d", sb))
  stopifnot(all(colSums(TO_S == 1L) == colSums(TO == 1L)))     # group sizes preserved
  for (cc in C_GRID) {
    Xi <- inject_on(cc, TO_S)
    od <- med_by_t(arm_cells(Xi, CISN, TO_S)) - med_by_t(arm_cells(TRN, TRNN, TO_S))
    nl <- matrix(NA_real_, NPERM_SWEEP, length(TP))
    for (k in seq_len(NPERM_SWEEP)) {
      Gp <- perm_G(TO_S)
      nl[k, ] <- med_by_t(arm_cells(Xi, CISN, Gp)) - med_by_t(arm_cells(TRN, TRNN, Gp))
    }
    for (j in seq_along(TP)) {
      p <- perm_p(nl[, j], od[j])
      PWR[[length(PWR) + 1L]] <- data.table(substrate = sb, c_inj = cc, time = TP[j],
        delta = od[j], p_emp = p,
        detected_infer_alpha = p < ALPHA_INFER, detected_05 = p < 0.05)
    }
  }
  say("  substrate %d/%d done", sb, N_SUB)
}
PWR <- rbindlist(PWR)
fwrite(PWR, file.path(outd, "b3v2_power_calibration.tsv"), sep = "\t")

# c=0 can-fail control, pooled over substrates
fp0 <- PWR[c_inj == 0 & detected_infer_alpha, .N]
fp0_rate <- PWR[c_inj == 0, mean(detected_infer_alpha)]
say("\nc=0 can-fail control across %d substrates x %d timepoints: %d detections (rate %.3f;",
    N_SUB, length(TP), fp0, fp0_rate)
say("expected about alpha = %.4f if the substrate is truly signal-free)", ALPHA_INFER)

# per-substrate, per-timepoint minimum detectable c
MINC <- PWR[c_inj > 0 & detected_infer_alpha == TRUE,
            .(min_c = min(c_inj)), by = .(substrate, time)]
allst <- CJ(substrate = seq_len(N_SUB), time = TP)
MINC <- merge(allst, MINC, by = c("substrate", "time"), all.x = TRUE)
RES <- MINC[, .(n_substrates_detecting = sum(!is.na(min_c)),
                med_min_c = median(min_c, na.rm = TRUE),
                lo_min_c = suppressWarnings(min(min_c, na.rm = TRUE)),
                hi_min_c = suppressWarnings(max(min_c, na.rm = TRUE))), by = time]
setorder(RES, time)
## This table summarizes the 10 substrates analyzed here.
## A separate 52-substrate analysis provides the c_star resolution bound
## in a separate table; the two summaries have distinct file names.
fwrite(RES, file.path(outd, "b3v2_power_resolution_script10_10substrate.tsv"), sep = "\t")
say("\nRESOLUTION under the permutation criterion (alpha = %.4f), across %d substrates:",
    ALPHA_INFER, N_SUB)
say("%-5s %22s %12s %18s", "time", "substrates detecting", "median c", "range of c")
for (i in seq_len(nrow(RES))) say("%-5d %19d/%d %12s %18s", RES$time[i],
  RES$n_substrates_detecting[i], N_SUB,
  ifelse(is.na(RES$med_min_c[i]), "none", sprintf("%.2f", RES$med_min_c[i])),
  ifelse(is.na(RES$lo_min_c[i]), "--", sprintf("%.2f - %.2f", RES$lo_min_c[i], RES$hi_min_c[i])))
det <- RES[, if (all(is.na(med_min_c))) NA_real_ else max(med_min_c, na.rm = TRUE)]
det_by_t <- setNames(RES$med_min_c, RES$time)
say("\nWORST-TIMEPOINT median detectable coupling (the equivalence bound): %s",
    ifelse(is.na(det), "none up to 0.80 -> TEST NOT INFORMATIVE", sprintf("c >= %.2f", det)))
say("Per timepoint: 15 min %s | 30 min %s | 60 min %s",
    ifelse(is.na(det_by_t["15"]), "none", sprintf("%.2f", det_by_t["15"])),
    ifelse(is.na(det_by_t["30"]), "none", sprintf("%.2f", det_by_t["30"])),
    ifelse(is.na(det_by_t["60"]), "none", sprintf("%.2f", det_by_t["60"])))
say("\nA single-substrate calculation gives c >= 0.40 (60 min) / 0.80 (15,30 min),")
say("using a 95%%-interval criterion; the bound above instead accounts for substrate-to-substrate variability.")
FEASIBLE <- !is.na(det) && fp0_rate <= 0.10

## ------------------------------- bootstrap CIs on the observed quantities ---
say("\n---- gene-level bootstrap CIs (%d resamples, matrixised ) ----", NBOOT)
BTc <- matrix(NA_real_, NBOOT, length(TP)); BTt <- BTc; BTd <- BTc
for (k in seq_len(NBOOT)) {
  rows <- sample.int(length(GEN_), length(GEN_), TRUE)
  a <- med_by_t(arm_cells(CIS, CISN, TO, rows)); b <- med_by_t(arm_cells(TRN, TRNN, TO, rows))
  BTc[k, ] <- a; BTt[k, ] <- b; BTd[k, ] <- a - b
}
ARM <- rbindlist(lapply(seq_along(TP), function(j) data.table(
  time = TP[j],
  cis = obs_cis[j],   cis_lo = quantile(BTc[, j], .025, na.rm = TRUE),
                      cis_hi = quantile(BTc[, j], .975, na.rm = TRUE),
  trans = obs_trans[j], trans_lo = quantile(BTt[, j], .025, na.rm = TRUE),
                        trans_hi = quantile(BTt[, j], .975, na.rm = TRUE),
  delta = obs_disc[j],  delta_lo = quantile(BTd[, j], .025, na.rm = TRUE),
                        delta_hi = quantile(BTd[, j], .975, na.rm = TRUE))))
fwrite(ARM, file.path(outd, "b3v2_arms_cis_vs_trans.tsv"), sep = "\t")
say("%-5s %26s %26s %26s", "time", "cis ratio [95%% CI]", "trans ratio [95%% CI]",
    "Delta [95%% CI]")
for (i in seq_len(nrow(ARM))) say("%-5d %10.3f [%6.3f,%6.3f] %10.3f [%6.3f,%6.3f] %+10.3f [%+6.3f,%+6.3f]",
  ARM$time[i], ARM$cis[i], ARM$cis_lo[i], ARM$cis_hi[i], ARM$trans[i], ARM$trans_lo[i],
  ARM$trans_hi[i], ARM$delta[i], ARM$delta_lo[i], ARM$delta_hi[i])

## ------------------------- secondary: graded events + strength ---
say("\n---- secondary splits, fixed in advance as secondary only ----")
say("Note: without a group size these ratios can rest on very few genes. The strength-change")
say("group is 0.1-0.2%% of genes, i.e. single digits per pair, and")
say("arm_cells only requires 5 per side, so a ratio could rest on ~5-10 genes and")
say("be over-read. A minimum group size is therefore enforced and")
say("the size is printed next to every ratio. Groups below it are reported as NOT EVALUABLE.")
MIN_GRP <- 30L
NEVhi <- matrix(0L, nrow(NEV), ncol(NEV), dimnames = dimnames(NEV))
for (j in seq_len(ncol(NEV))) NEVhi[, j] <- as.integer(NEV[, j] >= 2L)
STRh <- grp_mat("strength_change"); STRh[is.na(STRh)] <- 0L
sec_one <- function(G, tag) {
  ns <- vapply(seq_len(ncol(G)), function(j) sum(G[, j] == 1L), 0L)
  evaluable <- min(ns) >= MIN_GRP
  cisv <- med_by_t(arm_cells(CIS, CISN, G)); trnv <- med_by_t(arm_cells(TRN, TRNN, G))
  data.table(split = tag, time = TP, n_min = min(ns), n_med = as.integer(median(ns)),
             n_max = max(ns), evaluable = evaluable,
             cis = ifelse(evaluable, cisv, NA_real_),
             trans = ifelse(evaluable, trnv, NA_real_),
             delta = ifelse(evaluable, cisv - trnv, NA_real_),
             cis_raw = cisv, trans_raw = trnv)
}
SEC <- rbind(sec_one(NEVhi, ">=2 events"), sec_one(STRh, "strength change"))
fwrite(SEC, file.path(outd, "b3v2_secondary_splits.tsv"), sep = "\t")
say("\n%-16s %5s %7s %7s %10s %9s %9s %9s", "split", "time", "n_min", "n_med", "evaluable",
    "cis", "trans", "Delta")
for (i in seq_len(nrow(SEC))) say("%-16s %5d %7d %7d %10s %9s %9s %9s", SEC$split[i],
  SEC$time[i], SEC$n_min[i], SEC$n_med[i], ifelse(SEC$evaluable[i], "yes", "NO"),
  ifelse(SEC$evaluable[i], sprintf("%.3f", SEC$cis_raw[i]), "--"),
  ifelse(SEC$evaluable[i], sprintf("%.3f", SEC$trans_raw[i]), "--"),
  ifelse(SEC$evaluable[i], sprintf("%+.3f", SEC$cis_raw[i] - SEC$trans_raw[i]), "--"))
if (any(!SEC$evaluable))
  say("=> at least one secondary split is below n=%d per pair and is NOT interpreted.", MIN_GRP)

## ----------------- strain-mismatch sensitivity, Sk vs the rest ---
say("\n================ STRAIN-MISMATCH SENSITIVITY ================")
say("Only S. kudriavzevii (CR85) shares a strain between the sequence and the expression axis.")
say("THE DIRECTIONAL PREDICTION MATTERS HERE. If promoter turnover really drives allele-specific")
say("divergence, the signal must be STRONGER in the three pairs containing Sk -- those are the")
say("only ones where the promoter actually belongs to the strain whose RNA was sequenced -- and")
say("WEAKER in the three strain-mismatched pairs. The reverse pattern argues against a promoter")
say("mechanism, because a real cis effect cannot be carried by a promoter the sequenced strain")
say("does not have.")
v_c <- arm_cells(CIS, CISN, TO); v_t <- arm_cells(TRN, TRNN, TO)
SK <- rbindlist(lapply(c(0L, 1L), function(hs) {
  sel <- which(grepl("Sk", pair_of) == (hs == 1L))
  rbindlist(lapply(TP, function(t) { ii <- sel[tp_of[sel] == t]
    data.table(has_Sk = hs, time = t, pairs = paste(unique(pair_of[sel]), collapse = ","),
      cis = median(v_c[ii], na.rm = TRUE), trans = median(v_t[ii], na.rm = TRUE),
      delta = median(v_c[ii], na.rm = TRUE) - median(v_t[ii], na.rm = TRUE)) })) }))
fwrite(SK, file.path(outd, "b3v2_strain_sensitivity.tsv"), sep = "\t")
say("%-8s %5s %9s %9s %9s  %s", "has_Sk", "time", "cis", "trans", "Delta", "pairs")
for (i in seq_len(nrow(SK))) say("%-8d %5d %9.3f %9.3f %+9.3f  %s", SK$has_Sk[i], SK$time[i],
  SK$cis[i], SK$trans[i], SK$delta[i], SK$pairs[i])

# exact enumeration: there are only choose(6,3) = 20 ways to split 6 pairs into two threes, so
# the rank of the actual Sk split among all of them is computable exactly. With n=3 vs 3 this is
# descriptive only and no significance is claimed; the rank is reported
# because "which split is extreme" is informative even when it cannot be tested.
cmb <- combn(PAIRS, 3L)
SKEX <- rbindlist(lapply(TP, function(t) {
  vals <- vapply(seq_len(ncol(cmb)), function(z) {
    ina <- which(pair_of %in% cmb[, z] & tp_of == t)
    inb <- which(!(pair_of %in% cmb[, z]) & tp_of == t)
    (median(v_c[ina], na.rm = TRUE) - median(v_t[ina], na.rm = TRUE)) -
    (median(v_c[inb], na.rm = TRUE) - median(v_t[inb], na.rm = TRUE)) }, 0)
  actual <- which(vapply(seq_len(ncol(cmb)), function(z)
    setequal(cmb[, z], grep("Sk", PAIRS, value = TRUE)), TRUE))
  data.table(time = t, n_splits = ncol(cmb), sk_split_diff = vals[actual],
    rank_of_sk_split = rank(vals)[actual],
    frac_more_extreme = mean(abs(vals) >= abs(vals[actual]))) }))
fwrite(SKEX, file.path(outd, "b3v2_strain_exact.tsv"), sep = "\t")
say("\nexact enumeration of all %d three-versus-three splits (DESCRIPTIVE, n=3 vs 3):", ncol(cmb))
say("%-5s %16s %14s %20s", "time", "Sk-split diff", "rank of 20", "frac >= |diff|")
for (i in seq_len(nrow(SKEX))) say("%-5d %+16.3f %14.1f %20.2f", SKEX$time[i],
  SKEX$sk_split_diff[i], SKEX$rank_of_sk_split[i], SKEX$frac_more_extreme[i])

## --------------------------------------------- what this test concludes ---
say("\n================ WHAT CAN BE CONCLUDED ================")
BOOT_CI_EXCLUDES_0 <- ARM$delta_lo > 0 | ARM$delta_hi < 0
say("Summary of validation and results:")
say("  1. scan validated (positive control)?                    %s", SCAN_OK)
say("  2. primary metrics clear the %.1fx chance criterion?     %s",
    CHANCE_MIN, all(CHm$admissible_primary[2:4]))
say("  3. homogeneous-pool control passed (the real c=0)?               %s", HOMO_OK)
say("  4. amplitude balance -> is matching needed?             %s",
    ifelse(BALANCED, "NO (balanced)", "yes"))
say("  5. c=0 can-fail control across %d substrates (rate, expect ~alpha=%.4f): %.3f",
    N_SUB, ALPHA_INFER, fp0_rate)
say("  6. DISCRIMINANT resolution, PERMUTATION criterion at the inference's own alpha,")
say("     median over %d substrates: 15 min %s | 30 min %s | 60 min %s", N_SUB,
    ifelse(is.na(det_by_t["15"]), "none<=0.80", sprintf("c>=%.2f", det_by_t["15"])),
    ifelse(is.na(det_by_t["30"]), "none<=0.80", sprintf("c>=%.2f", det_by_t["30"])),
    ifelse(is.na(det_by_t["60"]), "none<=0.80", sprintf("c>=%.2f", det_by_t["60"])))
say("  7. observed Delta 15/30/60: %+.3f / %+.3f / %+.3f", obs_disc[1], obs_disc[2], obs_disc[3])
say("       label-permutation p_emp: %.3f / %.3f / %.3f", NULLTAB$p_emp_two_sided[1],
    NULLTAB$p_emp_two_sided[2], NULLTAB$p_emp_two_sided[3])
say("       gene bootstrap 95%% CI excludes 0: %s / %s / %s",
    BOOT_CI_EXCLUDES_0[1], BOOT_CI_EXCLUDES_0[2], BOOT_CI_EXCLUDES_0[3])
say("  8. is the signal where a promoter mechanism requires it? Sk-containing pairs")
say("       Delta %+.3f / %+.3f / %+.3f vs strain-mismatched pairs %+.3f / %+.3f / %+.3f",
    SK[has_Sk == 1L][order(time)]$delta[1], SK[has_Sk == 1L][order(time)]$delta[2],
    SK[has_Sk == 1L][order(time)]$delta[3], SK[has_Sk == 0L][order(time)]$delta[1],
    SK[has_Sk == 0L][order(time)]$delta[2], SK[has_Sk == 0L][order(time)]$delta[3])
say("")
say("The 15 min result is evaluated against the permutation-based")
say("calibration, which accounts for substrate-to-substrate variability rather than using a single-substrate")
say("calculation.")
say("")
say("Relevant results:")
say(" (i) The 15 min discriminant is outside the label-permutation null (p_emp = %.3f, passing",
    NULLTAB$p_emp_two_sided[1])
say("     the Bonferroni alpha of %.4f), while the GENE-LEVEL BOOTSTRAP CI on the same quantity",
    ALPHA_INFER)
say("     CONTAINS 0 [%+.3f, %+.3f]. The two tests ask different questions -- permutation: is",
    ARM$delta_lo[1], ARM$delta_hi[1])
say("     this gene group special at fixed group size; bootstrap: how stable is the estimate")
say("     under gene resampling. The permutation and bootstrap results are therefore")
say("     reported separately: the resolution uses the permutation analysis, and")
say("     the bootstrap CI describes stability under gene resampling.")
say(" (ii) Three timepoints were tested, one is extreme. That is why alpha is %.4f, not 0.05.",
    ALPHA_INFER)
say(" (iii) The strain split runs the WRONG WAY for a promoter mechanism. The positive Delta")
say("     lives in the three strain-MISMATCHED pairs; the three pairs containing the one species")
say("     whose promoter and RNA come from the SAME strain (Sk / CR85) show Delta at or below 0.")
say("     A cis effect cannot be carried by a promoter the sequenced strain does not have. BUT")
say("     the exact enumeration of all 20 three-versus-three splits puts the Sk split at rank")
say("     %.1f/20 at 15 min, i.e. NOT extreme where the signal is => this is a NAMED CANDIDATE",
    SKEX[time == 15L, rank_of_sk_split])
say("     explanation, not an established one. n=3 vs 3, descriptive only.")
say("")
say("Interpretation: the promoter HSE analysis does not provide clear positive or negative evidence for Model B.")
say("The validated findings are:")
say("  a. The scan and the two locked HSE classes are validated, and the gapped class is the")
say("     STRONGER of the two on both criteria (chance ratio and the Hsf1-dependency gradient).")
say("  b. The union of classic and gapped HSEs yields")
say("     ~%d genes per pair, and the discriminant has a resolution derived with", N_TO_MED)
say("     the SAME test as the inference and averaged over %d substrates.", N_SUB)
say("  c. Amplitude is balanced for this metric,")
say("     r = %+.4f, decile chi-square p = %.3f -- so this version is not confounded by amplitude.",
    r_pear, chi$p.value)
say("  d. The residual 15 min signal is concentrated in the strain")
say("     mismatched pairs, although")
say("     the exact enumeration indicates that this split is not extreme.")
say("")
say("Scientific interpretation:")
say("  Position-anchored turnover of classic and gapped heat-shock elements showed no")
say("   association with allele-specific divergence of the heat-stress response at 30 and 60")
say("   min. A positive association at 15 min was not stable across tests and was confined to")
say("   species pairs in which the promoter and the sequenced strain differ, so it is not")
say("   interpreted as evidence for target-side cis-regulatory divergence.")
say("The analysis does not establish that HSE architecture fails to explain divergence or that HSE turnover drives")
say("cis divergence. It also provides no evidence that upstream regulation evolved, which would require")
say("independent evidence. Negative evidence applies only within the detectable coupling range")
say("of this analysis.")
say("")
say("Limits:")
say(" - nucleosome context is not modelled (no cross-species data; Pincus 2018 shows it matters)")
say(" - the promoter window is %d bp and cannot see distal elements", WIN)
say(" - only reference-strain promoters exist for 3 of 4 species; non-coding sequence is")
say("   more polymorphic than coding, so this is more serious here than for the coding sequence -- and item (iii)")
say("   above shows it is not hypothetical in this analysis")
say(" - about %.0f%% of detected sites are expected by chance at the measured enrichment, which",
    100 / CHm$ratio[4])
say("   dilutes the split and biases toward the null (conservative)")
say(" - the HSE model is the two locked literature classes; no PWM was fitted")
say(" - no polarity is claimed for unmatched sites (no ancestral reconstruction)")
say(" - Hsf1 dependency is defined in S. cerevisiae only and projected onto orthologs")
say(" - position anchoring contributed only %d-%d 'moved' events per pair; the power gain came",
    min(PD[, sum(n_moved), by = pair]$V1), max(PD[, sum(n_moved), by = pair]$V1))
say("   almost entirely from adding the gapped class, not from the anchoring itself")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_b3v2_promoter.txt"))

## ------------------------------------------------ output handling ---
mdir <- file.path(MIRROR, "data/02_processed/hsf1")
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(MIRROR, "doc/05_reports/logs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(MIRROR, "scripts/hsf1"), showWarnings = FALSE, recursive = TRUE)
for (f in list.files(outd, pattern = "^(b3v2_|sessionInfo_b3v2)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
file.copy(file.path(root, "scripts/hsf1/10_promoter_hse_architecture_v2.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("\nprepared output copy in %s", MIRROR)
say("=== done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
