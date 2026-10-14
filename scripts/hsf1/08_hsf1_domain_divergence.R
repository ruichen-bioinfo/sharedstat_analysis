# 08_hsf1_domain_divergence.R -- per-domain sequence divergence of the four
# Saccharomyces Hsf1 orthologs extracted by 07_hsf1_locus_extract.R. Addresses the trans
# side of the sequence question (does Hsf1 itself carry sequence change that could underlie
# the response divergence?) at the only resolution n=4 permits: descriptive per-domain
# identity plus a flagged list of divergent residues. It is a candidate list, never a
# conclusion.
#
# Limits:
#   - No dN/dS, no selection test, no ancestral reconstruction. n=4 orthologs has no power for
#     any of them, and forcing them would produce a control that cannot fail.
#   - BLOSUM/identity is evidence about substitution similarity, not about functional effect.
#     A conserved residue is not proof of function and a divergent one is not proof of
#     adaptation. The output lists candidate positions, nothing
#     more.
#   - Species-reference level only: three of the four genomes are not the RNA-seq
#     strain, so no residue difference is attributed to the sequenced strain.
#
# Domain definitions come from an external source. Coordinates
# are taken from UniProt P10961 (HSF_YEAST), the reviewed S. cerevisiae Hsf1 entry, which is
# 833 aa and matches the Sc translation from 07_hsf1_locus_extract.R exactly -- so its
# coordinates map onto our Sc sequence with no offset. Regions as annotated in UniProt
# release 2026_03 (fetched 2026-09-24):
#     DNA binding                170-259
#     Flexible linker            260-280
#     Disordered                 277-309  (compositional bias "Low complexity" 277-296)
#     Involved in trimerization  350-403
# and the remainder as disordered / compositionally biased regions.
# Note: UniProt labels 277-309 "Disordered" and places "Involved in trimerization" at
# 350-403, so identity in 277-309 says nothing about the trimerization region. The
# alignment does not depend on these labels; only the slicing of its columns into regions
# does. The activation-region scheme of Choo & Kang 2021 (AR1/AR2/CE2/CTM, PMID 33527316)
# is not adopted here: UniProt does not carry those names, and their Sc boundaries would
# have to be guessed. The biologically decisive contrast -- structured DBD/trimerization
# vs intrinsically disordered activation regions -- is fully answerable from the UniProt
# annotation, and a sliding-window identity is reported beside the domain table so the
# conclusion does not rest on any single boundary choice.
#
# Alignment. The four proteins were aligned with a
# Needleman-Wunsch global aligner with affine gaps (Gotoh) and BLOSUM62 is
# implemented here. This is safe because the four are Saccharomyces sensu stricto
# orthologs, near-colinear (length spread 1.6%, 07_hsf1_locus_extract.R); this
# global aligner would be unreliable for distant sequences but is appropriate for these.
# All four are anchored on the Sc coordinate frame (pairwise to Sc) so UniProt's Sc
# coordinates address alignment columns directly; insertions relative to Sc are counted,
# not silently dropped.
#
# Usage: Rscript scripts/hsf1/08_hsf1_domain_divergence.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
seqd <- file.path(root, "data/02_processed/hsf1/hsf1_sequences")
outd <- file.path(root, "data/02_processed/hsf1")
repd <- file.path(root, "doc/05_reports")
say <- function(...) cat(sprintf(...), "\n", sep = "")

logf <- file.path(repd, "logs/hsf1_domain_divergence.log")
dir.create(dirname(logf), showWarnings = FALSE, recursive = TRUE)
sink(logf, split = TRUE)
say("=== 08_hsf1_domain_divergence | %s ===", format(Sys.time()))
say("Per-domain identity of 4 Saccharomyces Hsf1 orthologs. DESCRIPTIVE, candidate-level.")
say("Domain coordinates externally locked to UniProt P10961 (833 aa, matches the S. cerevisiae locus exactly).")

## -------------------------- read the proteins from 07_hsf1_locus_extract.R ---
read_fasta <- function(path) {
  ln <- readLines(path); h <- which(startsWith(ln, ">"))
  nm <- sub("^>(\\S+).*", "\\1", ln[h]); s <- h + 1L; e <- c(h[-1] - 1L, length(ln))
  setNames(vapply(seq_along(h), function(i) paste(ln[s[i]:e[i]], collapse = ""), ""), nm)
}
P <- read_fasta(file.path(seqd, "hsf1_protein.fasta"))
names(P) <- sub("_HSF1$", "", names(P))
SPECIES <- c("Sc", "Sp", "Sk", "Su")
stopifnot(all(SPECIES %in% names(P)))
P <- P[SPECIES]
say("\nprotein lengths (aa): %s", paste(sprintf("%s=%d", names(P), nchar(P)), collapse = " "))
stopifnot(nchar(P["Sc"]) == 833L)   # must equal the UniProt reference length or coords shift

## ------------------------------------------------------------- BLOSUM62 ---
AAO <- strsplit("ARNDCQEGHILKMFPSTWYVBZX*", "")[[1]]
b62 <- "
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1"
SM <- matrix(scan(text = b62, quiet = TRUE), length(AAO), length(AAO), byrow = TRUE,
             dimnames = list(AAO, AAO))
stopifnot(isSymmetric(SM), SM["W","W"] == 11, SM["C","C"] == 9)
sc <- function(a, b) { a <- ifelse(a %in% AAO, a, "X"); b <- ifelse(b %in% AAO, b, "X"); SM[cbind(a, b)] }

## -------------------------------- Needleman-Wunsch, affine gaps (Gotoh) ---
GO <- 11; GE <- 1   # BLOSUM62 standard affine penalties
nw <- function(s1, s2) {
  a <- strsplit(s1, "")[[1]]; b <- strsplit(s2, "")[[1]]; n <- length(a); m <- length(b)
  NEG <- -1e9
  M <- matrix(NEG, n + 1L, m + 1L); Ix <- matrix(NEG, n + 1L, m + 1L); Iy <- matrix(NEG, n + 1L, m + 1L)
  Tm <- matrix(0L, n + 1L, m + 1L); Tx <- matrix(0L, n + 1L, m + 1L); Ty <- matrix(0L, n + 1L, m + 1L)
  M[1, 1] <- 0
  for (i in 2:(n + 1L)) Ix[i, 1] <- -GO - (i - 2L) * GE
  for (j in 2:(m + 1L)) Iy[1, j] <- -GO - (j - 2L) * GE
  smat <- outer(a, b, sc)   # n x m substitution scores
  for (i in 2:(n + 1L)) for (j in 2:(m + 1L)) {
    s <- smat[i - 1L, j - 1L]
    cand <- c(M[i-1,j-1], Ix[i-1,j-1], Iy[i-1,j-1]) + s
    k <- which.max(cand); M[i,j] <- cand[k]; Tm[i,j] <- k
    # pointers store the actual previous-state code (1=M, 2=Ix, 3=Iy), not the which.max
    # index: Iy is state 3, so a "extend Iy" must be recorded as 3, not 2.
    cx <- c(M[i-1,j] - GO, Ix[i-1,j] - GE); kx <- which.max(cx); Ix[i,j] <- cx[kx]; Tx[i,j] <- if (kx == 1L) 1L else 2L
    cy <- c(M[i,j-1] - GO, Iy[i,j-1] - GE); ky <- which.max(cy); Iy[i,j] <- cy[ky]; Ty[i,j] <- if (ky == 1L) 1L else 3L
  }
  # traceback from the best of the three end states
  ends <- c(M[n+1,m+1], Ix[n+1,m+1], Iy[n+1,m+1]); st <- which.max(ends)
  i <- n + 1L; j <- m + 1L; o1 <- character(0); o2 <- character(0)
  while (i > 1L || j > 1L) {
    # boundaries: only one kind of move is possible, and the pointer arrays are not defined
    # there, so consume the remainder directly rather than reading a (zero) pointer.
    if (j == 1L) { o1 <- c(a[i-1], o1); o2 <- c("-", o2); i <- i-1L; next }
    if (i == 1L) { o1 <- c("-", o1); o2 <- c(b[j-1], o2); j <- j-1L; next }
    if (st == 1L) { o1 <- c(a[i-1], o1); o2 <- c(b[j-1], o2); ns <- Tm[i,j]; i <- i-1L; j <- j-1L; st <- ns }
    else if (st == 2L) { o1 <- c(a[i-1], o1); o2 <- c("-", o2); ns <- Tx[i,j]; i <- i-1L; st <- ns }
    else { o1 <- c("-", o1); o2 <- c(b[j-1], o2); ns <- Ty[i,j]; j <- j-1L; st <- ns }
  }
  list(a = o1, b = o2, score = max(ends))
}

## ------------------------- Sc-anchored columns: partner residue per Sc pos ---
# For each non-Sc species, align it to Sc and record, for every Sc residue position, the
# aligned partner residue ("-" where that species has a deletion). Insertions relative to Sc
# are counted per species but not placed in the Sc frame (they cannot be, by definition).
L <- nchar(P["Sc"]); sc_res <- strsplit(P["Sc"], "")[[1]]
COL <- matrix("-", L, length(SPECIES), dimnames = list(NULL, SPECIES))
COL[, "Sc"] <- sc_res
ins_count <- setNames(integer(length(SPECIES)), SPECIES)
## The whole-protein identities below are Sp 91.0, Sk 81.4, and Su 78.5.
## They are written to a table as well as to the log, so each figure can be checked against
## a column; nothing about how they are computed differs from the logged values.
WP <- list()
for (sp in setdiff(SPECIES, "Sc")) {
  al <- nw(P["Sc"], P[sp])
  pos <- 0L
  for (k in seq_along(al$a)) {
    if (al$a[k] != "-") { pos <- pos + 1L; COL[pos, sp] <- al$b[k] }
    else ins_count[sp] <- ins_count[sp] + 1L   # residue in sp with no Sc position
  }
  gid <- mean(al$a == al$b) * 100
  say("aligned %s vs Sc: score %.0f, whole-protein identity %.1f%%, %d insertions rel. to Sc",
      sp, al$score, gid, ins_count[sp])
  WP[[length(WP) + 1L]] <- data.table(species = sp, reference = "Sc", align_score = al$score,
                                      whole_protein_identity_pct = gid,
                                      aligned_columns = length(al$a),
                                      insertions_vs_Sc = ins_count[sp])
}
WP <- rbindlist(WP)
fwrite(WP, file.path(outd, "hsf1_whole_protein_identity.tsv"), sep = "\t")
say("wrote hsf1_whole_protein_identity.tsv (%d rows) so these figures are column-addressable", nrow(WP))

## ------------------------------------------ domain table (UniProt-locked) ---
UNIPROT_RELEASE <- "2026_03"
## Non-overlapping partition of the Sc coordinates. Where UniProt regions overlap (linker 260-280 and disordered
## 277-309), the named functional region keeps its full span and the overlap goes to it; 310-349 carries no UniProt
## region and is its own row, so no identity is attributed to a region that is not annotated.
DOM <- data.table(
  domain = c("N-terminal (disordered)", "DBD", "flexible linker", "disordered (low complexity)",
             "unannotated", "trimerization", "C-terminal (disordered/AR)"),
  start  = c(1L,   170L, 260L, 281L, 310L, 350L, 404L),
  end    = c(169L, 259L, 280L, 309L, 349L, 403L, 833L),
  kind   = c("disordered", "structured", "linker", "disordered", "unannotated", "structured",
             "disordered/regulatory"),
  source = paste("UniProt P10961 release", UNIPROT_RELEASE))
stopifnot(DOM[domain == "DBD", start == 170L & end == 259L],
          DOM[domain == "flexible linker", start == 260L & end == 280L],
          DOM[domain == "trimerization", start == 350L & end == 403L],
          all(DOM$start[-1] == DOM$end[-nrow(DOM)] + 1L), DOM$start[1] == 1L, DOM$end[nrow(DOM)] == 833L)
say("domain coordinates: UniProt P10961 release %s", UNIPROT_RELEASE)

pair_ids <- combn(SPECIES, 2, simplify = FALSE)
dom_identity <- function(rows) {
  # mean over the 6 species pairs of identity across the given Sc positions, counting only
  # columns where both species have a residue (gaps excluded from the denominator)
  ids <- sapply(pair_ids, function(pr) {
    x <- COL[rows, pr[1]]; y <- COL[rows, pr[2]]; ok <- x != "-" & y != "-"
    if (!any(ok)) NA_real_ else mean(x[ok] == y[ok]) * 100
  })
  c(mean = mean(ids, na.rm = TRUE), min = min(ids, na.rm = TRUE), max = max(ids, na.rm = TRUE))
}
# low-complexity fraction of the Sc segment: sum of the three commonest residues. A global
# aligner cannot be trusted inside a low-complexity repeat (it will pair up identical residues
# by position rather than by homology), so any identity reported for such a segment is an
# artefact, not a measurement. Flag it rather than interpret it.
lc_top3 <- function(rows) {
  r <- COL[rows, "Sc"]; r <- r[r != "-"]
  if (!length(r)) return(NA_real_); tb <- sort(table(r), decreasing = TRUE)
  sum(head(tb, 3)) / length(r)
}
say("\n================ PER-DOMAIN IDENTITY (mean over 6 species pairs) ================")
say("%-38s %5s %5s %8s %8s %8s %7s", "domain", "start", "end", "id_mean", "id_min", "id_max", "lowcplx")
RES <- list()
for (i in seq_len(nrow(DOM))) {
  rows <- DOM$start[i]:DOM$end[i]
  st <- dom_identity(rows); lc <- lc_top3(rows)
  RES[[i]] <- data.table(domain = DOM$domain[i], kind = DOM$kind[i],
    start = DOM$start[i], end = DOM$end[i], n_aa = length(rows),
    id_mean = st["mean"], id_min = st["min"], id_max = st["max"], lc_top3 = lc,
    align_reliable = lc < 0.5)
  say("%-38s %5d %5d %7.1f%% %7.1f%% %7.1f%% %6.2f%s", DOM$domain[i], DOM$start[i], DOM$end[i],
      st["mean"], st["min"], st["max"], lc, ifelse(lc >= 0.5, " *LOWCPLX*", ""))
}
RES <- rbindlist(RES)
whole <- dom_identity(1:L)
say("%-38s %5d %5d %7.1f%% %7.1f%% %7.1f%%", "WHOLE PROTEIN", 1L, L,
    whole["mean"], whole["min"], whole["max"])
fwrite(RES, file.path(outd, "hsf1_domain_identity.tsv"), sep = "\t")

say("\nNOTE: 277-309 is 'Disordered' in UniProt %s, not the trimerization region, which is 350-403.",
    UNIPROT_RELEASE)
say("277-309 is low-complexity (top-3 residues = %.0f%% of it), so its identity is flagged as unreliable.",
    100 * lc_top3(277:309))
say("trimerization region 350-403: %.1f%% mean identity, top-3 residue fraction %.2f",
    RES[domain == "trimerization", id_mean], RES[domain == "trimerization", lc_top3])
say("\nPRIOR (Choo & Kang 2021, PMID 33527316): even Ogataea Hsf1 at 26%% whole-protein")
say("identity to Sc retains all functional domains. The one structured domain this analysis")
say("can align reliably -- the DBD -- is near-invariant across the four species, so the trans")
say("(Hsf1-protein) candidate for the response divergence does NOT lie in the DNA-contacting")
say("core. ")

## ------------------------------------------ divergent-residue candidate list ---
# A Sc position is flagged when at least one species differs from the Sc residue. This is the
# candidate list of positions, not a mechanism.
say("\n================ DIVERGENT-RESIDUE CANDIDATES ================")
divg <- which(sapply(1:L, function(p) {
  r <- COL[p, ]; r <- r[r != "-"]; length(unique(r)) > 1L
}))
dom_of <- function(p) { d <- DOM[start <= p & end >= p, domain]; if (length(d)) d[1] else NA }
DV <- data.table(sc_pos = divg, sc_aa = COL[divg, "Sc"],
                 Sp = COL[divg, "Sp"], Sk = COL[divg, "Sk"], Su = COL[divg, "Su"],
                 domain = sapply(divg, dom_of))
fwrite(DV, file.path(outd, "hsf1_divergent_residues.tsv"), sep = "\t")
say("total variable Sc positions: %d of %d (%.1f%%)", length(divg), L, 100 * length(divg) / L)
say("by domain:")
tab <- DV[, .(n_variable = .N), by = domain]
tab <- merge(DOM[, .(domain, n_aa = end - start + 1L, kind)], tab, by = "domain", all.x = TRUE)
tab[is.na(n_variable), n_variable := 0L]
tab[, pct_variable := 100 * n_variable / n_aa]
tab <- tab[match(DOM$domain, tab$domain)]
for (i in seq_len(nrow(tab))) say("  %-38s %3d/%3d variable (%.1f%%) [%s]",
  tab$domain[i], tab$n_variable[i], tab$n_aa[i], tab$pct_variable[i], tab$kind[i])
say("\nDBD variable residues (the candidates that would matter most if any) -- the DBD is the")
say("one structured, reliably-aligned domain, so these are the only structured-core candidates:")
key <- DV[domain == "DBD"]
if (nrow(key) == 0L) {
  say("  NONE. The DNA-binding domain is identical across all four species.")
} else {
  for (i in seq_len(nrow(key))) {
    cons <- sc(key$sc_aa[i], c(key$Sp[i], key$Sk[i], key$Su[i]))
    say("  Sc %d (%s): Sp=%s Sk=%s Su=%s  [BLOSUM62 vs Sc: %s] -> candidate, NOT a mechanism",
        key$sc_pos[i], key$sc_aa[i], key$Sp[i], key$Sk[i], key$Su[i],
        paste(cons, collapse = "/"))
  }
  say("  (Positive BLOSUM = conservative substitution. these are candidate residues")
  say("   only, not evidence of functional change; no functional assay exists.)")
}
say("The low-complexity 277-309 segment is not a structured core and its residues are not candidates here.")

## ---------------------------------------------- sliding-window identity ---
# So the domain conclusion does not rest on the boundary choice: 30-aa windows, step 10,
# mean pairwise identity, printed coarsely.
say("\n---- sliding-window identity (30 aa, step 10; guards against boundary artefacts) ----")
W <- 30L; STEP <- 10L
starts <- seq(1L, L - W + 1L, by = STEP)
SW <- rbindlist(lapply(starts, function(s) {
  st <- dom_identity(s:(s + W - 1L))
  data.table(win_start = s, win_end = s + W - 1L, id_mean = st["mean"])
}))
fwrite(SW, file.path(outd, "hsf1_sliding_identity.tsv"), sep = "\t")
say("windows overlapping DBD (170-259) and trimerization (350-403):")
for (i in seq_len(nrow(SW))) if (SW$win_end[i] >= 170L && SW$win_start[i] <= 403L)
  say("  %3d-%3d : %5.1f%%%s", SW$win_start[i], SW$win_end[i], SW$id_mean[i],
      ifelse(SW$win_start[i] <= 259 & SW$win_end[i] >= 170, " <DBD>",
             ifelse(SW$win_start[i] <= 403 & SW$win_end[i] >= 350, " <trimer>", "")))
say("min-identity window over whole protein: %.1f%% at %d-%d (the most divergent stretch)",
    min(SW$id_mean), SW$win_start[which.min(SW$id_mean)], SW$win_end[which.min(SW$id_mean)])
say("max-identity window: %.1f%% at %d-%d",
    max(SW$id_mean), SW$win_start[which.max(SW$id_mean)], SW$win_end[which.max(SW$id_mean)])

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_b2_domain.txt"))
say("\n=== done ===")
