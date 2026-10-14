# 12_upstream_sensor_conservation.R -- a conservation check on two named upstream components
# across the four Fay reference genomes already in the repository:
#   SSA1 (YAL005C), cytosolic Hsp70 -- and specifically the residue corresponding to
#     S. cerevisiae Ssa1 T492, reported as a phosphorylation site in
#     PMID 41390490 (Nat Commun 2025, doi 10.1038/s41467-025-67204-7), where the
#     mechanosensor Mid2 senses heat-induced membrane stretch and triggers Ssa1 T492
#     phosphorylation, which rearranges the epichaperome and tunes HSF activity,
#     translational pausing, MAPK signalling and stress-granule resolution.
#   MID2 (YLR332W), the mechanosensor itself.
#
# Rationale. Divergence is not concentrated in the Hsf1-associated core (script 04), trans
# exceeds cis with the cis share falling with time (script 05), and the Hsf1 DBD is 97.2%
# identical (script 08), so the leading model is "conserved Hsf1 module with rewiring somewhere
# else". PMID 41390490 and Keyport Kik 2024 (PMID 38605014, which proposes that "upstream
# sensors may include thermally triggered biomolecular condensates which titrate
# Hsf1-repressive Hsp70") both point at the same layer, and two of its components can be
# checked in genomes already on disk.
#
# Decision rule, fixed before this analysis was run: the result is purely descriptive, and
# whichever way it comes out it is not evidence for upstream evolution:
#   - if T492 is conserved, that does not support Model A (a conserved residue is the null
#     expectation for a functional phosphosite and says nothing about regulatory divergence)
#   - if T492 is not conserved, that does not support Model A either, because there is no
#     phosphoproteomics, no kinase assay, no reporter and no cross-species perturbation, and
#     sequence alone does not support a claim about function
#   - Sequence conservation alone cannot establish upstream sensor
#     evolution.
# The output is a descriptive sequence observation and a candidate-list entry.
#
# Limitations, unchanged by adding these genes:
#   - strain mismatch: these are species reference strains DBVPG6765/CBS432/CR85/YJF1449,
#     while the expression data are T73/N17/CR85/BMV58. Only S. kudriavzevii matches. Every
#     statement is at the species-reference level.
#   - n = 4 species: no selection test, no dN/dS, no ancestral reconstruction.
#   - Only the abstract and data-availability statement were available for PMID 41390490, so its T492 claim is
#     evaluated at abstract-level detail.
#
# Method. The genetic code, FASTA reader, translate and revcomp are parsed out of
# scripts/hsf1/07_hsf1_locus_extract.R rather than copied, so they cannot drift from the
# functions that produced the HSF1 locus result, and script 07's own assertions are re-run here.
# Identity and the focal residue use no aligner: if the four proteins are not the same length
# the script reports that and does not compare positions, rather than guessing an offset,
# because position-by-position comparison of unequal-length proteins pairs non-homologous
# residues.
#
# Usage: Rscript scripts/hsf1/12_upstream_sensor_conservation.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
gdir <- file.path(root, "data/03_external/fay2023_genomes")
hsd  <- file.path(root, "data/02_processed/hsf1")
repd <- file.path(root, "doc/05_reports")
outd <- hsd
MIRROR <- "logs/output_copy"
SRC07 <- file.path(root, "scripts/hsf1/07_hsf1_locus_extract.R")
SRC08 <- file.path(root, "scripts/hsf1/08_hsf1_domain_divergence.R")   # reuse the script-08 aligner (MID2 only)

GEN <- data.table(
  species = c("Sc", "Sp", "Sk", "Su"),
  strain  = c("DBVPG6765", "CBS432", "CR85", "YJF1449"),
  matched_strain = c(FALSE, FALSE, TRUE, FALSE))
# systematic names; the reported residue is 1-based in the S. cerevisiae protein
TARGETS <- data.table(
  gene = c("SSA1", "MID2"),
  sys  = c("YAL005C", "YLR332W"),
  focal_residue = c(492L, NA_integer_),
  expected_aa = c("T", NA_character_),
  source = c("PMID 41390490 (Ssa1 T492 phosphosite)", "PMID 41390490 (mechanosensor)"))

say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, "logs/hsf1_upstream_sensor_conservation.log")
sink(logf, split = TRUE)
say("=== 12_upstream_sensor_conservation | %s ===", format(Sys.time()))
say("This analysis is purely descriptive. Neither conservation nor divergence of these")
say("residues establishes upstream regulatory evolution or supports Model A.")

## ------- transplant the code/translate/reader from script 07, then verify ---
say("\n---- transplanting genetic code + translate + readers out of %s ----", basename(SRC07))
ex <- parse(SRC07)
want_obj <- c("bases", "g", "codons", "aa", "CODE")
want_fun <- c("translate", "revcomp", "read_fasta")
got <- character(0)
for (e in ex) {
  if (is.call(e) && length(e) >= 3L && as.character(e[[1]]) %in% c("<-", "=") &&
      is.name(e[[2]]) && as.character(e[[2]]) %in% c(want_obj, want_fun)) {
    eval(e, envir = globalenv()); got <- c(got, as.character(e[[2]]))
  }
}
say("transplanted: %s", paste(got, collapse = ", "))
missing <- setdiff(c(want_obj, want_fun), got)
if (length(missing)) stop(sprintf("script 07 no longer defines: %s", paste(missing, collapse = ", ")))
# script 07's semantics, asserted here so the transplant is provably faithful
stopifnot(length(CODE) == 64L, CODE[["ATG"]] == "M", CODE[["TAA"]] == "*",
          CODE[["TGA"]] == "*", CODE[["TAG"]] == "*",
          translate("ATGGCTTAA") == "MA*",
          translate("atgGCTtaa") == "MA*",              # soft-masked lowercase must survive
          revcomp("ATGC") == "GCAT")
say("script 07 code/translate assertions re-run: PASS")

## ------------------------------------------------- extract + translate ---
# Exact attribute parsing. Note: matching CDS to their gene with
#   grepl(paste0("Parent=", gid), attr, fixed = TRUE)
# over-matches by prefix on numeric ids: Parent=SC_59 also matches Parent=SC_590 ... SC_599.
# For Sc SSA1 that assembles 3381 nt out of a gene whose span is only 1929 bp, giving a 1127 aa
# "protein" with an internal stop and no terminal stop. Attributes are therefore parsed into
# exact key=value pairs, Parent is matched exactly (allowing the comma-separated multi-parent
# form), and the assembled CDS is checked against the gene span, frame, terminal stop and
# internal stops before anything is reported. The length guard is a last check, not the fix.
attr_get <- function(a, key) {
  out <- rep(NA_character_, length(a))
  parts <- strsplit(a, ";", fixed = TRUE)
  for (i in seq_along(parts)) {
    kv <- parts[[i]]
    hit <- kv[startsWith(kv, paste0(key, "="))]
    if (length(hit)) out[i] <- sub(paste0("^", key, "="), "", hit[1])
  }
  out
}
parent_matches <- function(a, gid) {
  p <- attr_get(a, "Parent")
  vapply(p, function(x) !is.na(x) && gid %in% trimws(strsplit(x, ",", fixed = TRUE)[[1]]), TRUE)
}
say("\n---- extracting CDS for %s across 4 reference genomes ----",
    paste(TARGETS$gene, collapse = ", "))
PROT <- list(); QC <- list()
for (i in seq_len(nrow(GEN))) {
  sp <- GEN$species[i]; strain <- GEN$strain[i]
  gff <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(gdir, paste0(strain, ".final.gff")))),
               sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gff, c("chr","src","type","start","end","score","strand","frame","attr"))
  fa <- read_fasta(file.path(gdir, paste0(strain, ".fasta")))
  gff[, gname := attr_get(attr, "Gene_Name")]
  gff[, sysn := sub("_.*$", "", gname)]
  for (t in seq_len(nrow(TARGETS))) {
    sysn_t <- TARGETS$sys[t]; gname_t <- TARGETS$gene[t]
    gi <- gff[type == "gene" & sysn == sysn_t]
    if (nrow(gi) != 1L) {
      say("  %-3s %-5s: gene records found = %d -> SKIPPED (not a clean single copy)",
          sp, gname_t, nrow(gi))
      QC[[length(QC)+1L]] <- data.table(species = sp, strain = strain, gene = gname_t,
        sys = sysn_t, status = sprintf("gene_records=%d", nrow(gi)), n_cds = NA_integer_,
        cds_nt = NA_integer_, gene_span = NA_integer_, prot_aa = NA_integer_,
        start_M = NA, end_stop = NA, internal_stop = NA, in_frame = NA,
        within_span = NA, n_X = NA_integer_)
      next
    }
    gid <- attr_get(gi$attr[1], "ID")
    span <- gi$end[1] - gi$start[1] + 1L
    cds <- gff[type == "CDS" & chr == gi$chr[1] & strand == gi$strand[1] &
               parent_matches(attr, gid)][order(start)]
    if (nrow(cds) == 0L) {
      # fall back to containment, but only within the gene's own span
      cds <- gff[type == "CDS" & chr == gi$chr[1] & strand == gi$strand[1] &
                 start >= gi$start[1] & end <= gi$end[1]][order(start)]
      if (nrow(cds)) say("  %-3s %-5s: no Parent match, fell back to containment within the gene span",
                         sp, gname_t)
    }
    if (nrow(cds) == 0L) { say("  %-3s %-5s: no CDS -> SKIPPED", sp, gname_t); next }
    cs <- fa[[gi$chr[1]]]
    pieces <- vapply(seq_len(nrow(cds)), function(k) substr(cs, cds$start[k], cds$end[k]), "")
    cds_seq <- if (gi$strand[1] == "+") paste(pieces, collapse = "")
               else revcomp(paste(rev(pieces), collapse = ""))
    prot <- translate(cds_seq)
    prot_noStop <- sub("\\*$", "", prot)
    nX <- lengths(regmatches(prot_noStop, gregexpr("X", prot_noStop)))
    chk <- list(start_M = substr(prot, 1, 1) == "M",
                end_stop = grepl("\\*$", prot),
                internal_stop = grepl("\\*", prot_noStop),
                in_frame = nchar(cds_seq) %% 3L == 0L,
                within_span = nchar(cds_seq) <= span)
    QC[[length(QC)+1L]] <- data.table(species = sp, strain = strain, gene = gname_t,
      sys = sysn_t, status = "extracted", n_cds = nrow(cds), cds_nt = nchar(cds_seq),
      gene_span = span, prot_aa = nchar(prot_noStop), start_M = chk$start_M,
      end_stop = chk$end_stop, internal_stop = chk$internal_stop, in_frame = chk$in_frame,
      within_span = chk$within_span, n_X = nX)
    clean <- chk$start_M && chk$end_stop && !chk$internal_stop && chk$in_frame && chk$within_span
    say("  %-3s %-5s (%-9s) chr=%-10s %s pieces=%2d  %5d nt (span %5d)  %4d aa | M=%s stop=%s int*=%s frame=%s inspan=%s X=%d %s",
        sp, gname_t, strain, gi$chr[1], gi$strand[1], nrow(cds), nchar(cds_seq), span,
        nchar(prot_noStop), chk$start_M, chk$end_stop, chk$internal_stop, chk$in_frame,
        chk$within_span, nX, ifelse(clean, "", "<- REJECTED"))
    # only a clean single ORF is admitted; a failing extraction is dropped, not "used carefully"
    if (clean) PROT[[paste(gname_t, sp, sep = "|")]] <- prot_noStop
  }
}
QC <- rbindlist(QC, fill = TRUE)
fwrite(QC, file.path(outd, "upstream_sensor_locus_qc.tsv"), sep = "\t")

## ------------------------------- per gene: identity + the focal residue ---
RES <- list(); FOC <- list()
for (t in seq_len(nrow(TARGETS))) {
  gn <- TARGETS$gene[t]
  ks <- paste(gn, GEN$species, sep = "|")
  have <- ks[ks %in% names(PROT)]
  say("\n================ %s (%s) ================", gn, TARGETS$sys[t])
  if (length(have) < 4L) {
    say("only %d/4 species yielded a protein -> NOT COMPARED (no partial-set claims).", length(have))
    next
  }
  P <- unlist(PROT[ks]); names(P) <- GEN$species
  lens <- nchar(P)
  say("protein lengths: %s", paste(sprintf("%s %d", names(lens), lens), collapse = " | "))
  equal_len <- length(unique(lens)) == 1L
  say("all four the same length: %s", equal_len)
  sc <- strsplit(P[["Sc"]], "", fixed = TRUE)[[1]]
  if (equal_len) {
    idn <- vapply(c("Sp","Sk","Su"), function(s) {
      o <- strsplit(P[[s]], "", fixed = TRUE)[[1]]; mean(o == sc) }, 0)
    say("ungapped identity vs Sc (valid because all four lengths are equal):")
    say("  Sp %.1f%% | Sk %.1f%% | Su %.1f%%", 100*idn[1], 100*idn[2], 100*idn[3])
  } else {
    same <- names(lens)[lens == lens[["Sc"]]]
    idn <- setNames(rep(NA_real_, 3L), c("Sp","Sk","Su"))
    say("lengths differ => full-length ungapped identity is NOT computed for the differing")
    say("species, because position-based comparison of unequal-length proteins is the trap-8")
    say("failure mode. Reported only for the species matching Sc in length: %s",
        paste(setdiff(same, "Sc"), collapse = ", "))
    for (s in setdiff(same, "Sc")) {
      o <- strsplit(P[[s]], "", fixed = TRUE)[[1]]
      idn[s] <- mean(o == sc)
      say("  %s %.1f%% (equal length, ungapped)", s, 100*idn[s])
    }
    for (s in setdiff(c("Sp","Sk","Su"), same))
      say("  %s -- length %d vs Sc %d: requires an alignment, NOT REPORTED", s, lens[[s]], lens[["Sc"]])
  }
  # ---- focal residue by context anchoring, the same logic as for HSE sites in script 10 ----
  # No aligner: take the Sc window centred on the focal residue and require exactly one
  # occurrence in each other species. A unique hit makes the correspondence unambiguous; zero or
  # multiple hits are reported as unresolved rather than resolved by guesswork. This reads a
  # single coordinate across proteins of unequal length without a traceback-based aligner,
  # whose state encoding is error-prone (see script 10).
  fr <- TARGETS$focal_residue[t]
  focal_ok <- NA
  if (!is.na(fr)) {
    if (fr > lens[["Sc"]]) {
      say("focal residue %d is beyond the Sc protein length %d -> CANNOT CHECK", fr, lens[["Sc"]])
    } else {
      ANCH <- 7L                                 # +/-7 aa => a 15-mer anchor
      w1 <- max(1L, fr - ANCH); w2 <- min(lens[["Sc"]], fr + ANCH)
      anchor <- substr(P[["Sc"]], w1, w2)
      off_in_anchor <- fr - w1 + 1L
      say("focal residue %d (reported %s in Sc, %s)", fr, TARGETS$expected_aa[t], TARGETS$source[t])
      say("  Sc anchor (%d-mer, focal at anchor position %d): %s",
          nchar(anchor), off_in_anchor, anchor)
      aas <- setNames(rep(NA_character_, 4L), GEN$species)
      pos <- setNames(rep(NA_integer_, 4L), GEN$species)
      nhit <- setNames(rep(NA_integer_, 4L), GEN$species)
      aas[["Sc"]] <- substr(P[["Sc"]], fr, fr); pos[["Sc"]] <- fr; nhit[["Sc"]] <- 1L
      for (s in c("Sp","Sk","Su")) {
        hits <- gregexpr(anchor, P[[s]], fixed = TRUE)[[1]]
        h <- if (hits[1] == -1L) integer(0) else as.integer(hits)
        nhit[[s]] <- length(h)
        if (length(h) == 1L) {
          pos[[s]] <- h + off_in_anchor - 1L
          aas[[s]] <- substr(P[[s]], pos[[s]], pos[[s]])
        }
      }
      say("  %-3s %-6s %-9s %s", "sp", "anchor", "position", "residue")
      for (s in GEN$species) say("  %-3s %-6s %-9s %s", s,
        ifelse(nhit[[s]] == 1L, "unique", sprintf("%d hits", nhit[[s]])),
        ifelse(is.na(pos[[s]]), "--", as.character(pos[[s]])),
        ifelse(is.na(aas[[s]]), "unresolved", aas[[s]]))
      matches_expected <- aas[["Sc"]] == TARGETS$expected_aa[t]
      say("  Sc residue matches the published identity (%s): %s %s",
          TARGETS$expected_aa[t], matches_expected,
          ifelse(matches_expected, "", "<- MISMATCH: coordinate or annotation differs; do NOT report"))
      resolved <- !any(is.na(aas))
      say("  all four resolved by unique anchor: %s", resolved)
      if (resolved) {
        focal_ok <- length(unique(aas)) == 1L
        say("  identical in all four species: %s", focal_ok)
      } else {
        say("  => NOT RESOLVED for at least one species; no statement is made about this site.")
      }
      for (s in GEN$species) if (!is.na(pos[[s]])) {
        a <- max(1L, pos[[s]] - ANCH); b <- min(lens[[s]], pos[[s]] + ANCH)
        say("  %s  %s", s, substr(P[[s]], a, b))
      }
      FOC[[gn]] <- data.table(gene = gn, sc_position = fr,
        Sc = aas[["Sc"]], Sp = aas[["Sp"]], Sk = aas[["Sk"]], Su = aas[["Su"]],
        pos_Sp = pos[["Sp"]], pos_Sk = pos[["Sk"]], pos_Su = pos[["Su"]],
        anchor = anchor, anchor_unique_all = resolved,
        matches_published_Sc = matches_expected, identical_all4 = focal_ok)
    }
  }
  RES[[gn]] <- data.table(gene = gn, comparable_full_length = equal_len,
    len_Sc = lens[["Sc"]], len_Sp = lens[["Sp"]], len_Sk = lens[["Sk"]], len_Su = lens[["Su"]],
    ident_vs_Sc_Sp = idn[["Sp"]], ident_vs_Sc_Sk = idn[["Sk"]], ident_vs_Sc_Su = idn[["Su"]],
    focal_ok = focal_ok)
}
RES <- rbindlist(RES, fill = TRUE); fwrite(RES, file.path(outd, "upstream_sensor_identity.tsv"), sep = "\t")
if (length(FOC)) fwrite(rbindlist(FOC), file.path(outd, "upstream_sensor_focal_residue.tsv"), sep = "\t")

## ----------------------- MID2: descriptive comparison only ------
# MID2 is compared descriptively using script 08's
# Needleman-Wunsch affine-gap aligner. A
# single descriptive comparison is provided for a supplementary table and does not test
# upstream evolution because Mid2
# sequence identity is uninformative about regulatory divergence in either direction;
# (2) Mid2 is a membrane mechanosensor whose function depends on structure, O-glycosylation and
# a Ser/Thr-rich region, so percent identity is especially uninterpretable for it.
# The aligner is parsed out of script 08, so no alignment code is written here. What is written
# here is verification of the transplant, which is test code, not alignment code.
say("\n================ MID2: DESCRIPTIVE SEQUENCE COMPARISON ================")
ex8 <- parse(SRC08)
want8 <- c("AAO", "b62", "SM", "sc", "GO", "GE", "nw")
got8 <- character(0)
for (e in ex8) {
  if (is.call(e) && length(e) >= 3L && as.character(e[[1]]) %in% c("<-", "=") &&
      is.name(e[[2]]) && as.character(e[[2]]) %in% want8) {
    eval(e, envir = globalenv()); got8 <- c(got8, as.character(e[[2]]))
  }
}
say("transplanted from script 08 (zero new alignment code): %s", paste(got8, collapse = ", "))
if (length(setdiff(want8, got8))) {
  say("*** script 08 no longer defines %s -> MID2 not compared.",
      paste(setdiff(want8, got8), collapse = ", "))
} else {
  # verification of the transplant (test code, not alignment code)
  stopifnot(isSymmetric(SM), SM["W","W"] == 11, SM["C","C"] == 9, GO == 11, GE == 1)
  t_id <- nw("MKVLAAG", "MKVLAAG")
  stopifnot(length(t_id$a) == 7L, !any(t_id$a == "-"), !any(t_id$b == "-"))
  t_gap <- nw("MKVLAAG", "MKVAAG")                       # one deletion in the second sequence
  stopifnot(length(t_gap$a) == length(t_gap$b), sum(t_gap$b == "-") == 1L,
            paste(t_gap$a, collapse = "") == "MKVLAAG")
  t_sym <- nw("ACDEFGHIK", "ACDEGHIK")
  stopifnot(abs(t_sym$score - nw("ACDEGHIK", "ACDEFGHIK")$score) < 1e-9)
  say("aligner transplant verified: identity / single-deletion / score symmetry: PASS")

  gn <- "MID2"; ks <- paste(gn, GEN$species, sep = "|")
  if (!all(ks %in% names(PROT))) {
    say("MID2 missing for at least one species -> NOT COMPARED.")
  } else {
    P <- unlist(PROT[ks]); names(P) <- GEN$species
    say("lengths: %s", paste(sprintf("%s %d", names(P), nchar(P)), collapse = " | "))
    MRES <- rbindlist(lapply(c("Sp","Sk","Su"), function(s) {
      al <- nw(P[["Sc"]], P[[s]])
      a <- al$a; b <- al$b
      aligned <- a != "-" & b != "-"
      ident <- sum(a[aligned] == b[aligned]) / sum(aligned)
      # low-complexity guard, exactly as script 08 applies it: a global alignment inside a
      # low-complexity stretch pairs residues by position, not by homology
      w <- paste(a[aligned], collapse = "")
      tri <- if (nchar(w) >= 3L) substring(w, 1:(nchar(w)-2L), 3:nchar(w)) else character(0)
      tb <- sort(table(tri), decreasing = TRUE)
      lc <- if (length(tb)) as.numeric(sum(tb[seq_len(min(3L, length(tb)))]) / length(tri)) else 1
      data.table(gene = gn, vs = s, score = al$score, aligned_cols = sum(aligned),
                 gaps_in_Sc = sum(a == "-"), gaps_in_other = sum(b == "-"),
                 identity = ident, lc_top3 = lc, align_reliable = lc < 0.5)
    }))
    fwrite(MRES, file.path(outd, "upstream_sensor_mid2_alignment.tsv"), sep = "\t")
    say("%-4s %8s %13s %11s %13s %9s %10s", "vs", "score", "aligned cols", "gaps Sc",
        "gaps other", "identity", "lc_top3")
    for (i in seq_len(nrow(MRES))) say("%-4s %8.0f %13d %11d %13d %8.1f%% %10.3f%s",
      MRES$vs[i], MRES$score[i], MRES$aligned_cols[i], MRES$gaps_in_Sc[i],
      MRES$gaps_in_other[i], 100*MRES$identity[i], MRES$lc_top3[i],
      ifelse(MRES$align_reliable[i], "", "  <- low complexity, alignment NOT reliable"))
    say("")
    say("STATUS: SUPPLEMENTARY-TABLE DESCRIPTIVE RESULT. This does not establish a regulatory")
    say("mechanism in either direction. Mid2 identity is especially difficult to interpret because its")
    say("function depends on structure, O-glycosylation and a Ser/Thr-rich region rather than on")
    say("primary-sequence identity, and BLOSUM/identity is evidence about substitution")
    say("similarity, not about functional effect.")
  }
}
say("\n================ WHAT THIS DOES AND DOES NOT SUPPORT ================")
say("Supported, and only this: a descriptive statement of per-species identity for SSA1 and MID2")
say("at the species-reference level, and a descriptive statement of the residue at the position")
say("corresponding to the published Ssa1 T492 phosphosite.")
say("")
say("Not supported:")
say(" - sequence identity does not establish upstream sensor evolution or support Model A")
say(" - treating a conserved T492 as support for a conserved upstream layer: conservation is the")
say("   NULL expectation for a functional phosphosite and carries no information about")
say("   regulatory divergence")
say(" - treating a substituted T492 as support for upstream evolution: there is no")
say("   phosphoproteomics, no kinase assay, no reporter, no cross-species perturbation")
say(" - any dN/dS, selection test or ancestral reconstruction at n = 4")
say(" - the T492 evidence from PMID 41390490 is limited to its abstract and")
say("   data-availability statement")
say("")
say("Reference-strain caveat applies unchanged: DBVPG6765/CBS432/CR85/YJF1449 here vs")
say("T73/N17/CR85/BMV58 on the expression axis; only S. kudriavzevii matches.")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_upstream_sensor.txt"))
mdir <- file.path(MIRROR, "data/02_processed/hsf1")
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(outd, pattern = "^(upstream_sensor|sessionInfo_upstream)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
dir.create(file.path(MIRROR, "scripts/hsf1"), showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(root, "scripts/hsf1/12_upstream_sensor_conservation.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("
copied outputs to %s", MIRROR)
say("=== done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
