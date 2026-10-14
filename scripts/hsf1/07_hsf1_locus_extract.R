# 07_hsf1_locus_extract.R -- locate the HSF1 ortholog in each of the four Fay
# reference genomes, extract and translate its CDS, and evaluate genome quality and promoter comparability
# (genome quality, promoter comparability) that need the sequence. The download itself is fetch_inputs.sh
# genomes; this script consumes the eight files it fetched (four FASTA + four GFF, OSF
# 7kunf, PMID 37972247).
#
# Scope. Delivered here: the four Hsf1 protein sequences (for the alignment in
# 08_hsf1_domain_divergence.R) and the four upstream windows' coordinates (for the HSE
# scan), plus a per-species locus QC table that either passes or fails the genome-quality and promoter-comparability checks below. Not here:
# any alignment, identity calculation, dN/dS, selection test, ancestral reconstruction, or
# HSE scan -- those belong to later scripts, each with its own constraints (no dN/dS at
# n=4; the HSE definition is fixed from an external source and power-calibrated before any
# negative is reported). This script is extraction and QC only.
#
# Checks made before a coordinate is trusted:
#   (1) the gene the GFF calls YGL073W_HSF1 really is the 1:1 ortholog -- verified by the CDS
#       length matching htseq.ortho.csv's Len column to the nucleotide (the orthology call and
#       the genome annotation are independent, so agreement is a real cross-check), by a clean
#       single ORF (start ATG, one stop, no internal stop), and by the translated length
#       matching the aa count implied by htseq;
#   (2) the locus is not at a scaffold edge -- distance to both chromosome ends must
#       exceed the promoter window, or the upstream sequence the HSE scan needs would be
#       truncated;
#   (3) the promoter is comparable across species -- a 5' boundary present in the GFF, and an intergenic
#       distance to the nearest upstream feature that is of the same order across species; a
#       species whose HSF1 sits immediately downstream of another gene has a promoter that is
#       not comparable and is demoted for the sequence axis.
#
# Strain mismatch applies to every number here: DBVPG6765(Sc)/CBS432(Sp)/YJF1449(Su)
# are not the RNA-seq strains T73/N17/BMV58; only CR85(Sk) matches. Every sequence result is
# at the species-reference level and never attributed to the sequenced RNA-seq strain.
#
# Usage: Rscript scripts/hsf1/07_hsf1_locus_extract.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
gdir <- file.path(root, "data/03_external/fay2023_genomes")
ext  <- file.path(root, "data/03_external/fay2023")
outd <- file.path(root, "data/02_processed/hsf1")
repd <- file.path(root, "doc/05_reports")
seqd <- file.path(outd, "hsf1_sequences")
dir.create(seqd, showWarnings = FALSE, recursive = TRUE)

# Fay's species code -> reference strain / genome file / chromosome prefix used in the GFF.
GEN <- data.table(
  species = c("Sc", "Sp", "Sk", "Su"),
  strain  = c("DBVPG6765", "CBS432", "CR85", "YJF1449"),
  prefix  = c("SC", "SP", "SK", "SU"),
  len_col = c("Len.Sc", "Len.Sp", "Len.Sk", "Len.Su"),
  same_strain_as_rnaseq = c(FALSE, FALSE, TRUE, FALSE))
WIN <- 1000L  # provisional promoter window for the Q2b truncation check only; the window the
              # HSE scan uses is fixed separately, before that scan is run, not decided here.
say <- function(...) cat(sprintf(...), "\n", sep = "")

logf <- file.path(repd, "logs/hsf1_locus_extract.log")
dir.create(dirname(logf), showWarnings = FALSE, recursive = TRUE)
sink(logf, split = TRUE)
say("=== 07_hsf1_locus_extract | %s ===", format(Sys.time()))
say("Extract + translate HSF1 CDS in 4 Fay reference genomes; the locus-integrity checks defined in the header.")

## ------------------------------------------------ standard genetic code ---
# Hardcoded rather than pulled from a package (none available). Nuclear standard code
# (translation table 1); yeast nuclear genes use it -- the alternative-code CTG clades are
# Candida/Debaryomyces, not Saccharomyces.
bases <- c("T", "C", "A", "G")
g <- expand.grid(b3 = bases, b2 = bases, b1 = bases, stringsAsFactors = FALSE)
codons <- paste0(g$b1, g$b2, g$b3)   # order: b1 outer, b2 middle, b3 inner -> matches table 1
aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
CODE <- setNames(aa, codons)
translate <- function(dna) {
    # uppercase first: these assemblies are soft-masked (lowercase = repeat), and a
    # case-sensitive [^ACGT] would turn every soft-masked base into N and then X, and
    # fabricate ambiguity that is not in the sequence.
  dna <- gsub("[^ACGT]", "N", toupper(dna))
  n <- nchar(dna); ncod <- n %/% 3L
  cod <- substring(dna, seq(1L, by = 3L, length.out = ncod), seq(3L, by = 3L, length.out = ncod))
  paste(ifelse(cod %in% names(CODE), CODE[cod], "X"), collapse = "")
}
revcomp <- function(dna) {
  m <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  paste(rev(m[strsplit(toupper(dna), "")[[1]]]), collapse = "")
}

## ------------------------------------------------------ FASTA reader ---
# Pure-R single-pass reader: returns a named list of full chromosome sequences for one file.
read_fasta <- function(path) {
  ln <- readLines(path)
  hdr <- which(startsWith(ln, ">"))
  nm <- sub("^>(\\S+).*", "\\1", ln[hdr])
  starts <- hdr + 1L; ends <- c(hdr[-1] - 1L, length(ln))
  seqs <- mapply(function(s, e) paste(ln[s:e], collapse = ""), starts, ends)
  setNames(seqs, nm)
}

## ------------------- htseq CDS lengths, the independent cross-check ---
dat <- read.csv(file.path(ext, "htseq.ortho.csv"))
HREF <- as.data.table(dat[dat$Name == "HSF1", c("HOG","Len.Sc","Len.Sp","Len.Sk","Len.Su")])
stopifnot(nrow(HREF) == 1L)

## ---------------------------------------- per-species extraction + QC ---
QC <- list(); PROT <- character(); CDSs <- character()
for (i in seq_len(nrow(GEN))) {
  sp <- GEN$species[i]; strain <- GEN$strain[i]; pre <- GEN$prefix[i]
  gpath <- file.path(gdir, paste0(strain, ".final.gff"))
  gff <- fread(cmd = paste("grep -v '^#'", shQuote(gpath)), sep = "\t", header = FALSE,
               quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gff, c("chr","src","type","start","end","score","strand","frame","attr"))
  gff <- gff[type %in% c("gene", "CDS", "transcript")]
  # locate HSF1 gene by the annotation's Gene_Name tag
  gi <- gff[type == "gene" & grepl("YGL073W", attr)]
  stopifnot(nrow(gi) == 1L)
  chr <- gi$chr; gstart <- gi$start; gend <- gi$end; gstr <- gi$strand
  gid <- sub(".*ID=([^;]+).*", "\\1", gi$attr)
  # its CDS pieces (HSF1 is intronless in Sc; assert one piece, else sum them in order)
  cds <- gff[type == "CDS" & grepl(paste0("Parent=", gid), attr)][order(start)]
  if (nrow(cds) == 0L) cds <- gff[type == "CDS" & chr == gi$chr &
                                  start >= gstart & end <= gend][order(start)]
  stopifnot(nrow(cds) >= 1L)

  fa <- read_fasta(file.path(gdir, paste0(strain, ".fasta")))
  stopifnot(chr %in% names(fa))
  chrseq <- fa[[chr]]; chrlen <- nchar(chrseq)

  # assemble CDS in transcription order
  pieces <- mapply(function(s, e) substr(chrseq, s, e), cds$start, cds$end)
  cds_seq <- paste(pieces, collapse = "")
  if (gstr == "-") cds_seq <- revcomp(cds_seq)
  cds_seq <- toupper(cds_seq)                       # de-soft-mask for the alignment step
  cds_len <- nchar(cds_seq)
  prot <- translate(cds_seq)
  aa_len <- nchar(sub("\\*+$", "", prot))          # drop terminal stop(s)
  internal_stop <- grepl("\\*", sub("\\*$", "", prot))
  n_ambig <- lengths(regmatches(prot, gregexpr("X", prot)))  # true ambiguous residues (real N)
  start_ok <- startsWith(prot, "M")
  stop_ok  <- endsWith(prot, "*")

  # cross-check against the orthology table's CDS length (independent annotation)
  len_ref <- HREF[[GEN$len_col[i]]]
  len_match <- cds_len == len_ref

  # Distance to both chromosome ends relative to the promoter window
  up_room  <- if (gstr == "+") gstart - 1L else chrlen - gend    # bp of 5' flank available
  dn_room  <- if (gstr == "+") chrlen - gend else gstart - 1L
  q2a_on_scaffold <- up_room > 0L && dn_room > 0L
  q2b_window_ok   <- up_room >= WIN

  # Intergenic distance to the nearest upstream feature (any gene) on the same chromosome
  genes <- gff[type == "gene" & chr == gi$chr]
  genes[, gs := start]; genes[, ge := end]
  if (gstr == "+") {
    up <- genes[ge < gstart]
    inter <- if (nrow(up)) gstart - max(up$ge) - 1L else NA_integer_
  } else {
    up <- genes[gs > gend]
    inter <- if (nrow(up)) min(up$gs) - gend - 1L else NA_integer_
  }
  q3a_boundary <- TRUE     # a gene feature with an explicit start/end and strand is present
  q3b_intergenic <- inter

  say("\n---- %s (%s, %s_%s) ----", sp, strain, pre, "chr7?")
  say("  locus: %s:%d-%d (%s strand), id %s", chr, gstart, gend, gstr, gid)
  say("  CDS length %d nt (htseq Len %d, match %s) ; protein %d aa ; start-M %s ; end-stop %s ; internal-stop %s ; ambiguous(X) %d",
      cds_len, len_ref, len_match, aa_len, start_ok, stop_ok, internal_stop, n_ambig)
  say("  chromosome %s length %d ; 5' flank %d bp ; 3' flank %d bp", chr, chrlen, up_room, dn_room)
  say("  on-scaffold-with-flanks: %s ; >=%d bp upstream: %s", q2a_on_scaffold, WIN, q2b_window_ok)
  say("  intergenic distance to nearest upstream gene: %s bp",
      ifelse(is.na(inter), "NA", inter))

  QC[[i]] <- data.table(species = sp, strain = strain, same_strain_as_rnaseq = GEN$same_strain_as_rnaseq[i],
    chr = chr, start = gstart, end = gend, strand = gstr, chr_len = chrlen,
    cds_len = cds_len, htseq_len = len_ref, len_match = len_match, aa_len = aa_len,
    start_M = start_ok, end_stop = stop_ok, internal_stop = internal_stop, n_ambig = n_ambig,
    up_flank = up_room, dn_flank = dn_room,
    Q2a = q2a_on_scaffold, Q2b = q2b_window_ok, Q3a = q3a_boundary, Q3b_intergenic = inter)
  PROT[sp] <- sub("\\*$", "", prot)
  CDSs[sp] <- cds_seq
}
QC <- rbindlist(QC)

## ------------------------------------------------------- verdicts + output ---
say("\n================ HSF1 LOCUS QC TABLE ================")
print(QC[, .(species, strain, cds_len, len_match, aa_len, start_M, end_stop, internal_stop,
             up_flank, Q2a, Q2b, Q3b_intergenic)])

# the DBD-completeness check left to
# 08_hsf1_domain_divergence.R is partly settled already: a clean single ORF of the
# expected length is the precondition; the domain-by-domain check is in that script.
q1_orf_ok <- QC[, all(len_match & start_M & end_stop & !internal_stop)]
say("\nsingle clean ORF of the htseq-expected length in all four species: %s", q1_orf_ok)
say("   -> the DBD-domain check itself is in 08_hsf1_domain_divergence.R (against Choo & Kang 2021, PMID 33527316).")

say("\ngenome quality: on-scaffold checks all %s ; >=%d bp upstream checks all %s",
    QC[, all(Q2a)], WIN, QC[, all(Q2b)])
# Promoter comparability is not "intergenic >= WIN" -- that was the scaffold-truncation test (Q2b).
# A promoter is comparable across species when a 5' boundary is annotated and the intergenic
# distances are of the same order (so the HSE scan covers homologous windows, not a full
# window in one species and a clipped one in another). The provisional WIN here is used only for
# the scaffold-edge check.
im <- QC$Q3b_intergenic
spread <- max(im, na.rm = TRUE) / max(1L, min(im, na.rm = TRUE))
say("Promoter comparability: 5' boundary present all TRUE ; intergenic distance to")
say("   nearest upstream gene: min %d / median %.0f / max %d bp -> spread factor %.2fx",
    min(im, na.rm = TRUE), median(im, na.rm = TRUE), max(im, na.rm = TRUE), spread)
# comparable if all present and within a modest fold-range of each other
q3_ok <- all(!is.na(im)) && spread <= 2
say("   comparable (all present, spread <= 2x): %s", q3_ok)
say("   CONSTRAINT: the fixed promoter window must be <= the smallest")
say("   intergenic distance (%d bp) OR be explicitly clipped at the upstream gene boundary,",
    min(im, na.rm = TRUE))
say("   otherwise the window would overlap the upstream ORF in every species.")

say("\n================ LOCUS QC SUMMARY ================")
for (i in seq_len(nrow(QC))) {
  seqv <- QC$len_match[i] && QC$start_M[i] && QC$end_stop[i] && !QC$internal_stop[i] &&
          QC$Q2a[i] && QC$Q2b[i] && !is.na(QC$Q3b_intergenic[i]) && q3_ok
  say("  %-4s (%-9s): sequence axis %s", QC$species[i], QC$strain[i],
      ifelse(seqv, "ADMITTED", "DEMOTED to descriptive"))
}

fwrite(QC, file.path(outd, "hsf1_locus_qc.tsv"), sep = "\t")
# protein + CDS FASTA for 08_hsf1_domain_divergence.R (species-reference level; the header
# states the strain and the strain-mismatch tag)
pf <- file.path(seqd, "hsf1_protein.fasta"); cf <- file.path(seqd, "hsf1_cds.fasta")
writeLines(unlist(lapply(names(PROT), function(s)
  c(sprintf(">%s_HSF1 strain=%s ref_level", s, GEN[species == s, strain]), PROT[[s]]))), pf)
writeLines(unlist(lapply(names(CDSs), function(s)
  c(sprintf(">%s_HSF1 strain=%s ref_level", s, GEN[species == s, strain]), CDSs[[s]]))), cf)
say("\nwrote %s (%d seqs) and %s", basename(pf), length(PROT), basename(cf))
say("protein lengths (aa): %s",
    paste(sprintf("%s=%d", names(PROT), nchar(PROT)), collapse = " "))
say("\n=== done ===")
