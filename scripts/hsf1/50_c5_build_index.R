#!/usr/bin/env Rscript
## 50_c5_build_index.R -- extract CDS from the two parental references and build one salmon index
## containing both, so that a hybrid read is assigned to a parental allele rather than to a species genome.
##
## Strains, from the primary records. All eighteen BioSample records for PRJNA576452 give the same three strain
## strings: the hybrid is S288C x CBS7001, the S. cerevisiae parent is S288C, the S. uvarum parent is CBS7001.
## Note: the hybrid is not YPS128 x UWOPS99-807. It matters because CBS7001 is the S. uvarum type strain and
## UWOPS99-807 is a more distant isolate, which would change how far the reads sit from the reference.
##
## Why the P3 references and not the matching strains. The reads come from S288C and CBS7001; the references are
## Fay et al.'s DBVPG6765 and YJF1449, which is what the companion cross-species analysis uses. Keeping them
## means every gene here carries the same HOG identifier and the same Hsf1-dependence stratum as in that
## analysis. Building a new orthology map against S288C and CBS7001 assemblies would be a larger unverified step
## than the mapping asymmetry it removes.
##
## Why the strain mismatch largely cancels. The estimand is a response, the 12 C minus 30 C difference within one
## allele. A mapping penalty from strain distance is a constant for a fixed combination of allele, reference and
## read type, so it cancels to first order when the two temperatures are subtracted. It distorts the allele ratio
## at each temperature, not the change in that ratio. What survives is second order: expression-dependent
## mappability, small against three-replicate noise. This is the same reason the companion works in response
## space throughout and the reason its leakage term is a leakage and not a distance.
##
## Decision rule for this analysis (scripts 50 to 52), written before any quantification:
##   What is tested: whether the cis/trans decomposition shape of the cross-species analysis reproduces (trans
##     exceeding cis, cis share below 0.5, a negative cross-replicate compensation) in another laboratory,
##     another strain pair, and the opposite temperature direction.
##   What is not tested: time dynamics. This design has two levels and no time course, so "cis front-loaded,
##     trans accumulating" cannot be examined here and cannot be inferred from this design.
##   Expected: cis share < 0.5 genome-wide and per stratum; trn_tru > cis_tru; cross-replicate cor(cis, trans)
##     negative and smaller in magnitude than the same-replicate value; the Hsf1_core to other_induced cis-share
##     ratio with an interval covering 1, as in the cross-species analysis.
##   Failable: cis share >= 0.5, or a non-negative cross-replicate correlation. Either is reported as the cold
##     response in this species pair decomposing differently from the Fay data.
##   MINCPM, the filter, and replicate pairing are fixed before analysis; a third
##     decomposition is not selected after observing results.
##   Cross-mapping leakage, per species, larger of the two: lambda < 0.10 use as is; 0.10 <= lambda < 0.20 report
##     with the caveat attached to every cis estimate; lambda >= 0.20 stop and retain the analysis that the decomposition
##     cannot be done on this data. Diagnostics include each parent's own mapping rate against its own reference, and if those two
##     differ by more than 3 percentage points the direction of the asymmetry is reported.
##   Baseline diagnostic: the genome-wide median allelic log ratio in the hybrid at 30 C, and the same at 12 C.
##     A systematic mapping asymmetry appears as a departure from zero; the estimand is the difference of the
##     two, so the diagnostic states what the constant is and shows it cancels.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
gdir <- file.path(root, "data/03_external/fay2023_genomes")
ext  <- file.path(root, "data/03_external/hov2020")
outd <- file.path(root, "data/02_processed/hov2020"); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

REFS <- list(Sc = "DBVPG6765", Su = "YJF1449")

read_fasta <- function(p) {
  x <- readLines(p, warn = FALSE)
  h <- grep("^>", x)
  nm <- sub("^>\\s*", "", sub("\\s.*$", "", x[h]))
  st <- h + 1L; en <- c(h[-1] - 1L, length(x))
  setNames(vapply(seq_along(h), function(i) paste0(x[st[i]:en[i]], collapse = ""), ""), nm)
}
revcomp <- function(s) chartr("ACGTacgtNn", "TGCAtgcaNn", paste(rev(strsplit(s, NULL)[[1]]), collapse = ""))

## ------------------------------------------------------------------ CDS extraction ----
## One record per gene, built from its CDS blocks in genomic order and reverse-complemented on the minus strand.
## Naming is <systematic>|<allele> so that a salmon quantification can be split by allele with one string
## operation and joined to the companion's HOG table by the other half.
FA <- list(); N <- list()
for (al in names(REFS)) {
  strain <- REFS[[al]]
  gf <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(gdir, paste0(strain, ".final.gff")))),
              sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
  setnames(gf, c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
  ## The GFF nests attributes over three levels: `gene` carries Gene_Name, `transcript` points at its gene
  ## through Parent, and `CDS` points at its transcript. A CDS line therefore has no gene name of its own and
  ## the map takes two hops. Note: Gene_Name looked up on the CDS lines finds nothing and gives an empty
  ## FASTA, which a consistency-only size assertion accepts (0 == 0 + 0); the assertions below have a floor.
  gg <- gf[type == "gene" & grepl("Gene_Name=", attr, fixed = TRUE)]
  gg[, gid  := sub("^.*ID=([^;]+).*$", "\\1", attr)]
  gg[, sysn := sub("_.*$", "", sub("^.*Gene_Name=([^;]+).*$", "\\1", attr))]
  gg <- gg[grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", sysn)]
  tr <- gf[type == "transcript"]
  tr[, tid := sub("^.*ID=([^;]+).*$", "\\1", attr)]
  tr[, gid := sub("^.*Parent=([^;]+).*$", "\\1", attr)]
  t2s <- merge(tr[, .(tid, gid)], gg[, .(gid, sysn)], by = "gid")
  ## Some systematic names appear at more than one locus in these assemblies (YAL067C sits on three
  ## chromosomes here, being a subtelomeric SEO1/PAU family member), so they are multiple genomic copies given
  ## the same best-hit name, not isoforms of one gene. Allele-specific quantification needs exactly one sequence
  ## per gene per allele, and a read landing in a three-copy family cannot be assigned to one of them, so those
  ## names are dropped. This removes the same kind of gene the companion's one-to-one orthologue filter removes,
  ## and it is reported rather than applied silently.
  multi <- t2s[, .(n_loci = uniqueN(gid)), by = sysn][n_loci > 1L]
  say("  systematic names at more than one locus, dropped: %d (%d transcripts)",
      nrow(multi), t2s[sysn %in% multi$sysn, uniqueN(tid)])
  t2s <- t2s[!sysn %in% multi$sysn]
  stopifnot(t2s[, uniqueN(tid), by = sysn][, max(V1)] == 1L)
  cds <- gf[type == "CDS"]
  cds[, tid := sub("^.*Parent=([^;]+).*$", "\\1", attr)]
  cds <- merge(cds, t2s[, .(tid, sysn)], by = "tid")
  stopifnot(nrow(cds) > 4000L, uniqueN(cds$sysn) > 4000L)
  fa <- read_fasta(file.path(gdir, paste0(strain, ".fasta")))
  setorder(cds, sysn, chr, start)
  say("%s (%s): %d CDS blocks over %d genes on %d sequences",
      al, strain, nrow(cds), uniqueN(cds$sysn), uniqueN(cds$chr))
  seqs <- cds[, {
    ch <- chr[1]
    if (!ch %in% names(fa)) list(s = NA_character_) else {
      pieces <- substring(fa[[ch]], start, end)
      s <- paste0(pieces, collapse = "")
      if (strand[1] == "-") s <- revcomp(paste0(rev(pieces), collapse = ""))
      list(s = s)
    }
  }, by = sysn]
  seqs <- seqs[!is.na(s) & nchar(s) >= 150L]
  say("  usable transcripts (>= 150 nt): %d ; median length %d", nrow(seqs), median(nchar(seqs$s)))
  FA[[al]] <- seqs[, .(name = paste0(sysn, "|", al), s)]
  N[[al]] <- seqs$sysn
}
shared <- intersect(N$Sc, N$Su)
say("\ngenes present in both references: %d (Sc only %d, Su only %d)",
    length(shared), length(setdiff(N$Sc, N$Su)), length(setdiff(N$Su, N$Sc)))

ALL <- rbindlist(FA)
## Assert a floor as well as consistency: consistency alone is satisfied by an empty file.
stopifnot(!any(duplicated(ALL$name)), nrow(ALL) == length(N$Sc) + length(N$Su),
          length(N$Sc) > 4000L, length(N$Su) > 4000L, length(shared) > 4000L)
fp <- file.path(outd, "combined_cds.fa")
con <- file(fp, "wt")
for (i in seq_len(nrow(ALL))) { writeLines(paste0(">", ALL$name[i]), con); writeLines(ALL$s[i], con) }
close(con)
say("wrote %s: %d transcripts, %.1f MB", basename(fp), nrow(ALL), file.size(fp) / 1e6)

## Actual CDS lengths, from the extracted sequences. These are needed downstream because an allelic ratio in
## TPM is only a count ratio when the two alleles have the same length, and salmon's EffectiveLength is not the
## sequence length: it folds in the library's fragment-length distribution and the transcript's mappability, so
## it varies between libraries for the same gene. Used as a length criterion it excludes 41 to 61% of genes,
## with a count that moves from 2140 to 3148 across six libraries of the same hybrid, so it does not measure
## sequence length.
LEN <- merge(data.table(sys = N$Sc, len_Sc = nchar(FA$Sc$s)),
             data.table(sys = N$Su, len_Su = nchar(FA$Su$s)), by = "sys")
LEN[, `:=`(len_diff = abs(len_Su - len_Sc), len_frac = abs(len_Su - len_Sc) / pmax(len_Sc, len_Su))]
stopifnot(nrow(LEN) == length(shared))
fwrite(LEN, file.path(outd, "cds_lengths.tsv"), sep = "\t")
say("\nCDS length agreement between the two alleles, from the sequences themselves:")
for (thr in c(0.01, 0.02, 0.05, 0.10)) say("  differing by >= %4.0f%%: %4d of %d genes (%.1f%%)",
    100 * thr, LEN[len_frac >= thr, .N], nrow(LEN), 100 * LEN[len_frac >= thr, .N] / nrow(LEN))
say("  median absolute difference %d nt, median relative %.3f%%",
    as.integer(median(LEN$len_diff)), 100 * median(LEN$len_frac))

fwrite(data.table(sys = shared), file.path(outd, "shared_genes.tsv"), sep = "\t")
say("wrote shared_genes.tsv (%d genes usable for allelic comparison)", length(shared))

## The design table, taken from the ENA report rather than from run order, which does not follow the design.
rep <- fread(file.path(ext, "ena_filereport.tsv"))
D <- rep[, .(run = run_accession, alias = sample_alias)]
D[, c("geno", "temp", "rep") := tstrsplit(alias, "_", fixed = TRUE)]
D[, temp := as.integer(sub("C$", "", temp))]
stopifnot(nrow(D) == 18L, uniqueN(D$geno) == 3L, uniqueN(D$temp) == 2L,
          D[, .N, by = .(geno, temp)][, all(N == 3L)])
fwrite(D, file.path(outd, "design.tsv"), sep = "\t")
say("wrote design.tsv: %d runs, %d genotypes x %d temperatures x 3 replicates",
    nrow(D), uniqueN(D$geno), uniqueN(D$temp))
