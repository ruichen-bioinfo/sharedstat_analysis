#!/usr/bin/env Rscript
## 49_p2x_full_scan.R -- compute the sequence measures on every S. cerevisiae gene, not on the
## cross-species ortholog subset.
##
## Background. The cross-species promoter set, where a promoter must exist in all four
## species, have a one-to-one orthologue, and have a length within 2x of the others. Those are the right filters
## for a four-species comparison and the wrong ones for P2, which is a single-species analysis. The cost is
## specific and checkable: of fourteen canonical Hsf1 targets, the three that fall in the `Hsf1_only` class
## include SSA1 and HSP82, and neither is in that scan; only HSP60 is. Genes with paralogues are exactly what a
## one-to-one filter removes, and the HSP families are the most paralogous Hsf1 targets there are, so the 21.3%
## HSE carriage estimated for `Hsf1_only` on that subset is expected to be an underestimate.
##
## The single-species analysis selects all genes with a systematic name. Script 10 selects genes with
##     tg <- genes[sysn %in% sysmap$sys & !is.na(sysn)]
## where `sysmap` is the one-to-one HOG universe. Here the selection is every gene carrying a systematic name.
## The reference, the 600 bp window, the clipping at the nearest upstream feature, the HSE scanner and the STRE
## count are transplanted from script 10 unchanged. S. cerevisiae only: the other three strains are irrelevant
## to P2.
##
## The R64 GTF is not used for extraction. Its coordinates are S288C's and the sequence on disk is DBVPG6765's,
## so applying one to the other would silently shift every window. It is used only to report how much of the
## annotated genome the DBVPG6765 GFF accounts for.
##
## On genes shared by the single-species and cross-species sets, every promoter and every HSE and STRE count must be identical.
## A mismatch means that the single-species set is not a strict extension, and output is not produced.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
refd <- file.path(root, "data/03_reference")
gdir <- file.path(root, "data/03_external/fay2023_genomes")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

STRE <- "CCCCT"
CANON_HSF1 <- c(SSA1 = "YAL005C", HSP82 = "YPL240C", HSP60 = "YLR259C")

## ------------------------------------------------------------------ transplant from script 10 ----
SRC <- file.path(root, "scripts/hsf1/10_promoter_hse_architecture_v2.R")
EX  <- parse(SRC)
need <- c("count_hse", "read_fasta", "revcomp", "WIN", "PRc")
have <- function() all(vapply(need, exists, TRUE, envir = globalenv()))
n_eval <- 0L
for (e in EX) {
  s <- paste(deparse(e), collapse = " ")
  if (grepl("^sink\\(|^fwrite\\(|quit\\(|^writeLines\\(", s)) next
  if (grepl("^say\\(", s)) next
  try(eval(e, envir = globalenv()), silent = TRUE)
  n_eval <- n_eval + 1L
  if (have()) break
}
say("transplant: evaluated %d of %d expressions from %s", n_eval, length(EX), basename(SRC))
stopifnot(have(), is.function(count_hse), is.function(read_fasta), WIN == 600L)
OLD <- PRc[species == "Sc", .(sys, prom_old = prom, plen_old = plen)]
say("cross-species Sc promoter set (orthologue-filtered): %d genes, window %d bp", nrow(OLD), WIN)

## ------------------------------------------------------------------ extraction, filter removed ----
STRAIN <- "DBVPG6765"
gf <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(gdir, paste0(STRAIN, ".final.gff")))),
            sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
setnames(gf, c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
gf[, gname := NA_character_]
hasn <- grepl("Gene_Name=", gf$attr, fixed = TRUE)
gf[hasn, gname := sub("^.*Gene_Name=([^;]+).*$", "\\1", attr)]
gf[, sysn := sub("_.*$", "", gname)]
genes <- gf[type == "gene"]
feats <- gf[type %in% c("gene", "tRNA", "rRNA", "ncRNA")][, .(chr, start, end)]
setorder(feats, chr, start)
fa <- read_fasta(file.path(gdir, paste0(STRAIN, ".fasta")))

## single-species selection: every gene with a systematic name, no orthology or length condition
tg <- genes[!is.na(sysn) & grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", sysn)]
tg <- tg[!duplicated(sysn)]
say("genes in the %s GFF: %d total, %d with a systematic name after de-duplication",
    STRAIN, nrow(genes), nrow(tg))

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
  out[[k]] <- data.table(sys = tg$sysn[k], prom = seqp, plen = nchar(seqp))
}
PN <- rbindlist(out)
say("promoters extracted: %d genes, median length %d, %.1f%% at the full %d bp, %.1f%% under 200 bp",
    nrow(PN), median(PN$plen), 100 * mean(PN$plen >= WIN), WIN, 100 * mean(PN$plen < 200L))

## ------------------------------------------------------------------ the three canonical targets ----
say("\n=== canonical targets absent from the orthologue-filtered set ===")
for (nm in names(CANON_HSF1)) {
  s <- CANON_HSF1[[nm]]
  say("  %-6s %-10s orthologue-filtered %-8s full-genome %s", nm, s,
      if (s %in% OLD$sys) "present" else "ABSENT", if (s %in% PN$sys) "present" else "ABSENT")
}
stopifnot(all(unlist(CANON_HSF1) %in% PN$sys))

## ------------------------------------------------------------------ fidelity check ----
say("\n=== validation: shared genes are byte-identical ===")
SH <- merge(OLD, PN, by = "sys")
say("  genes in both scans: %d", nrow(SH))
stopifnot(nrow(SH) > 4000L)
same_seq <- SH[, sum(prom_old != prom)]
same_len <- SH[, sum(plen_old != plen)]
say("  promoters differing in sequence: %d ; in length: %d", same_seq, same_len)
stopifnot(same_seq == 0L, same_len == 0L)
HO <- cbind(SH[, .(sys)], rbindlist(lapply(SH$prom_old, function(s) as.data.table(count_hse(s)))))
HN <- cbind(SH[, .(sys)], rbindlist(lapply(SH$prom,     function(s) as.data.table(count_hse(s)))))
stopifnot(identical(HO, HN))
say("  HSE counts identical on all %d shared genes. The full-genome set is a superset of the orthologue-filtered set.", nrow(SH))

## ------------------------------------------------------------------ measures on the full set ----
revcomp_s <- function(s) chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
count_fixed <- function(s, motif) {
  f <- function(m) { if (regexpr(m, s, fixed = TRUE) == -1L) 0L else length(gregexpr(m, s, fixed = TRUE)[[1]]) }
  f(motif) + f(revcomp_s(motif))
}
H <- cbind(PN[, .(sys, plen)], rbindlist(lapply(PN$prom, function(s) as.data.table(count_hse(s)))))
H[, n_stre := vapply(PN$prom, count_fixed, 0L, motif = STRE)]
H[, `:=`(hse_any = n_classic > 0 | n_gapped > 0, stre_per_kb = 1000 * n_stre / plen)]
fwrite(H, file.path(hsd, "p2x_full_scan_per_gene.tsv"), sep = "\t")

Y <- fread(file.path(refd, "yeastract2019-flat-file.tsv"), header = FALSE, quote = "", sep = "\t",
           select = 1:9, col.names = c("tf_sys", "tf_name", "tgt_sys", "tgt_name", "pmid", "date",
                                       "condition", "effect", "evidence"))
Y[, tf_name := toupper(tf_name)]
tgt <- function(tf, ev = NULL) { z <- Y[tf_name == tf]; if (!is.null(ev)) z <- z[evidence == ev]; unique(z$tgt_sys) }
MSN_D <- union(tgt("MSN2", "Direct"), tgt("MSN4", "Direct")); MSN_A <- union(tgt("MSN2"), tgt("MSN4"))
HSF_D <- tgt("HSF1", "Direct"); HSF_A <- tgt("HSF1")

CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
UN <- CL[, .(nd = uniqueN(class), cls = class[1]), by = gene][nd == 1L, .(sys = gene, cls)]
D <- merge(UN, H, by = "sys")
D[, `:=`(msn_direct = sys %in% MSN_D, msn_any = sys %in% MSN_A,
         hsf_direct = sys %in% HSF_D, hsf_any = sys %in% HSF_A)]
say("\nP2 classes with a scanned promoter: %d of %d unanimous genes (was 4374)", nrow(D), nrow(UN))

ORD <- c("Hsf1_only", "both_required", "redundant", "Msn24_only", "TF_independent", "unassigned", "not_induced")
D[, cls := factor(cls, levels = intersect(ORD, unique(cls)))]
TAB <- D[, .(n = .N, hse_pct = 100 * mean(hse_any), stre_med = as.numeric(median(n_stre)),
             stre_per_kb = mean(stre_per_kb), msn_direct_pct = 100 * mean(msn_direct),
             msn_any_pct = 100 * mean(msn_any), hsf_direct_pct = 100 * mean(hsf_direct),
             hsf_any_pct = 100 * mean(hsf_any)), by = cls]
setorder(TAB, cls)
fwrite(TAB, file.path(hsd, "p2x_full_class_sequence_curation.tsv"), sep = "\t")

say("\n=== full-genome scan, by P2 class ===")
print(TAB[, .(cls, n, hse_pct = round(hse_pct, 1), stre_med, stre_per_kb = round(stre_per_kb, 2),
              msn_direct_pct = round(msn_direct_pct, 1), msn_any_pct = round(msn_any_pct, 1),
              hsf_direct_pct = round(hsf_direct_pct, 1))], row.names = FALSE)

say("\n=== comparison with the orthologue-filtered set ===")
OT <- fread(file.path(hsd, "p2x_class_sequence_curation.tsv"))
CM <- merge(OT[, .(cls, n_old = n, hse_old = hse_pct, stre_old = stre_per_kb, msnd_old = msn_direct_pct)],
            TAB[, .(cls = as.character(cls), n_new = n, hse_new = hse_pct, stre_new = stre_per_kb,
                    msnd_new = msn_direct_pct)], by = "cls")
CM[, cls := factor(cls, levels = intersect(ORD, cls))]
setorder(CM, cls)
print(CM[, .(cls, n_old, n_new, hse_old = round(hse_old, 1), hse_new = round(hse_new, 1),
             d_hse = round(hse_new - hse_old, 1), stre_old = round(stre_old, 2),
             stre_new = round(stre_new, 2), msnd_old = round(msnd_old, 1),
             msnd_new = round(msnd_new, 1))], row.names = FALSE)

say("\n=== does the enlarged scan change the two clean conclusions? ===")
o <- TAB[order(-hse_pct)]
say("  HSE carriage ranking: %s", paste(sprintf("%s %.1f%%", o$cls, o$hse_pct), collapse = " > "))
say("  Hsf1_only still top on HSE: %s", identical(as.character(o$cls[1]), "Hsf1_only"))
s <- TAB[order(-stre_per_kb)]
say("  STRE ranking: %s", paste(sprintf("%s %.2f", s$cls, s$stre_per_kb), collapse = " > "))
say("  Msn24_only still top on STRE: %s", identical(as.character(s$cls[1]), "Msn24_only"))
r <- TAB[cls == "redundant"]; m <- TAB[cls == "Msn24_only"]; b <- TAB[cls == "both_required"]
say("\n  mixture check, 60/40 Msn2/4 to other: STRE expected %.2f, observed %.2f",
    0.6 * m$stre_per_kb + 0.4 * b$stre_per_kb, r$stre_per_kb)
say("  HSE expected %.1f%%, observed %.1f%%", 0.6 * m$hse_pct + 0.4 * b$hse_pct, r$hse_pct)

## coverage of the annotated genome, reported not assumed
G <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(refd, "Saccharomyces_cerevisiae.R64-1-1.75.gtf"))),
           sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE, select = c(3, 9),
           col.names = c("type", "attr"))
r64 <- unique(sub('^.*gene_id "([^"]+)".*$', "\\1", G[type == "gene", attr]))
r64 <- r64[grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", r64)]
say("\n  R64 annotates %d nuclear protein-coding systematic names; this scan covers %d of them (%.1f%%)",
    length(r64), length(intersect(r64, PN$sys)), 100 * length(intersect(r64, PN$sys)) / length(r64))
say("  scanned genes absent from R64 (strain-specific or renamed): %d", length(setdiff(PN$sys, r64)))
say("")
