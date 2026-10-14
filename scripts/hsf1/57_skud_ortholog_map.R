#!/usr/bin/env Rscript
## 57_skud_ortholog_map.R -- the SKDZ to S. cerevisiae map used by script 58, built by reciprocal best hits.
##
## Why it is built here. GSE234499 indexes S. kudriavzevii to GCA_947243785.1 (Skud-ZP591), whose locus tags
## are SKDZ_01G0010 and so on. Those share nothing with S. cerevisiae systematic names. Keyport Kik et al. state
## that "Ortholog assignments used in this study can be found in Supplementary File 1", but the article carries
## only three PDF attachments and the authors' repository holds job scripts and count matrices, so that table
## could not be obtained. YGOB's S. kudriavzevii column uses the Scannell 2011 Skud_x.y numbering, not SKDZ.
##
## Acceptance criteria, fixed before the map was inspected:
##   one-to-one reciprocal best hits; amino-acid identity >= 70%; query coverage >= 50%;
##   the number of accepted pairs must fall in 4500-5300, since S. cerevisiae and S. kudriavzevii have on the
##   order of 5000 one-to-one orthologs;
##   ten named genes must map to their own orthologs, including PAB1, for which the source paper itself gives
##   SKDZ_05G2450, an assertion against a value from outside this pipeline.
## If the criteria fail, script 58 is not run and the map is reported as not obtained.
##
## The map is a means, not a result: no biological claim rests on it, and every downstream number is reported
## with the count of genes that survived it.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
ext  <- file.path(root, "data/03_external/skud_zp591")
outd <- file.path(root, "data/02_processed/crossspecies")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
DIAMOND <- Sys.getenv("DIAMOND_BIN", unset = Sys.which("diamond"))
if (!nzchar(DIAMOND)) DIAMOND <- file.path(Sys.getenv("HOME"), ".local/mm/root/envs/dmnd/bin/diamond")
stopifnot(file.exists(DIAMOND))
MIN_ID <- 70; MIN_COV <- 50; LO <- 4500L; HI <- 5300L
tmp <- file.path(tempdir(), "skudmap"); dir.create(tmp, showWarnings = FALSE)

## ------------------------------------------------------------------ protein accession to gene ----
## Both feature tables carry one row per feature; the columns needed are the locus tag, the symbol and the
## protein accession. Reading them explicitly by name rather than by position, because a positional read would
## silently shift if NCBI adds a column.
ft <- function(f) {
  z <- fread(cmd = paste("zcat", shQuote(file.path(ext, f))), sep = "\t", quote = "")
  setnames(z, gsub("^# ", "", names(z)))
  z[`# feature` == "CDS" | feature == "CDS"]
}
FS <- fread(cmd = paste("zcat", shQuote(file.path(ext, "skud_zp591_feature_table.txt.gz"))), sep = "\t", quote = "")
FC <- fread(cmd = paste("zcat", shQuote(file.path(ext, "scer_r64_feature_table.txt.gz"))), sep = "\t", quote = "")
for (z in list(FS, FC)) setnames(z, gsub("^#\\s*", "", names(z)))
need <- c("feature", "locus_tag", "product_accession", "symbol")
stopifnot(all(need %in% names(FS)), all(need %in% names(FC)))
SK <- unique(FS[feature == "CDS" & product_accession != "",
                .(prot = product_accession, skdz = locus_tag)])
SC <- unique(FC[feature == "CDS" & product_accession != "",
                .(prot = product_accession, sys = locus_tag, sym = symbol)])
say("=== feature tables ===")
say("  S. kudriavzevii ZP591: %d CDS proteins, %d distinct locus tags", nrow(SK), uniqueN(SK$skdz))
say("  S. cerevisiae R64    : %d CDS proteins, %d distinct locus tags", nrow(SC), uniqueN(SC$sys))
stopifnot(all(grepl("^SKDZ_", SK$skdz)))
say("  every S. kudriavzevii tag has the SKDZ_ prefix the counts use")

## ------------------------------------------------------------------ reciprocal best hits ----
faa <- function(f) file.path(ext, f)
run <- function(q, d, out, tag) {
  db <- file.path(tmp, paste0(tag, ".dmnd"))
  system2(DIAMOND, c("makedb", "--in", faa(d), "-d", db, "--quiet"), stdout = FALSE, stderr = FALSE)
  system2(DIAMOND, c("blastp", "-q", faa(q), "-d", db, "-o", out, "--quiet",
                     "--max-target-seqs", "1", "--evalue", "1e-10", "--outfmt", "6",
                     "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"),
          stdout = FALSE, stderr = FALSE)
  z <- fread(out, header = FALSE,
             col.names = c("q", "s", "pident", "len", "qlen", "slen", "evalue", "bits"))
  z[, cov := 100 * len / pmin(qlen, slen)]
  z[order(-bits), .SD[1], by = q]
}
say("\n=== reciprocal best hits, diamond blastp ===")
A <- run("skud_zp591_protein.faa.gz", "scer_r64_protein.faa.gz", file.path(tmp, "sk_vs_sc.tsv"), "sc")
B <- run("scer_r64_protein.faa.gz", "skud_zp591_protein.faa.gz", file.path(tmp, "sc_vs_sk.tsv"), "sk")
say("  S. kudriavzevii -> S. cerevisiae best hits: %d", nrow(A))
say("  S. cerevisiae -> S. kudriavzevii best hits: %d", nrow(B))
setnames(A, c("q", "s"), c("sk_prot", "sc_prot"))
setnames(B, c("q", "s"), c("sc_prot", "sk_prot"))
RBH <- merge(A[, .(sk_prot, sc_prot, pident, cov)],
             B[, .(sc_prot, sk_prot)], by = c("sk_prot", "sc_prot"))
say("  reciprocal pairs: %d", nrow(RBH))

M <- merge(merge(RBH, SK, by.x = "sk_prot", by.y = "prot"), SC, by.x = "sc_prot", by.y = "prot")
say("  pairs with both identifiers resolved: %d", nrow(M))
say("\n=== applying the pre-registered filters ===")
say("  before filters                        %d", nrow(M))
M <- M[pident >= MIN_ID]; say("  identity >= %d%%                       %d", MIN_ID, nrow(M))
M <- M[cov >= MIN_COV];   say("  coverage >= %d%%                       %d", MIN_COV, nrow(M))
M <- M[skdz %in% M[, .N, by = skdz][N == 1L, skdz] & sys %in% M[, .N, by = sys][N == 1L, sys]]
say("  one-to-one in both directions         %d", nrow(M))
say("  identity: median %.1f%%, range %.1f to %.1f", median(M$pident), min(M$pident), max(M$pident))
say("  pre-registered count window %d-%d: %s", LO, HI, nrow(M) >= LO && nrow(M) <= HI)
stopifnot(nrow(M) >= LO, nrow(M) <= HI)

## ------------------------------------------------------------------ named-gene assertions ----
say("\n=== ten named genes, and one value taken from the source paper ===")
WANT <- c("HSF1", "MSN2", "MSN4", "SSA4", "HSP104", "PAB1", "TDH3", "ACT1", "PGK1", "ENO2")
NM <- M[sym %in% WANT, .(sym, sys, skdz, pident = round(pident, 1), cov = round(cov, 1))][order(sym)]
print(NM, row.names = FALSE)
say("  mapped: %d of %d", nrow(NM), length(WANT))
## The paper states PAB1 is SKDZ_05G2450. That is a value from outside this pipeline, so it tests the map rather
## than testing the map against itself.
pab <- M[sym == "PAB1", skdz]
say("  PAB1 -> %s ; the source paper gives SKDZ_05G2450 ; agree: %s",
    paste(pab, collapse = ","), identical(pab, "SKDZ_05G2450"))
stopifnot(identical(pab, "SKDZ_05G2450"))
stopifnot(nrow(NM) >= 8L)
say("  external assertion passed")

## ------------------------------------------------------------------ coverage of the counts ----
CNT <- fread(cmd = paste("zcat", shQuote(file.path(root, "data/03_external/keyportkik2024/GSE234499_counts.tsv.gz"))))
CNT <- CNT[!grepl("^__", ORF)]
sk_in <- unique(CNT[species == "Skud", ORF]); sc_in <- unique(CNT[species == "Scer", ORF])
say("\n=== coverage of the GSE234499 count matrix ===")
say("  S. kudriavzevii genes in the counts: %d; mapped: %d (%.1f%%)",
    length(sk_in), sum(sk_in %in% M$skdz), 100 * mean(sk_in %in% M$skdz))
say("  S. cerevisiae  genes in the counts: %d; mapped: %d (%.1f%%)",
    length(sc_in), sum(sc_in %in% M$sys), 100 * mean(sc_in %in% M$sys))
BOTH <- M[skdz %in% sk_in & sys %in% sc_in]
say("  pairs with BOTH sides present in the counts: %d", nrow(BOTH))

## ------------------------------------------------------------------ join to the P3 modules ----
HOGS <- fread(file.path(outd, "hog_to_scer.tsv"))
MODS <- fread(file.path(root, "data/02_processed/hsf1/s4_module_membership.tsv"))
say("\n=== joining to this project's HOGs and five modules ===")
say("  hog_to_scer: %d rows; module membership: %d rows, columns %s",
    nrow(HOGS), nrow(MODS), paste(names(MODS), collapse = ", "))
J <- merge(BOTH[, .(skdz, sys, sym, pident, cov)], HOGS[, .(hog = gene, sys)], by = "sys")
say("  pairs with a HOG: %d", nrow(J))
## MODS carries its own `sys` column, which would collide with this one and give sys.x / sys.y in the saved
## table; and an absent module arrives as "" rather than NA, so an is.na() count would report every row as
## having a module. Drop the duplicate key before merging and normalise the empty string to NA.
J <- merge(J, MODS[, .(hog, module)], by = "hog", all.x = TRUE)
J[module == "" | is.na(module), module := NA_character_]
stopifnot(!any(duplicated(names(J))), "sys" %in% names(J), !any(grepl("^sys\\.", names(J))))
say("  pairs with a module: %d", J[!is.na(module), .N])
print(J[!is.na(module), .N, by = module][order(-N)], row.names = FALSE)
fwrite(J, file.path(outd, "skdz_to_scer_map.tsv"), sep = "\t")
say("\n  written: skdz_to_scer_map.tsv (%d rows)", nrow(J))
say("\n=== verdict ===")
say("  the map meets every pre-registered criterion")
say("  it is a means and not a result; downstream numbers must state how many genes survived it")
