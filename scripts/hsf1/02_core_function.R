# 02_core_function.R -- what biological processes does the Hsf1-associated core
# consist of?
#
# What this analysis can and cannot show.
#
# 01_epistasis_decomposition.R produces two gene sets from the Mühlhofer data (GSE179258): a unanimous
# Hsf1-associated core (111 genes) and a joint-requirement candidate stratum (51 genes).
# 04_divergence_by_dependency.R projects the orthologs of these sets onto Fay et al. 2023
# and asks whether cross-species divergence concentrates in them. Before that projection
# is meaningful the biological content of the sets must be established -- otherwise a divergence
# result cannot be interpreted.
#
# Caveat: one of the contrasts below cannot fail.
#
# Testing the Hsf1 core against the whole tested universe and finding "protein folding"
# is very close to a tautology. The core was defined as genes whose heat induction
# requires Hsf1; Hsf1 is textbook-known as the chaperone regulator; and the GO
# annotations for those chaperones are curated from the same literature that established
# Hsf1 as their regulator, including Solis & Denic 2016 itself, which
# 01_epistasis_decomposition.R already uses as a positive control. So that contrast
# re-uses its own evidence and is used only as a labelled sanity check and
# not interpreted as a finding.
#
# The contrast that carries information is the second one: Hsf1-associated core versus
# the other heat-induced genes. Both arms are heat-induced, so "induced by heat" cancels,
# and what remains is what distinguishes Hsf1-dependent induction from Hsf1-independent
# induction. That is the contrast the cross-species test needs, and its answer is not
# known in advance.
#
# Definitions come from external sources: GO Consortium SGD GAF and go-basic.obo,
# both downloaded by fetch_inputs.sh with URL, date and sha256 recorded in
# manifests/EXTERNAL_INPUTS.tsv. No hand-built gene lists, no self-made
# categories. Annotations are propagated to ancestors through is_a/part_of, because a
# gene annotated to a specific term is by definition also annotated to its parents; a
# hypergeometric test on direct annotations only would under-count every general term.
#
# Statistics. Hypergeometric (one-sided over-representation) per GO term, BH within
# aspect. This is an enrichment test, not a correction to any
# differential-expression estimate, and it is run on gene-set membership, not on
# log2FCs. Terms are restricted to a testable size range in the relevant background;
# terms outside it are reported as untestable rather than as negative.
#
# Usage: Rscript scripts/hsf1/02_core_function.R [--anchor=37_30]

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
proc <- file.path(root, "data/02_processed")
outd <- file.path(proc, "hsf1")
ref  <- file.path(root, "data/03_reference")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

ALPHA    <- 0.05
MIN_TERM <- 5      # below this the hypergeometric cannot reach significance at n~100
MAX_TERM <- 500    # above this the term is too general to be informative
say <- function(...) cat(sprintf(...), "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
h <- grep("^--anchor=", args, value = TRUE)
ANCHOR <- if (length(h)) sub("^--anchor=", "", h[1]) else "37_30"

say("=== 02_core_function | %s | anchor %s C ===", format(Sys.time()), ANCHOR)
say("S. cerevisiae only. Enrichment test, not a DE correction.")

## ------------------------------------------------------- source versions ---
man <- file.path(ref, "../../manifests/EXTERNAL_INPUTS.tsv")
if (!file.exists(man)) stop("manifests/EXTERNAL_INPUTS.tsv is missing; run ./fetch_inputs.sh")
MAN <- fread(man)
for (f in c("data/03_reference/sgd.gaf.gz", "data/03_reference/go-basic.obo")) {
  r <- MAN[file == f]
  if (!nrow(r)) stop("not in manifest, refusing to proceed: ", f)
  say("locked: %-32s sha256 %s  fetched %s", basename(f), substr(r$sha256[1], 1, 16), r$fetched_utc[1])
}
OBO_VER <- sub("^data-version: ", "",
               grep("^data-version:", readLines(file.path(ref, "go-basic.obo"), n = 40), value = TRUE)[1])
say("ontology release: %s", OBO_VER)

## ---------------------------------------------------------------- name bridge -----
gtf  <- file.path(ref, "Saccharomyces_cerevisiae.R64-1-1.75.gtf")
gl   <- grep("\tgene\t", readLines(gtf, warn = FALSE), value = TRUE)
sysn <- sub('.*gene_id "([^"]+)".*', "\\1", gl)
stdn <- ifelse(grepl('gene_name "', gl), sub('.*gene_name "([^"]+)".*', "\\1", gl), NA_character_)
nmap <- unique(data.table(sys = sysn, std = stdn))
lut  <- unique(rbind(nmap[, .(key = toupper(sys), sys)], nmap[!is.na(std), .(key = toupper(std), sys)]))
lut  <- lut[!duplicated(key)]

## ------------------------------------------------------------------- ontology -----
## Parse go-basic.obo. Only is_a and part_of are followed: go-basic is the filtered
## release that guarantees these two relations form a DAG, so upward closure terminates.
say("\nparsing ontology ...")
ob <- readLines(file.path(ref, "go-basic.obo"), warn = FALSE)
starts <- which(ob == "[Term]")
## vectorised block assignment: every line gets the index of the [Term] it belongs to
blk <- cumsum(ob == "[Term]")
keep <- blk > 0
D <- data.table(blk = blk[keep], line = ob[keep])
D <- D[line != "[Term]" & line != ""]
D[, tag := sub(":.*$", "", line)]
D[, val := sub("^[^:]+: ", "", line)]

ids  <- D[tag == "id",        .(blk, id  = sub(" .*$", "", val))]
nms  <- D[tag == "name",      .(blk, name = val)]
nsp  <- D[tag == "namespace", .(blk, ns   = val)]
obso <- D[tag == "is_obsolete" & val == "true", unique(blk)]
TERM <- Reduce(function(x, y) merge(x, y, by = "blk", all.x = TRUE), list(ids, nms, nsp))
TERM <- TERM[!blk %in% obso]
say("  terms: %d non-obsolete", nrow(TERM))

## edges: child -> parent
e1 <- D[tag == "is_a", .(blk, parent = sub(" .*$", "", val))]
e2 <- D[tag == "relationship" & grepl("^part_of ", val), .(blk, parent = sub("^part_of ([^ ]+).*$", "\\1", val))]
ED <- merge(rbind(e1, e2), ids, by = "blk")[, .(child = id, parent)]
ED <- ED[child %in% TERM$id & parent %in% TERM$id]
say("  edges (is_a + part_of): %d", nrow(ED))

## ----------------------------------------------------------------- annotations ----
## SGD GAF 2.2. Column 4 is the qualifier: rows carrying NOT are negative annotations
## and must be dropped, otherwise a gene explicitly stated not to do something would be
## counted as doing it. ND (no biological data) is dropped for the same reason.
say("\nparsing annotations ...")
## The GAF carries a variable-length block of `!` comment lines before the data, so it
## has to be stripped in the pipe; fread's own header detection cannot handle it.
GAF <- fread(cmd = paste("zcat", shQuote(file.path(ref, "sgd.gaf.gz")), "| grep -v '^!'"),
             sep = "\t", quote = "", header = FALSE, fill = TRUE,
             select = c(3, 4, 5, 7, 9, 11),
             col.names = c("symbol", "qualifier", "go", "evidence", "aspect", "synonyms"))
n0 <- nrow(GAF)
GAF <- GAF[!grepl("\\bNOT\\b", qualifier) & evidence != "ND" & go != ""]
say("  GAF rows: %d, kept after dropping NOT and ND: %d", n0, nrow(GAF))

## Bridge GAF identifiers to our systematic IDs. SGD puts the systematic name in the
## synonym column, so try the symbol first and then each synonym.
GAF[, sys := lut[match(toupper(symbol), key), sys]]
un <- GAF[is.na(sys)]
if (nrow(un)) {
  sp <- strsplit(un$synonyms, "|", fixed = TRUE)
  cand <- vapply(sp, function(v) {
    if (!length(v)) return(NA_character_)
    m <- lut[match(toupper(v), key), sys]
    m <- m[!is.na(m)]
    if (length(m)) m[1] else NA_character_
  }, character(1))
  GAF[is.na(sys), sys := cand]
}
say("  mapped to systematic IDs: %d of %d rows (%.1f%%), %d distinct genes",
    GAF[!is.na(sys), .N], nrow(GAF), 100 * GAF[, mean(!is.na(sys))], GAF[!is.na(sys), uniqueN(sys)])
A <- unique(GAF[!is.na(sys) & go %in% TERM$id, .(sys, go)])

## Propagate to ancestors. Iterate child->parent joins until no new (gene, term) pair
## appears; that is the ancestor closure of the annotation set, which is what the
## ontology's true-path rule requires.
say("  propagating to ancestors ...")
frontier <- copy(A); it <- 0L
repeat {
  it <- it + 1L
  nxt <- merge(frontier, ED, by.x = "go", by.y = "child", allow.cartesian = TRUE)[
    , .(sys, go = parent)]
  nxt <- unique(nxt)
  new <- fsetdiff(nxt, A)
  if (!nrow(new)) break
  A <- unique(rbind(A, new)); frontier <- new
  if (it > 50L) stop("closure did not terminate -- ontology is not acyclic as assumed")
}
say("  annotations: %d direct+propagated pairs over %d genes and %d terms (%d iterations)",
    nrow(A), A[, uniqueN(sys)], A[, uniqueN(go)], it)
A <- merge(A, TERM[, .(go = id, name, ns)], by = "go")

## ------------------------------------------------------------------ gene sets -----
CL   <- fread(file.path(outd, sprintf("hsf1_epistasis_classes_%s.tsv", ANCHOR)))[method == "limma"]
CORE <- fread(file.path(outd, sprintf("hsf1_core_set_%s.tsv", ANCHOR)))
UNIV      <- CL$gene
core      <- CORE[set == "Hsf1_core_unanimous", gene]
joint     <- CORE[set == "Hsf1_joint_candidate", gene]
induced   <- CL[class != "not_induced", gene]
other_ind <- setdiff(induced, c(core, joint))
say("\ngene sets: universe %d | induced %d | Hsf1 core %d | joint candidates %d | other induced %d",
    length(UNIV), length(induced), length(core), length(joint), length(other_ind))

## -------------------------------------------------------------------- testing -----
## Hypergeometric over-representation of `set` within `background`. BH within aspect.
enrich <- function(set, background, label) {
  set <- intersect(set, background)
  AB  <- A[sys %in% background]
  ## term sizes in this background, not genome-wide: the background defines the draw
  sz  <- AB[, .(m = uniqueN(sys)), by = .(go, name, ns)]
  hit <- AB[sys %in% set, .(k = uniqueN(sys)), by = go]
  T2  <- merge(sz, hit, by = "go", all.x = TRUE)
  T2[is.na(k), k := 0L]
  n_bg <- length(background); n_set <- length(set)
  T2[, testable := m >= MIN_TERM & m <= MAX_TERM]
  T2[, expected := m * n_set / n_bg]
  T2[, odds := (k / max(n_set - k, 1e-9)) / (pmax(m - k, 1e-9) / max(n_bg - n_set - (m - k), 1e-9))]
  ## P(X >= k) with X ~ Hypergeometric(m successes in population n_bg, n_set draws)
  T2[, p := phyper(k - 1L, m, n_bg - m, n_set, lower.tail = FALSE)]
  T2[testable == TRUE, padj := p.adjust(p, method = "BH"), by = ns]
  T2[, `:=`(set = label, n_set = n_set, n_bg = n_bg)]
  T2[order(p)]
}

## ------------------------------------------------------ robustness of each term ---
## Two ways a hypergeometric enrichment can be significant without being a real
## statement about a biological category:
##   (1) one gene carries it -- leave-one-out exposes this;
##   (2) the hits are members of one expanded paralog family, so they are one
##       evolutionary event counted several times. The hypergeometric assumes genes are
##       exchangeable draws, which paralogs are not.
## (2) applies here: "fungal-type cell wall" reaches padj 4.1e-03 with 12 hits, of which
## TIP1, TIR1, TIR2, TIR3 and TIR4 are one mannoprotein family. Collapsing them to a
## single representative moves the raw p from 4.8e-05 to 6.2e-03, which does not survive
## BH over the testable terms in that aspect.
## Family membership is approximated by the gene-name stem (standard name with trailing
## digits removed). This is a crude proxy: it correctly groups TIR1-TIR4 and SSA1-SSA4,
## but it over-collapses HSP26/HSP42/HSP82/HSP104, which are not one family. Over-
## collapsing is conservative -- it can only make a term look weaker -- so a term that
## survives stem collapse is robust to this check, while a term that fails it is marked as not established.
## Two levels of grouping are used, and they disagree, so both are reported.
##
## Level 1, automatic: gene-name stem (standard name minus trailing digits). Cheap and
## reproducible, but wrong in both directions -- it splits TIP1 from TIR1-TIR4 although
## they are one family, and merges HSP26/HSP42/HSP82/HSP104 although they are not.
##
## Level 2, externally verified families. Required because level 1's verdict on the
## cell-wall terms depends entirely on whether TIP1 is grouped with the TIRs. The single
## family entered here is verified, not assumed:
##   Zhou X, Fay JC et al. -- anaerobically induced cell wall mannoproteins (anCWMPs)
##   TIP1, TIR1, TIR2, TIR3, TIR4, DAN1, DAN4 "likely derive evolutionarily from a single
##   ancestral anCWMP locus, which was duplicated and translocated to new genomic
##   contexts several times", Mol Biol Evol 2022, PMID 36134526, doi 10.1093/molbev/msac202.
##   The same family is also known for anaerobic and cold-shock induction rather than heat
##   (PMID 11292809, doi 10.1128/JB.183.9.2881-2887.2001), which is why its appearance in a
##   heat-induced Hsf1-dependent set is worth stating at gene level.
## Because these genes descend from one ancestral locus, their joint presence in the core
## is one evolutionary event, and a hypergeometric test that counts it five times is
## overstating the evidence for a cell-wall category.
FAMILIES <- list(
  anCWMP = c("TIP1", "TIR1", "TIR2", "TIR3", "TIR4", "DAN1", "DAN4")   # PMID 36134526
)
stem_of <- local({
  s  <- nmap[!is.na(std)]
  st <- toupper(sub("[0-9]+$", "", s$std))
  st[nchar(st) < 2] <- toupper(s$std[nchar(st) < 2])
  names(st) <- s$sys
  ## externally verified families override the stem heuristic
  for (fam in names(FAMILIES)) {
    ids <- lut[match(FAMILIES[[fam]], key), sys]
    ids <- ids[!is.na(ids)]
    st[names(st) %in% ids] <- paste0("FAM_", fam)
  }
  st
})

robust <- function(T2, set, background, topn = 12) {
  AB <- A[sys %in% background]
  s0 <- intersect(set, background)
  out <- T2[testable == TRUE & padj < ALPHA][order(p)][seq_len(min(topn, .N))]
  if (!nrow(out)) return(out)
  res <- rbindlist(lapply(seq_len(nrow(out)), function(i) {
    tg  <- AB[go == out$go[i], sys]
    hits <- intersect(s0, tg)
    ## leave-one-out
    loo <- max(vapply(hits, function(gi) {
      s <- setdiff(s0, gi); m <- length(tg)
      phyper(length(intersect(s, tg)) - 1L, m, length(background) - m, length(s), lower.tail = FALSE)
    }, numeric(1)))
    ## paralog-stem collapse: keep one representative per stem, drop the rest from both
    ## the set and the background so the term size shrinks consistently
    stems <- stem_of[hits]; stems[is.na(stems)] <- hits[is.na(stems)]
    keep  <- hits[!duplicated(stems)]
    drop  <- setdiff(hits, keep)
    s2 <- setdiff(s0, drop); bg2 <- setdiff(background, drop); tg2 <- setdiff(tg, drop)
    m2 <- length(tg2)
    p2 <- phyper(length(intersect(s2, tg2)) - 1L, m2, length(bg2) - m2, length(s2), lower.tail = FALSE)
    data.table(go = out$go[i], name = out$name[i], ns = out$ns[i], k = out$k[i], m = out$m[i],
               p = out$p[i], padj = out$padj[i], loo_worst_p = loo,
               n_stems = uniqueN(stems), families = paste(sort(unique(stems[duplicated(stems)])), collapse = "/"))
  }))
  res[order(p)]
}

## Full paralog-collapsed re-run of a contrast. Whole background is collapsed to one
## representative per name stem, then the identical enrich() + BH procedure is applied.
## Comparing padj to padj is the like-for-like comparison; comparing a collapsed raw p
## to a rank-1 Bonferroni bar is a different and much stricter test that flags even the
## known-true folding signal.
collapse_rerun <- function(set, background, label) {
  st <- stem_of[background]; st[is.na(st)] <- background[is.na(st)]
  keep <- background[!duplicated(st)]
  say("  paralog collapse: background %d -> %d genes (%d stems merged away)",
      length(background), length(keep), length(background) - length(keep))
  enrich(intersect(set, keep), keep, paste0(label, "_paralog_collapsed"))
}

report_robust <- function(R, EC, label) {
  if (!nrow(R)) return(invisible(NULL))
  R2 <- merge(R, EC[, .(go, k_coll = k, m_coll = m, p_coll = p, padj_coll = padj)], by = "go", all.x = TRUE)[order(p)]
  R2[, verdict := fifelse(is.na(padj_coll), "untestable after collapse",
                  fifelse(padj_coll < ALPHA, "survives", "FLAGGED"))]
  say("\nrobustness of the significant terms in %s", label)
  say("(same enrich+BH procedure re-run on a paralog-collapsed background):")
  say("%-38s %4s %6s %11s %7s %11s %-26s", "term", "k", "stems", "worst LOO p", "k_coll", "padj_coll", "verdict")
  for (i in seq_len(nrow(R2)))
    say("%-38s %4d %6d %11.2e %7s %11s %-26s", substr(R2$name[i], 1, 38), R2$k[i], R2$n_stems[i],
        R2$loo_worst_p[i],
        if (is.na(R2$k_coll[i])) "-" else as.character(R2$k_coll[i]),
        if (is.na(R2$padj_coll[i])) "-" else sprintf("%.2e", R2$padj_coll[i]),
        R2$verdict[i])
  fl <- R2[verdict == "FLAGGED"]
  if (nrow(fl)) {
    say("Terms marked as not established lose significance once each paralog family counts once. Families")
    say("responsible: %s", paste(unique(fl$families[fl$families != ""]), collapse = "; "))
    say("These do not establish category-level enrichment; the per-gene results are unchanged.")
  }
  say("NOTE: the stem heuristic over-collapses non-families (HSP26/HSP42/HSP82/HSP104 all")
  say("become 'HSP'), so this check is conservative: 'survives' is robust to this check, while 'not established'")
  say("means 'not established', not 'refuted'.")
  invisible(R2)
}

report <- function(T2, label, topn = 12) {
  say("\n--- %s ---", label)
  say("terms in background: %d, testable (%d<=size<=%d): %d, significant at BH<%.2f: %d",
      nrow(T2), MIN_TERM, MAX_TERM, T2[testable == TRUE, .N], ALPHA,
      T2[testable == TRUE & padj < ALPHA, .N])
  by_ns <- T2[testable == TRUE, .(testable = .N, sig = sum(padj < ALPHA)), by = ns]
  for (i in seq_len(nrow(by_ns)))
    say("  %-20s testable %5d  significant %4d", by_ns$ns[i], by_ns$testable[i], by_ns$sig[i])
  x <- T2[testable == TRUE & padj < ALPHA][order(p)][seq_len(min(topn, .N))]
  if (!nrow(x)) { say("  (nothing significant)"); return(invisible(NULL)) }
  say("%-42s %-14s %4s %5s %8s %7s %10s", "term", "aspect", "k", "size", "expected", "OR", "padj")
  for (i in seq_len(nrow(x)))
    say("%-42s %-14s %4d %5d %8.2f %7.1f %10.2e",
        substr(x$name[i], 1, 42), sub("^biological_", "bio_", sub("^molecular_", "mol_", sub("^cellular_", "cell_", x$ns[i]))),
        x$k[i], x$m[i], x$expected[i], x$odds[i], x$padj[i])
}

say("\n================ CONTRAST 1 (SANITY CHECK) ================")
say("Hsf1 core vs the full tested universe. This contrast partly re-uses its own")
say("evidence: the core was defined by Hsf1 dependence, and the GO annotations of the")
say("chaperones in it were curated from the same literature that established Hsf1 as")
say("their regulator (Solis & Denic 2016 included, also used as a positive control in the epistasis decomposition).")
say("A chaperone/folding signal here is expected by construction and is not interpreted")
say("as a new result. It checks that the analysis recovers a known signal.")
E1 <- enrich(core, UNIV, "core_vs_universe")
report(E1, "Hsf1 core (n=111) vs tested universe")
R1 <- robust(E1, core, UNIV)
EC1 <- collapse_rerun(core, UNIV, "core_vs_universe")
V1 <- report_robust(R1, EC1, "contrast 1")

say("\n================ CONTRAST 2 (CORE VERSUS OTHER HEAT-INDUCED GENES) ================")
say("Hsf1-associated core vs the OTHER heat-induced genes. Both arms are heat-induced,")
say("so induction cancels and what remains distinguishes Hsf1-DEPENDENT induction from")
say("Hsf1-INDEPENDENT induction. This is the relevant cross-species contrast and its answer was not")
say("known in advance.")
bg2 <- c(core, other_ind)
E2 <- enrich(core, bg2, "core_vs_other_induced")
report(E2, "Hsf1 core (n=111) vs other heat-induced genes")
RB2 <- robust(E2, core, bg2)
EC2 <- collapse_rerun(core, bg2, "core_vs_other_induced")
V2 <- report_robust(RB2, EC2, "contrast 2")

say("\n================ CONTRAST 3 (the joint-requirement stratum) ================")
say("The 51 joint-requirement candidates, same background. Kept separate from the core")
say("on purpose: their Hsf1 contribution is resolved at class level but only 13 of 51")
say("are significant per gene, so they remain separate from the core.")
bg3 <- c(joint, other_ind)
E3 <- enrich(joint, bg3, "joint_vs_other_induced")
report(E3, "joint-requirement candidates (n=51) vs other heat-induced genes")
## A negative result at n=51 is interpretable only with its detection limit:
## a negative is only reportable once the test is shown able to return a
## positive at some stated effect size.
say("\ndetection limit of contrast 3 (n=%d, background %d), for a term of median testable size %d:",
    length(intersect(joint, bg3)), length(bg3), as.integer(E3[testable == TRUE, median(m)]))
{
  n3 <- length(intersect(joint, bg3)); N3 <- length(bg3)
  m3 <- as.integer(E3[testable == TRUE, median(m)])
  nt <- E3[testable == TRUE, .N, by = ns][, max(N)]
  bar <- ALPHA / nt
  ks <- 1:min(n3, m3)
  ps <- phyper(ks - 1L, m3, N3 - m3, n3, lower.tail = FALSE)
  kmin <- ks[which(ps < bar)][1]
  if (is.na(kmin)) {
    say("  no attainable k clears alpha/%d => the test CANNOT return a positive here.", nt)
    say("  The test is not feasible at this effect size and does not provide a negative result.")
  } else {
    say("  needs k >= %d of %d set members in the term (OR ~ %.1f) to clear alpha/%d = %.2e",
        kmin, n3, (kmin / (n3 - kmin)) / ((m3 - kmin) / (N3 - n3 - (m3 - kmin))), nt, bar)
    say("  observed maximum k over all testable terms: %d", E3[testable == TRUE, max(k)])
    say("  => a negative result is informative only for effects at or above that odds ratio.")
  }
}

## ------------------------------------------------- is contrast 2 just chaperones? --
## If every term significant in contrast 2 is a chaperone/folding term, then contrast 2
## carries the same near-tautological information as contrast 1 and nothing more. Check
## explicitly rather than assuming either way.
say("\n=== does contrast 2 say anything beyond 'chaperones'? ===")
CHAP <- "chaperon|protein folding|unfolded|heat shock|proteasom|ubiquitin|protein refold|misfold"
s2 <- E2[testable == TRUE & padj < ALPHA]
say("significant terms in contrast 2: %d; matching a folding/degradation keyword: %d; other: %d",
    nrow(s2), s2[grepl(CHAP, name, ignore.case = TRUE), .N],
    s2[!grepl(CHAP, name, ignore.case = TRUE), .N])
oth <- s2[!grepl(CHAP, name, ignore.case = TRUE)][order(p)]
if (nrow(oth)) {
  say("the non-folding terms, WITH their paralog-collapse verdict:")
  for (i in seq_len(min(15L, nrow(oth)))) {
    v <- V2[name == oth$name[i]]
    say("  %-40s %-20s k=%3d/%3d  OR %5.1f  padj %.2e  %s",
        substr(oth$name[i], 1, 40), sub("_", " ", oth$ns[i]), oth$k[i], oth$m[i], oth$odds[i], oth$padj[i],
        if (!nrow(v)) "(not in top set)" else if (v$verdict[1] == "survives") "SURVIVES collapse" else "FLAGGED: one paralog family")
  }
  nsurv <- sum(vapply(oth$name, function(nn) { v <- V2[name == nn]; nrow(v) > 0 && v$verdict[1] == "survives" }, logical(1)))
  say("non-folding terms surviving paralog collapse: %d of %d", nsurv, nrow(oth))
  if (nsurv == 0) {
    say("Result: beyond the expected folding signal, contrast 2 yields no")
    say("category-level enrichment robust to paralog collapse. The per-gene")
    say("observations remain valid but do not establish a category-level effect.")
  }
} else {
  say("none: contrast 2 reduces to the folding/degradation signal, i.e. it adds nothing")
  say("beyond contrast 1, so it provides no additional biological result.")
}

## ------------------------------------------------------------------- outputs ------
ALL <- rbindlist(list(E1, E2, E3, EC1, EC2), fill = TRUE)
fwrite(ALL, file.path(outd, sprintf("hsf1_core_go_enrichment_%s.tsv", ANCHOR)), sep = "\t")
fwrite(A[sys %in% UNIV], file.path(outd, "go_annotations_propagated_universe.tsv.gz"), sep = "\t")
say("\nwrote enrichment table (%d rows, 3 contrasts) and the propagated annotation set", nrow(ALL))
say("ontology release %s; GAF sha256 %s", OBO_VER, substr(MAN[file == "data/03_reference/sgd.gaf.gz", sha256[1]], 1, 16))
say("=== done | %s ===", format(Sys.time()))
