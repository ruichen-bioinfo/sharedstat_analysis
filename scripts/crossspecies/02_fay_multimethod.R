# 02_fay_multimethod.R — does the "HSR has diverged between species" conclusion
# survive (a) a change of DE method and (b) a correctly paired model?
#
# Background. Fay et al. 2023 (GBE, PMID 37972247) assay 4 Saccharomyces species
# and their 6 pairwise hybrids under a +/-12C shift at 4 timepoints, 2 replicates.
# In the hybrids, reads are assigned to each parental allele, so a single library
# yields two count columns. Their model is
#     ~ Time * Allele          (Time a 4-level factor)
# fitted per strain x condition, with three nested DESeq2 LRTs:
#     time   : reduced ~ Allele          (response to the shift)
#     allele : reduced ~ Time            (cis divergence in overall level)
#     inter  : reduced ~ Time + Allele   (cis divergence in the response)
# `inter` is the contrast their evolutionary claim rests on. This script varies one
# thing at a time against that published baseline.
#
# Two things are varied here.
#
# (1) Method. DESeq2 LRT / edgeR quasi-likelihood F / limma-voom moderated F, on
#     one shared gene universe per group x condition so set sizes are comparable.
#     `inter` is a difference-of-differences with n=2 replicates: the lowest-powered
#     contrast in the design, and therefore where method choice should matter most.
#
# (2) Model specification. In design.csv the numeric suffix of `Name` is the library:
#     ScxSp/heat/T=0 gives Sc5 and Sp5 — the same library, sequenced once, split by
#     allele. Model A (published) treats those two rows as independent samples.
#     Model B blocks on library:
#         ~ Lib + Time * Allele      (Time main effects are aliased with Lib and drop)
#     Library then absorbs the shared RNA extraction, shift and sequencing variance,
#     and the allele and interaction terms become within-library contrasts. This is
#     the standard paired form for allele-specific expression. Model B costs residual
#     df (4 vs 8 here), so which specification wins is an empirical question, not a
#     matter of principle — that is what is quantified below.
#
#     Model B applies to hybrids only. In the `Parents` group each species is its own
#     library, so "Allele" is a between-library species contrast with nothing to pair.
#     That asymmetry is itself reported: Parents.inter and hybrid.inter are different
#     estimands (cis+trans vs cis alone, unpaired vs paired, 4 allele levels vs 2),
#     yet the published table places them in adjacent columns of one file.
#
# Not a benchmark: no method is treated as ground truth. The question is whether the
# published conclusion is stable, and where it is not.
#
# Usage: Rscript scripts/crossspecies/02_fay_multimethod.R [--groups=..] [--conditions=..]

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma); library(data.table)
})

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
ext  <- file.path(root, "data/03_external/fay2023")
outd <- file.path(root, "data/02_processed/crossspecies")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

SCRIPT_VERSION <- "02.1"
ALPHA    <- 0.05
MIN_CNT  <- 10   # shared universe: rowSums(count >= MIN_CNT) >= MIN_SMP
MIN_SMP  <- 3
SEED     <- 20260905
set.seed(SEED)

args <- commandArgs(trailingOnly = TRUE)
argv <- function(f, d) { h <- grep(paste0("^", f, "="), args, value = TRUE)
  if (length(h)) strsplit(sub(paste0("^", f, "="), "", h[1]), ",")[[1]] else d }

say <- function(...) cat(sprintf(...), "\n", sep = "")
say("=== 02_fay_multimethod v%s | %s | %s ===", SCRIPT_VERSION, format(Sys.time()), R.version.string)

## ------------------------------------------------------------------ inputs ---
dat    <- read.csv(file.path(ext, "htseq.ortho.csv"))
design <- read.csv(file.path(ext, "design.csv"))
cts <- as.matrix(dat[, -(1:14)]); rownames(cts) <- dat$HOG
stopifnot(nrow(design) == ncol(cts))

drop <- which(colnames(cts) == "Sp26")
if (length(drop) == 1) { cts <- cts[, -drop, drop = FALSE]; design <- design[-drop, , drop = FALSE] }
cts[is.na(cts)] <- 0
rownames(design) <- colnames(cts)

design$Time   <- factor(design$Time)
design$Allele <- factor(design$Allele)
## library id: the numeric suffix of Name. Sc5 and Sp5 are one library split by allele.
design$Lib    <- factor(sub("^[A-Za-z]+", "", design$Name))

## restrict to 1:1 single-copy orthologs present in all four species: the only gene
## set where an allele-level count is comparable across species without copy-number
## and mapping-length caveats
keep_hog <- dat$one == 1 & dat$all4 == 1
cts <- cts[keep_hog, , drop = FALSE]
hog_name <- dat$Name[keep_hog]
say("gene universe: %d 1:1 orthologs (of %d HOGs)", nrow(cts), nrow(dat))

## sanity check on the pairing claim
hyb <- c("ScxSp","ScxSk","ScxSu","SpxSk","SpxSu","SkxSu")
chk <- as.data.table(design)[Strain %in% hyb, .(n_allele = uniqueN(Allele), n_row = .N),
                             by = .(Strain, Condition, Lib)]
say("pairing check on hybrid libraries: %d libraries, %d with exactly 2 alleles, %d with 1",
    nrow(chk), sum(chk$n_allele == 2), sum(chk$n_allele == 1))
if (any(chk$n_allele > 2)) stop("a hybrid library carries more than two alleles")

## ------------------------------------------------------------ design matrices ---
## Full-rank model matrix, aliased columns removed by QR pivoting. Coefficient sets
## for the three nested tests are then read off the column names, so the same term
## definitions serve DESeq2, edgeR and limma.
build_mm <- function(dz, model) {
  f <- if (model == "A") ~ Time * Allele else ~ Lib + Time * Allele
  mm <- model.matrix(f, dz)
  q  <- qr(mm)
  if (q$rank < ncol(mm)) mm <- mm[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  cn    <- colnames(mm)
  inter <- which(grepl(":", cn, fixed = TRUE))
  tmain <- which(grepl("^Time", cn) & !grepl(":", cn, fixed = TRUE))
  amain <- which(grepl("^Allele", cn) & !grepl(":", cn, fixed = TRUE))
  ## the three nested tests, as coefficient sets to drop (matching Fay's reduced models)
  terms <- list(time = c(tmain, inter), allele = c(amain, inter), inter = inter)
  ## under model B the Time main effects are aliased with Lib and were dropped by the
  ## QR step, which would make `time` collapse onto `inter`; the shared response is
  ## simply not estimable once library is blocked, so the test is not offered
  if (length(tmain) == 0) terms$time <- NULL
  terms <- terms[vapply(terms, length, 1L) > 0]
  list(mm = mm, terms = terms, n_par = ncol(mm), resid_df = nrow(mm) - ncol(mm))
}

## ---------------------------------------------------------------- methods ------
## Each returns a data.table(gene, term, stat, pval, padj, eff) over the shared
## universe. `eff` is a comparable scalar effect size: for `inter` it is the
## interaction coefficient at the peak-response timepoint (15 min), which is the
## difference-of-differences the evolutionary claim is about; for `allele` the
## allele main effect; for `time` the 15 min time coefficient as well: eff_col()
## returns pick[1], and with Time a 4-level factor (levels 0/15/30/60) the first Time
## main-effect column is Time15.
## Note: the reference-level bias argument requires eff[time] and eff[inter] to be the
## same timepoint's coefficients, because that is what makes them share the reference
## measurement error with opposite sign. Column names are (Intercept), Time15, Time30,
## Time60, Alleleb, Time15:Alleleb, ..., so `time` picks Time15 and `inter` picks
## Time15:Alleleb.
eff_col <- function(cn, term) {
  has_int <- grepl(":", cn, fixed = TRUE)
  pick <- switch(term,
    inter  = which(has_int & grepl("^Time15:", cn)),
    allele = which(!has_int & grepl("^Allele", cn)),
    time   = which(!has_int & grepl("^Time", cn)))
  if (length(pick) == 0) NA_integer_ else pick[1]
}

run_deseq2 <- function(cw, dz, D) {
  dds <- suppressMessages(DESeqDataSetFromMatrix(cw, dz, design = ~ 1))
  dds <- estimateSizeFactors(dds)
  dds <- suppressMessages(estimateDispersions(dds, modelMatrix = D$mm, quiet = TRUE))
  out <- list()
  for (tm in names(D$terms)) {
    red <- D$mm[, -D$terms[[tm]], drop = FALSE]
    d <- suppressMessages(nbinomLRT(dds, full = D$mm, reduced = red, quiet = TRUE))
    r <- results(d, independentFiltering = FALSE, cooksCutoff = FALSE)
    ec <- eff_col(colnames(D$mm), tm)
    e  <- if (is.na(ec)) NA_real_ else coef(d)[, ec]
    out[[tm]] <- data.table(gene = rownames(cw), term = tm, stat = r$stat,
                            pval = r$pvalue, padj = r$padj, eff = as.numeric(e))
  }
  rbindlist(out)
}

run_edger <- function(cw, dz, D) {
  y <- calcNormFactors(DGEList(counts = cw))
  y <- estimateDisp(y, D$mm, robust = TRUE)
  fit <- glmQLFit(y, D$mm, robust = TRUE)
  out <- list()
  for (tm in names(D$terms)) {
    q  <- glmQLFTest(fit, coef = D$terms[[tm]])
    tt <- q$table
    ec <- eff_col(colnames(D$mm), tm)
    out[[tm]] <- data.table(gene = rownames(tt), term = tm, stat = tt$F,
                            pval = tt$PValue, padj = p.adjust(tt$PValue, "BH"),
                            eff = if (is.na(ec)) NA_real_ else fit$coefficients[, ec] / log(2))
  }
  rbindlist(out)
}

run_limma <- function(cw, dz, D) {
  y <- calcNormFactors(DGEList(counts = cw))
  v <- voom(y, D$mm, plot = FALSE)
  efit <- eBayes(lmFit(v, D$mm), robust = TRUE)
  out <- list()
  for (tm in names(D$terms)) {
    tt <- topTable(efit, coef = D$terms[[tm]], number = Inf, sort.by = "none")
    ## multi-coefficient topTable returns F/P.Value; single-coefficient returns t/logFC
    st <- if ("F" %in% names(tt)) tt$F else tt$t
    ec <- eff_col(colnames(D$mm), tm)
    out[[tm]] <- data.table(gene = rownames(tt), term = tm, stat = st,
                            pval = tt$P.Value, padj = p.adjust(tt$P.Value, "BH"),
                            eff = if (is.na(ec)) NA_real_ else efit$coefficients[rownames(tt), ec])
  }
  rbindlist(out)
}

runners <- list(DESeq2 = run_deseq2, edgeR = run_edger, limma = run_limma)

## ----------------------------------------------------------------- driver ------
groups <- argv("--groups", c(hyb, "Parents"))
conds  <- argv("--conditions", c("heat", "cold"))
say("groups: %s | conditions: %s", paste(groups, collapse = ","), paste(conds, collapse = ","))

res <- list(); summ <- list()
for (cond in conds) for (g in groups) {
  idx <- if (g == "Parents") which(!design$Strain %in% hyb) else which(design$Strain == g)
  idx <- intersect(idx, which(design$Condition == cond))
  dz  <- droplevels(design[idx, , drop = FALSE])
  cw  <- cts[, idx, drop = FALSE]

  ## one shared universe per group x condition, used by every method and both models
  uni <- rowSums(cw >= MIN_CNT) >= MIN_SMP
  cw  <- cw[uni, , drop = FALSE]
  storage.mode(cw) <- "integer"

  models <- if (g == "Parents") "A" else c("A", "B")
  say("\n[%s / %-7s] n=%2d samples | %d genes pass filter | libraries %d | models %s",
      cond, g, length(idx), nrow(cw), nlevels(dz$Lib), paste(models, collapse = "+"))

  for (model in models) {
    D <- build_mm(dz, model)
    say("  model %s: %d parameters, %d residual df, terms %s", model, D$n_par, D$resid_df,
        paste(names(D$terms), collapse = "/"))
    if (D$resid_df < 2) { say("  [skip] too few residual df"); next }
    for (m in names(runners)) {
      t0 <- Sys.time()
      r <- tryCatch(runners[[m]](cw, dz, D),
                    error = function(e) { say("    [FAIL] %s/%s: %s", m, model, conditionMessage(e)); NULL })
      if (is.null(r)) next
      r[, `:=`(condition = cond, group = g, model = model, method = m)]
      res[[length(res) + 1]] <- r
      nsig <- r[, .(n = sum(padj < ALPHA, na.rm = TRUE)), by = term]
      summ[[length(summ) + 1]] <- data.table(condition = cond, group = g, model = model,
        method = m, n_gene = nrow(cw), n_par = D$n_par, resid_df = D$resid_df,
        term = nsig$term, n_sig = nsig$n,
        secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      say("    %-7s #padj<0.05: %s (%.1fs)", m,
          paste(sprintf("%s=%d", nsig$term, nsig$n), collapse = " "),
          as.numeric(difftime(Sys.time(), t0, units = "secs")))
    }
  }
}

R <- rbindlist(res, use.names = TRUE)
S <- rbindlist(summ, use.names = TRUE)
R[, Name := hog_name[match(gene, rownames(cts))]]
fwrite(R, file.path(outd, "fay_multimethod.tsv.gz"), sep = "\t")
fwrite(S, file.path(outd, "fay_multimethod_summary.tsv"), sep = "\t")
say("\nwrote fay_multimethod.tsv.gz (%d rows) and fay_multimethod_summary.tsv (%d rows)",
    nrow(R), nrow(S))
writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_fay_multimethod.txt"))
say("=== done ===")
