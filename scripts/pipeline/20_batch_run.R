# 20_batch_run.R — batch DE run over the comparison registry, fixed before the run
#
# Design choices implemented here:
#   * shared gene universe across all methods and comparisons (if each method applied
#     its own filter, set sizes would not be comparable)
#   * dual reporting: padj on the shared universe (primary) + padj under each
#     method's native filtering (padj_native) so the effect of the choice is visible
#   * per-comparison two-group fit; the full interaction design is reserved for the
#     interaction analysis, not for these pairwise calls
#   * limma lfcSE taken from the moderated eBayes output rather than back-derived
#   * test/reference arms are built from the registry columns, then cross-checked
#     against comp_id, so the direction convention is checked rather than assumed
#
# Usage: Rscript scripts/pipeline/20_batch_run.R [--only=ID1,ID2] [--limit=N] [--force]
#   --only   run just these comp_ids (smoke-test the script before the full run)
#   --limit  run only the first N registered comparisons
#   --force  recompute comparisons whose output files already exist

suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(edgeR); library(limma)
})

root <- Sys.getenv("PBL_ROOT", unset = NA_character_)
if (is.na(root)) root <- if (dir.exists(".")) "." else normalizePath(getwd())
proc <- file.path(root, "data/02_processed")
repd <- file.path(root, "doc/05_reports")
logd <- file.path(root, "doc/05_reports/logs")
dir.create(logd, showWarnings = FALSE, recursive = TRUE)

args    <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) NA_character_ else sub(paste0("^", flag, "="), "", hit[1])
}
only  <- arg_val("--only")
limit <- suppressWarnings(as.integer(arg_val("--limit")))
force <- "--force" %in% args

PREFILTER_X   <- 3      # rowSums(count >= 10) >= x, frozen on the full 60-sample matrix
ALPHA         <- 0.05
SEED          <- 20260905
## Bump whenever the method code changes: it is stamped into every output file and
## cached results with a different method version are recomputed.
PIPELINE_VERSION <- "20.1"
METHOD_SCHEMA <- c("gene", "baseMean", "log2FC", "lfcSE", "stat", "pval", "padj", "direction")
set.seed(SEED)

logfile <- file.path(logd, sprintf("batch_run_%s.log", format(Sys.Date())))
say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = logfile, append = TRUE)
}
cat("", file = logfile)   # truncate previous log for this date
say("=== 20_batch_run.R start %s ===", format(Sys.time()))
say("root=%s | R=%s", root, R.version.string)

## ---------------------------------------------------------------- inputs ----
cntdf <- fread(file.path(proc, "count_matrix_aligned.tsv"), check.names = FALSE)
cnt <- as.matrix(cntdf[, -1]); rownames(cnt) <- cntdf$gene_id
storage.mode(cnt) <- "integer"
map <- read.delim(file.path(root, "data/01_metadata/sample_mapping.tsv"), stringsAsFactors = FALSE)
reg <- fread(file.path(proc, "comparisons.csv"))
stopifnot(identical(colnames(cnt), map$sample), nrow(reg) > 0)
say("counts %d x %d | registry %d comparisons", nrow(cnt), ncol(cnt), nrow(reg))

## EmpiReS is optional; without it the consensus is downgraded to silver = >=2/3 majority
have_empires <- requireNamespace("EmpiReS", quietly = TRUE)
methods_used <- c("DESeq2", "edgeR", "limma")
say("EmpiReS available: %s -> %d methods (%s)", have_empires, length(methods_used),
    if (have_empires) "four-method gold/silver" else "EmpiReS not available: silver = >=2/3 majority")
for (m in methods_used) dir.create(file.path(proc, m), showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------- shared universe ----
## One gene set for every method and every comparison. Independent filtering and
## edgeR's filterByExpr are switched off below so that padj is computed over this
## same set; each method's native filter is still recorded per comparison.
universe <- rownames(cnt)[rowSums(cnt >= 10) >= PREFILTER_X]
say("shared universe (rowSums(count>=10)>=%d on all 60 samples): %d genes of %d",
    PREFILTER_X, length(universe), nrow(cnt))
cnt_u <- cnt[universe, , drop = FALSE]

## --------------------------------------------------------- arm resolution ----
## Every family reduces to the same shape: a test arm and a reference arm defined
## by (knockout, temperature, time). Resolving from the registry columns and then
## checking comp_id is what checks the direction convention.
arm_samples <- function(g, t, s) map$sample[map$knockout == g & map$temperature == t & map$time == s]

resolve_arms <- function(r) {
  test <- arm_samples(r$genotype,     r$temperature,     r$time)
  ref  <- arm_samples(r$genotype_ref, r$temperature_ref, r$time_ref)
  list(test = test, ref = ref)
}

## comp_id encodes the same contrast; disagreement means the registry and the id
## drifted apart and the run must stop rather than silently mislabel a direction.
check_comp_id <- function(r) {
  id <- r$comp_id
  expect <- switch(r$family,
    time_series = sprintf("%s_%d_%d_vs_%d_%d", r$genotype, r$temperature, r$time,
                          r$temperature_ref, r$time_ref),
    temperature = sprintf("%s_%d_%d_vs_%d_%d", r$genotype, r$temperature, r$time,
                          r$temperature_ref, r$time_ref),
    genotype    = sprintf("%s_vs_%s_at_%d_%d", r$genotype, r$genotype_ref,
                          r$temperature, r$time),
    NA_character_)
  if (is.na(expect) || !identical(id, expect))
    stop(sprintf("comp_id/registry mismatch: id='%s' but columns imply '%s'", id, expect))
  invisible(TRUE)
}

dir_from <- function(lfc, padj) {
  fifelse(is.na(padj) | padj >= ALPHA, "ns", fifelse(lfc > 0, "up", "down"))
}

## ------------------------------------------------------------- methods ------
run_deseq2 <- function(cw, arm) {
  cd <- data.frame(arm = arm, row.names = colnames(cw))
  dds <- DESeqDataSetFromMatrix(cw, colData = cd, design = ~ arm)
  dds <- DESeq(dds, quiet = TRUE)
  cf  <- resultsNames(dds)[2]                      # arm_test_vs_ref
  ## primary: shared universe -> no independent filtering, no Cook's p-value removal
  res <- results(dds, name = cf, alpha = ALPHA,
                 independentFiltering = FALSE, cooksCutoff = FALSE)
  shr <- tryCatch(lfcShrink(dds, coef = cf, type = "apeglm", res = res, quiet = TRUE),
                  error = function(e) { say("    [warn] apeglm shrink failed: %s", conditionMessage(e)); res })
  ## native: DESeq2 defaults, for the dual report
  nat <- results(dds, name = cf, alpha = ALPHA)
  out <- data.table(gene = rownames(shr), baseMean = shr$baseMean,
                    log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE,
                    stat = res$stat, pval = res$pvalue, padj = res$padj)
  out[, direction := dir_from(log2FC, padj)]
  out[, padj_native := nat$padj[match(gene, rownames(nat))]]
  out[, filter_native := !is.na(padj_native)]
  out[]
}

run_edger <- function(cw, arm) {
  des <- model.matrix(~ arm)
  y <- DGEList(counts = cw, group = arm)
  native_keep <- filterByExpr(y, group = arm)      # recorded, not applied
  y <- calcNormFactors(y)
  y <- estimateDisp(y, des)
  fit <- glmQLFit(y, des)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- qlf$table
  out <- data.table(gene = rownames(tt), baseMean = 2^(y$AveLogCPM),
                    log2FC = tt$logFC,
                    lfcSE = NA_real_,              # QL framework exposes no logFC SE
                    stat = tt$F, pval = tt$PValue,
                    padj = p.adjust(tt$PValue, "BH"))
  out[, direction := dir_from(log2FC, padj)]
  out[, filter_native := native_keep[match(gene, names(native_keep))]]
  out[, padj_native := NA_real_]
  out[filter_native == TRUE, padj_native := p.adjust(pval, "BH")]
  out[]
}

run_limma <- function(cw, arm) {
  des <- model.matrix(~ arm)
  native_keep <- filterByExpr(cw, group = arm)
  ## TMM before voom: without it voom normalises on library size alone, and the
  ## heat-shock libraries are compositionally skewed (a few highly induced genes
  ## take a large share of reads), which biases log-fold-changes directionally
  y <- calcNormFactors(DGEList(counts = cw, group = arm))
  v <- voom(y, des, plot = FALSE)
  efit <- eBayes(lmFit(v, des), robust = TRUE)
  tt <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
  ## moderated standard error straight from eBayes, not logFC/t
  se <- sqrt(efit$s2.post) * efit$stdev.unscaled[, 2]
  out <- data.table(gene = rownames(tt), baseMean = 2^tt$AveExpr, log2FC = tt$logFC,
                    lfcSE = se[rownames(tt)], stat = tt$t,
                    pval = tt$P.Value, padj = tt$adj.P.Val)
  out[, direction := dir_from(log2FC, padj)]
  out[, filter_native := native_keep[match(gene, names(native_keep))]]
  out[, padj_native := NA_real_]
  out[filter_native == TRUE, padj_native := p.adjust(pval, "BH")]
  out[]
}

runners <- list(DESeq2 = run_deseq2, edgeR = run_edger, limma = run_limma)

## --------------------------------------------------------------- driver -----
todo <- reg
if (!is.na(only)) todo <- todo[comp_id %in% trimws(strsplit(only, ",")[[1]])]
if (!is.na(limit)) todo <- head(todo, limit)
say("running %d of %d registered comparisons%s", nrow(todo), nrow(reg),
    if (force) " (--force: recomputing existing outputs)" else "")
if (nrow(todo) == 0) stop("nothing to run: --only/--limit selected no registered comparison")

summ <- list(); failures <- character()

for (i in seq_len(nrow(todo))) {
  r <- todo[i]
  say("[%d/%d] %s (%s)", i, nrow(todo), r$comp_id, r$family)
  ok <- tryCatch({
    check_comp_id(r)
    arms <- resolve_arms(r)
    if (length(arms$test) != r$n_rep_per_arm || length(arms$ref) != r$n_rep_per_arm)
      stop(sprintf("arm size mismatch: test=%d ref=%d, registry expects %d each",
                   length(arms$test), length(arms$ref), r$n_rep_per_arm))
    if (length(intersect(arms$test, arms$ref)) > 0)
      stop("test and reference arms share samples")

    smp <- c(arms$ref, arms$test)
    arm <- factor(c(rep("ref", length(arms$ref)), rep("test", length(arms$test))),
                  levels = c("ref", "test"))
    cw <- cnt_u[, smp, drop = FALSE]

    for (m in methods_used) {
      outf <- file.path(proc, m, paste0(r$comp_id, ".tsv"))
      ## a cached file is valid only under the current
      ## universe: outputs with a different gene set
      ## or without the dual-report columns are incompatible
      cached <- FALSE
      if (file.exists(outf) && !force) {
        res <- fread(outf)
        cached <- all(c(METHOD_SCHEMA, "padj_native", "filter_native", "pipeline_version") %in% names(res)) &&
                  nrow(res) == length(universe) &&
                  identical(as.character(res$pipeline_version[1]), PIPELINE_VERSION)
        if (!cached) say("    %-7s stale output (different config) -> recomputing", m)
      }
      if (cached) {
        say("    %-7s cached", m)
      } else {
        res <- runners[[m]](cw, arm)
        stopifnot(all(METHOD_SCHEMA %in% names(res)))
        res[, pipeline_version := PIPELINE_VERSION]
        setcolorder(res, METHOD_SCHEMA)
        fwrite(res, outf, sep = "\t")
      }
      nde <- sum(res$padj < ALPHA, na.rm = TRUE)
      summ[[length(summ) + 1]] <- data.table(
        comp_id = r$comp_id, family = r$family, method = m,
        n_tested = sum(!is.na(res$padj)),
        n_universe = length(universe),
        n_native_filter = sum(res$filter_native, na.rm = TRUE),
        n_de = nde,
        n_up = sum(res$direction == "up", na.rm = TRUE),
        n_down = sum(res$direction == "down", na.rm = TRUE),
        n_de_native = sum(res$padj_native < ALPHA, na.rm = TRUE),
        median_abs_lfc_de = if (nde > 0) median(abs(res$log2FC[which(res$padj < ALPHA)])) else NA_real_,
        test_arm = paste(arms$test, collapse = ";"),
        ref_arm  = paste(arms$ref,  collapse = ";"))
      say("    %-7s #DE=%d (up %d / down %d), native-filter #DE=%d",
          m, nde, summ[[length(summ)]]$n_up, summ[[length(summ)]]$n_down,
          summ[[length(summ)]]$n_de_native)
    }
    TRUE
  }, error = function(e) { say("    [FAIL] %s", conditionMessage(e)); FALSE })
  if (!isTRUE(ok)) failures <- c(failures, r$comp_id)
}

S <- rbindlist(summ)
fwrite(S, file.path(proc, "batch_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(logd, sprintf("sessionInfo_%s.txt", format(Sys.Date()))))

## ---------------------------------------------------------------- report ----
wide <- dcast(S, comp_id + family ~ method, value.var = "n_de")
L <- c(sprintf("# Batch run report (%s, scripts/pipeline/20_batch_run.R)", format(Sys.Date())), "",
  sprintf("- registry %d comparisons, run %d, failed %d", nrow(reg), nrow(todo), length(failures)),
  sprintf("- methods available: %d (%s)%s", length(methods_used), paste(methods_used, collapse = "/"),
          if (have_empires) "" else " -- EmpiReS not available: silver = >=2/3 majority"),
  sprintf("- **common universe**: `rowSums(count>=10)>=%d` (all 60 samples) -> %d genes, shared by every method and comparison",
          PREFILTER_X, length(universe)),
  "- two padj columns: the primary padj on the common universe (DESeq2 without independent filtering or Cook's cutoff, edgeR/limma without filterByExpr); padj under each method's own filtering is kept as `padj_native`",
  sprintf("- random seed %d; runtime details and per-comparison logs are recorded in `%s`",
          SEED, basename(logfile)),
  "- direction checked per comparison: test/reference are read from the registry and checked against comp_id; a mismatch stops the run",
  "- edgeR `lfcSE` is NA: the quasi-likelihood framework gives no logFC standard error (limma: eBayes-moderated SE; DESeq2: apeglm-shrunken SE)", "")

if (length(failures) > 0)
  L <- c(L, "## Failed comparisons", paste0("- `", failures, "`"), "")

L <- c(L, "## DE genes per comparison (padj < 0.05, common universe)", "",
  paste0("| comp_id | family | ", paste(methods_used, collapse = " | "), " |"),
  paste0("|---|---|", paste(rep("---", length(methods_used)), collapse = "|"), "|"))
for (i in seq_len(nrow(wide)))
  L <- c(L, sprintf("| %s | %s | %s |", wide$comp_id[i], wide$family[i],
                    paste(unlist(wide[i, ..methods_used]), collapse = " | ")))

fam <- S[, .(med = as.integer(median(n_de))), by = .(family, method)]
famw <- dcast(fam, family ~ method, value.var = "med")
L <- c(L, "", "## By family (median DE genes per method)", "",
  paste0("| family | ", paste(methods_used, collapse = " | "), " |"),
  paste0("|---|", paste(rep("---", length(methods_used)), collapse = "|"), "|"))
for (i in seq_len(nrow(famw)))
  L <- c(L, sprintf("| %s | %s |", famw$family[i], paste(unlist(famw[i, ..methods_used]), collapse = " | ")))
L <- c(L, "", "")

writeLines(L, file.path(repd, sprintf("batch_run_%s.md", format(Sys.Date()))))

say("=== done: %d comparisons x %d methods, %d failures ===", nrow(todo), length(methods_used), length(failures))
say("summary -> data/02_processed/batch_summary.tsv | dated run report: batch_run_%s.md", format(Sys.Date()))
if (length(failures) > 0) quit(status = 1)
