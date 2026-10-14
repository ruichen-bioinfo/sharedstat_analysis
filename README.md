# Analysis code and derived data for three manuscripts on the yeast heat-shock response

Rui Chen, Ludwig-Maximilians-Universität München.

- P1. Shared-component error can reverse comparative gains in hybrid allele-specific expression: reliability
  diagnostics and correction.
- P2. Msn2/4-dependent expression is sustained at 37 °C after Hsf1 anchor-away with rapamycin: a four-genotype
  reanalysis of the *Saccharomyces cerevisiae* heat-shock response.
- P3. Allele-linked divergence falls while the parental remainder rises during heat shock in six *Saccharomyces*
  hybrids.

The three share inputs (P1 and P3 use the same hybrid data; P3 uses the gene classes defined in P2), so one record
serves all three. It holds only what reproduces the reported numbers and figures.

## Contents

| path | contents |
|---|---|
| `scripts/hsf1/`, `scripts/crossspecies/`, `scripts/pipeline/` | the analysis and figure scripts behind the reported numbers and figures |
| `scripts/figures/` | the shared figure style sourced by the figure scripts |
| `fetch_inputs.sh` | downloads the public inputs of `manifests/EXTERNAL_INPUTS.tsv` and checks their sha256 |
| `run_all.sh` | checks the inputs, then reruns every script in dependency order and stops at the first failure |
| `PARAMETERS.tsv` | the fixed settings of each script (thresholds, seeds, numbers of resamples, grids): value, meaning and the output it affects |
| `data/02_processed/` | the derived tables used to build the reported numbers and figures and read by downstream scripts; a full run also writes further diagnostic tables |
| `data/01_metadata/` | sample sheet |
| `doc/04_figures/final/` | the 25 figures as submitted (P2 main figures and P3 figures also as numbered copies), and `figure_map_P1.tsv`, `figure_map_P2.tsv` (main figures), `figure_map_P2_supplement.tsv` and `figure_map_P3.tsv` from figure numbers to files |
| `manifests/P1_NUMBER_PROVENANCE.tsv`, `P2_…`, `P3_…` | for each number in a journal manuscript, the tables and columns in which that value appears |
| `manifests/VALUES_NOT_IN_TABLES.tsv` | the two reported values that no table stores, and where they are printed |
| `manifests/FIGURES.tsv` | each figure and the script that draws it |
| `manifests/TABLES.tsv` | the derived tables with size, sha256, rows and columns (per-comparison differential-expression tables excluded) |
| `manifests/EXTERNAL_INPUTS.tsv` | the public input files the scripts read and that are not in this record, with source URL, fetch date, size and sha256 |
| `ENVIRONMENT/` | R `sessionInfo()`, installed R package versions and the conda package list |
| `proofs/lean/` | Lean 4 (v4.14.0, Mathlib) checks of the algebra in P1 Supplementary Note S2 |
| `REPRODUCTION.md` | the result of rerunning everything |
| `SHA256SUMS.txt` | checksums of every file |

The number tables are candidate lists: a value can also appear by coincidence in an unrelated column. Sections that report no results (data availability, acknowledgements, references, disclosures) are not included.

In the P3 tables, columns named `trans` or `trn` hold the parental remainder, the parental difference minus *cis*. The manuscript calls it the remainder because it contains more than regulatory *trans* effects.

## Public data used

No raw sequencing reads (FASTQ) are included. The record does include derived count and quantification files:
`data/02_processed/count_matrix_aligned.tsv` (GSE179258) and the salmon quantifications in `data/02_processed/hov2020/quant/`
(PRJNA576452).

| dataset | accession | used in |
|---|---|---|
| Hsf1 anchor-away × *msn2Δ msn4Δ* heat time course (Mühlhofer et al. 2024) | GEO GSE179258 (counts; and the authors' DESeq2 table, used to check our fold changes) | P2; gene classes for P1 and P3 |
| Six *Saccharomyces* hybrids and their parents, heat time course (Fay et al. 2023) | osf.io/7kunf; BioProject PRJNA909640 | P1, P3 |
| *S. cerevisiae* × *S. uvarum* cold shift (Hovhannisyan et al. 2020) | BioProject PRJNA576452 | P1, P3 |
| Two-replicate yeast data (Keyport Kik et al. 2024) | GEO GSE234499 | P1, P3 |
| Nascent transcription after Hsf1 anchor-away, NAC-seq (Pincus et al. 2018) | PMID 30332327, PMC6340206 (supplementary data) | P3 |
| SEQC mixtures | GEO GSE47774 | P1 |
| HEK293 cells (Ronnebaum and Schisler 2025) | GEO GSE310700 | P1 |

The GSE179258 counts enter as `data/02_processed/count_matrix_aligned.tsv`, assembled from the 60 per-sample count
files deposited in GEO. All other public inputs (the Fay et al. counts and genomes, the GSE234499, GSE47774 and
GSE310700 count files, Gene Ontology and SGD annotation, KEGG, YEASTRACT, the S. kudriavzevii ZP591 assembly, the
ERCC concentrations, the Pincus et al. NAC-seq table and the Ensembl annotation) are listed in
`manifests/EXTERNAL_INPUTS.tsv` with URL, fetch date and sha256. `./fetch_inputs.sh` downloads them to the listed
paths and checks each file. Two limits: the YEASTRACT+ 2019 flat file has no stable public URL and must be obtained from
YEASTRACT+ and placed at its path by hand; KEGG REST and the GEO text pages are live services, so a later download can
differ from the listed checksum, which the script reports.

## Running

R 4.4.3 with data.table, limma, edgeR and DESeq2 (all packages and versions in `ENVIRONMENT/`), and DIAMOND for the
*S. kudriavzevii* ortholog map (script 57; set `DIAMOND_BIN` to its path if it is not on the `PATH`). The scripts load the R package
sharedstat 0.1.3 from source: download it (doi:10.5281/zenodo.23070215) and unpack it into `pkg/sharedstat`. Then

    ./fetch_inputs.sh
    ./run_all.sh

`run_all.sh` first checks that sharedstat and every listed input are in place with the listed checksum and runs
nothing otherwise. Each step writes `logs/NNN_<script>.log`; the run stops at the first failing step, names it, writes
it to `logs/FAILED` and exits with a non-zero status, because later steps read its outputs.

Each script resolves its paths from `PBL_ROOT`, which `run_all.sh` sets to this directory. A single script runs the
same way, for example `PBL_ROOT=$PWD Rscript scripts/hsf1/37_p1_a1_provenance_lock.R`, which re-derives the P1
headline numbers from the tables and stops on a mismatch.

The Lean file: install elan, then in `proofs/lean/` run `lake exe cache get`, `lake build` and
`lake env lean check_axioms.lean`. Each of the 14 theorems depends only on `propext`, `Classical.choice` and
`Quot.sound`.

## Licence

Code: MIT (`LICENSE`). Derived tables: the author's contribution is released under CC0 1.0. Tables derived from the
public datasets above remain subject to the terms of those sources; please cite the original data producers. The
external inputs are not redistributed here.
