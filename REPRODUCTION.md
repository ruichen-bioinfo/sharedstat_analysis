# Reproduction

On 2026-10-12 this record was copied into an empty directory, sharedstat 0.1.3 was unpacked into `pkg/sharedstat`, and
`./fetch_inputs.sh` placed the public inputs listed in `manifests/EXTERNAL_INPUTS.tsv` and checked each against its
sha256. The YEASTRACT+ file was placed by hand, because it has no stable public URL. Every input matched its listed
checksum except the ENA run report, whose rows came back in a different order with the same content. `./run_all.sh`
then ran all 116 steps in one pass, from step 1, without a failure, with the environment recorded in `ENVIRONMENT/`.
The salmon quantifications of the cold-shift data and the GSE179258 count matrix shipped here were used as inputs.

The rerun tables were compared with the shipped ones after numbers were rounded to six significant figures.

| result | tables |
|---|---|
| identical | 236 |
| different | 5 |

All 25 figures used in the manuscripts were regenerated (`doc/04_figures/final/figure_map_P1.tsv`, `_P2.tsv`,
`_P2_supplement.tsv`, `_P3.tsv`), and P1 Table 1 (`doc/04_figures/final/TABLE1_per_hybrid.md`) was identical. Figures
were checked for existence only, because a PDF embeds its creation time.

The 5 differences:

- `crossspecies/fay_multimethod.tsv.gz`: 20 of 995,601 effect values differ by more than 10⁻⁶ relative, the largest
  by 30%; test statistics and p-values differ by at most 0.9% and 0.06%. The tables computed from it are identical.
- Labels or order only: `hsf1/a5_cis_dod_37_30.tsv.gz` (row order), `hsf1/p1_intergene_correlation.tsv` (genes
  without a class are labelled `unassigned` instead of left empty) and `hsf1/s4_module_stratum_split3_beta.tsv`
  (column order). No value differs.
- `comparisons.csv`: the shipped file lists 12 time-change comparisons that the registry script does not write; no
  script in this record uses their result tables.

After that run, comments and printed messages were edited in eleven scripts, and the printed result summary of
`25_b4_reference_class_simulation.R` was rewritten to read on the responding-class quantity; no computation changed.
In the same directory the edited scripts before step 95 (25, 32, 27, 37, 36, 54, 41, 55, 64) were rerun, and
`run_all.sh` was run again from step 95 to step 116, which includes the other two (15, 19). Every table was identical
to the one written by the full run.
A later edit of printed messages in scripts 25, 54 and 55 was checked the same way: those scripts and script 41
were rerun and every table was unchanged.

`scripts/hsf1/37_p1_a1_provenance_lock.R`, which re-derives the P1 headline numbers and stops on a mismatch, passed
(29 of 29); `scripts/hsf1/111_p3_composition_check.R` reproduced its result; and
`scripts/hsf1/29_geo_published_lfc_check.R` reproduced the agreement with the published GSE179258 fold changes
(`hsf1/geo_published_lfc_agreement.tsv`, identical).

Two steps were added to `run_all.sh` after that run: `scripts/hsf1/112_p3_constant_share_null.R` (the constant-share
null of P3 Supplementary Note N2, section 17) and `scripts/hsf1/85_p3_parental_axis_tests.R`, whose two tables had been
written by an earlier run only. Run on their own from the same inputs, both reproduced their tables exactly.

## Changes in version 1.0.1

The full rerun above was made for version 1.0.0 and was not repeated for version 1.0.1. Version 1.0.1 adds three
steps to `run_all.sh`: `scripts/crossspecies/02b_p1order_refit.R` after step 1, and two runs of
`scripts/crossspecies/11_gain_iv.R` with `P1_ORDER=name` (without and with `GAIN_IV_SELECT=padj`) after its two
alphabetical runs, so `run_all.sh` now has 121 steps. The step numbers above refer to `run_all.sh` as it was at
the time of that run.

For version 1.0.1 these scripts were run in a copy of the analysis tree, on the derived tables of version 1.0.0:
the refit; script 11 in its two alphabetical runs, which reproduced their tables byte for byte, and in its two
name-order runs; then scripts 94, 12, 102, 102c, 102d, 92, 41, 110 and 37. The refit's checks passed: with the original
allele levels it reproduced the stored fits, and its effects equal the stored effects re-expressed with
*S. paradoxus* as reference in the heat, model A fits that P1 uses. In the other condition and model cells, which P1
does not use, that comparison is reported but does not stop the script; there four model-B cells differ by up to
1.8 × 10⁻⁴. `scripts/hsf1/37_p1_a1_provenance_lock.R` passed (29 of 29). The tables and the figure written by these
runs are the ones shipped here.

The record itself was then checked on 2026-10-17. In a copy of this record, with sharedstat 0.1.3 unpacked into
`pkg/sharedstat` and the public inputs placed at their listed paths (taken from the analysis tree, not downloaded
again), the same scripts were run in the same order: the refit, the two name-order runs of script 11, then scripts
94, 12, 102, 102c, 102d, 92, 41, 110 and 37. Every shipped table they write was identical to the shipped file
(sha256): `fay_multimethod_p1order.tsv.gz`, `gain_iv_p1order.tsv`, `gain_iv_padj_filter_p1order.tsv`,
`gain_components_p1order.tsv`, the six P1 tables listed in `README.md` and P1 Table 1. The PNG of Figure 5 was
identical; its PDF differs only in the creation date, and both PDFs render to identical pages. No other file of the
record changed. The refit's checks passed and script 37 passed (29 of 29).
