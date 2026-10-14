#!/usr/bin/env bash
# run_all.sh -- rerun every analysis and figure script of this record, in dependency order.
# Needs: R 4.4 with the packages in ENVIRONMENT/; the sharedstat 0.1.3 source unpacked in pkg/sharedstat
# (doi:10.5281/zenodo.23070215); the public inputs of manifests/EXTERNAL_INPUTS.tsv at the listed paths (./fetch_inputs.sh).
# Outputs are written to data/02_processed and doc/04_figures/final, replacing the files shipped here.
# Several scripts take minutes; 13, 33, 105 and 105b take about 20 to 40 minutes each.
set -u -o pipefail
cd "$(dirname "$0")"; export PBL_ROOT="$PWD"
LOG="$PWD/logs"; mkdir -p "$LOG" doc/05_reports/logs # run logs and the scripts' report folders
## preflight: R, the sharedstat source, and every public input at its listed path with the listed sha256
command -v Rscript > /dev/null || { echo "Rscript not found"; exit 1; }
[ -d pkg/sharedstat/R ] || { echo "unpack sharedstat 0.1.3 into pkg/sharedstat first"; exit 1; }
missing=0
while IFS=$'\t' read -r file url fetched bytes sha note; do
  [ "$file" = "file" ] && continue
  if [ ! -f "$file" ]; then echo "  missing input: $file"; missing=1
  elif [ "$(sha256sum "$file" | cut -c1-64)" != "$sha" ]; then echo "  input differs from the manifest: $file"; missing=1; fi
done < <(sed -e 's/\t\t/\t-\t/g; s/\t\t/\t-\t/g' manifests/EXTERNAL_INPUTS.tsv)
[ $missing -eq 0 ] || { echo "Inputs missing or changed; run ./fetch_inputs.sh (see README). Nothing was run."; exit 1; }
## each step logs to logs/NN_<script>.log; the first failure stops the run, because later steps read its outputs
step=0
run() {
  step=$((step + 1)); local name; name=$(printf "%03d_%s" "$step" "$(basename "$1" .R)")
  echo "== [$step] $*"
  if ! Rscript "$@" > "$LOG/$name.log" 2>&1; then
    echo "   FAILED at step $step: $*   (log: logs/$name.log)"; tail -5 "$LOG/$name.log"
    echo "$step $*" > "$LOG/FAILED"; exit 1
  fi
}
rm -f "$LOG/FAILED"
run scripts/crossspecies/02_fay_multimethod.R
run scripts/pipeline/03_comparisons_registry.R
run scripts/hsf1/07_hsf1_locus_extract.R
run scripts/hsf1/101_p1_reliability_vs_iv_sim.R
run scripts/hsf1/109_p1_arm_error_covariance.R
run scripts/hsf1/12_upstream_sensor_conservation.R
run scripts/hsf1/24_fig_a2_ladder.R
run scripts/hsf1/25_b4_reference_class_simulation.R
run scripts/hsf1/31_p1_design_grid.R
run scripts/hsf1/35_p1_replicate_axis.R
run scripts/hsf1/40_fig_p1_architecture.R
run scripts/hsf1/43_fig_p1_robustness.R
run scripts/hsf1/50_c5_build_index.R
run scripts/hsf1/55_p1_fourth_dataset_kk2024.R
run scripts/hsf1/83_p1_iv_coverage_correlated.R 0.017
run scripts/hsf1/83_p1_iv_coverage_correlated.R 0.05
run scripts/hsf1/89_p1_seqc_known_truth.R BGI
run scripts/hsf1/89_p1_seqc_known_truth.R MAY
run scripts/hsf1/98_fig_p1_schematic.R
run scripts/hsf1/99_fig_p2_schematic.R
run scripts/hsf1/08_hsf1_domain_divergence.R
run scripts/crossspecies/102d_p1_joint_cluster_bootstrap.R
run scripts/hsf1/103b_p2_torc1_panel.R
run scripts/hsf1/104b_p2_class_split.R
run scripts/crossspecies/11_gain_iv.R
GAIN_IV_SELECT=padj run scripts/crossspecies/11_gain_iv.R
run scripts/pipeline/20_batch_run.R
run scripts/hsf1/29_geo_published_lfc_check.R
run scripts/hsf1/32_fig_p1_grid.R
run scripts/hsf1/39_fig_p1_attenuation.R
run scripts/hsf1/75_p1_iv_coverage.R
run scripts/hsf1/90_fig_p1_seqc.R
run scripts/hsf1/92_p1_iv_selection_readings.R
run scripts/hsf1/01_epistasis_decomposition.R
run scripts/crossspecies/06_module_bridge.R
run scripts/hsf1/110_fig_p1_hybrid_gains.R
run scripts/hsf1/02_core_function.R
run scripts/hsf1/03_fay_parental_dimensions.R
run scripts/hsf1/05_cis_trans_decomposition.R
run scripts/crossspecies/08_thermophile_allele_polarity.R
run scripts/hsf1/104_p2_class_split.R
run scripts/hsf1/107_p2_per_library_score.R
run scripts/crossspecies/14_mechanism_tests.R
run scripts/hsf1/17_p2_interaction_se_and_four_conditions.R
run scripts/hsf1/20_p2_offset_origin.R
run scripts/hsf1/27_pincus_external_core_control.R
run scripts/hsf1/91_p1_reference_arm_reliability.R
run scripts/hsf1/93_p1_reference_arm_reliability_ci.R
run scripts/hsf1/94_p1_anchor_diagnosis.R
run scripts/hsf1/95_p2_shutoff_timecourse.R
run scripts/hsf1/96_p2_time_contrasts.R
run scripts/hsf1/04_divergence_by_dependency.R
run scripts/hsf1/09_promoter_hse_architecture.R
run scripts/crossspecies/102_p1_table1_uncertainty.R
run scripts/hsf1/103_p2_torc1_readout.R
run scripts/hsf1/108_p3_time_invariant_null.R
run scripts/hsf1/112_p3_constant_share_null.R
run scripts/hsf1/10_promoter_hse_architecture_v2.R
run scripts/hsf1/11_chance_calibration_highres.R
run scripts/crossspecies/12_iv_breakdown.R
run scripts/hsf1/21_p2_exclusive_classes.R
run scripts/hsf1/22_fig_f2d.R
run scripts/hsf1/33_core_definition_robustness.R
run scripts/hsf1/37_p1_a1_provenance_lock.R
run scripts/hsf1/53_p3_r3_per_hybrid.R
run scripts/hsf1/60_p2_fig3_two_backgrounds.R
run scripts/hsf1/87b_p3_cisshare_adjusted.R
run scripts/hsf1/87c_p3_composition_adjusted.R
run scripts/hsf1/97_fig_p2_time_course.R
run scripts/hsf1/111_p3_composition_check.R
for sh in 1 2 3; do SHARD=$sh N_SUB_SHARD=14 run scripts/hsf1/13_resolution_shard.R; done
run scripts/hsf1/14_background_amplification_diagnostic.R
run scripts/hsf1/34_core_definition_strict.R
run scripts/hsf1/36_p1_measured_reliability.R
run scripts/hsf1/54_p1_hovhannisyan_dataset.R
run scripts/hsf1/41_table1_per_hybrid.R
run scripts/hsf1/46_p2_sequence_orthogonal_validation.R
run scripts/hsf1/52_c5_decomposition.R
run scripts/hsf1/53b_p3_r3_per_hybrid_nooverlap.R
run scripts/hsf1/79_p2_replicate_jackknife.R
run scripts/hsf1/11_b3v2_resolution_bound.R
run scripts/hsf1/45_core_definition_at_headline_beta.R
run scripts/hsf1/49_p2x_full_scan.R
run scripts/hsf1/64_b5_centred_reliability_with_ci.R
run scripts/hsf1/23_l6_bound_conversion.R
run scripts/hsf1/38_fig_p1_theory_and_reliability.R
run scripts/hsf1/74_p2_joint_candidates.R
run scripts/hsf1/82_p1_h_correlation.R
run scripts/hsf1/87a_p3_core_adjusted.R
run scripts/crossspecies/102c_p1_iv_widened.R
run scripts/hsf1/100_p2_temperature_difference.R
run scripts/hsf1/105_p3_null_calibration.R
run scripts/hsf1/105b_p3_power_injection.R
run scripts/hsf1/106_fig_supplement.R
run scripts/hsf1/15_module_regulatory_profile.R
run scripts/hsf1/57_skud_ortholog_map.R
run scripts/hsf1/16_s4_stem_collapse.R
run scripts/hsf1/18_p2_amplitude_and_cistrans_identity.R
run scripts/hsf1/26_p2_amplitude_bypair_readback.R
run scripts/hsf1/72_p2_i25.R
run scripts/hsf1/76_p2_i25_resolvability.R
run scripts/hsf1/77_p2_sensitivity.R
run scripts/hsf1/78_p2_detectability.R
run scripts/hsf1/80_p2_adjusted_uncertainty.R
run scripts/hsf1/84_p2_rho2_sensitivity.R
run scripts/hsf1/19_figures.R
run scripts/hsf1/28_cold_arm_specificity.R
run scripts/hsf1/30_amplitude_matched_module_test.R
run scripts/hsf1/44_p3_module_stratum_crosstab.R
run scripts/hsf1/47_cold_arm_regulon_split.R
run scripts/hsf1/48_flag_degenerate_cells.R
run scripts/hsf1/56_residual_polarity.R
run scripts/hsf1/58_species_matched_shock.R
run scripts/hsf1/59_name_residual_genes.R
run scripts/hsf1/85_p3_parental_axis_tests.R
run scripts/hsf1/86_p3_parental_axis_scope.R
echo "All $step steps finished; logs in logs/."
