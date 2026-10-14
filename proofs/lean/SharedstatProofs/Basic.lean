/-
  SharedstatProofs — machine-checked derivations for the shared-component estimator.

  WHAT IS AND IS NOT PROVED HERE. Covariance is modelled as a symmetric bilinear form `B` on a real
  vector space, with `var x = B x x`. Every orthogonality the paper assumes is an explicit hypothesis.
  What the file establishes is therefore the ALGEBRA: given that second-order structure, each stated
  identity follows. It does not establish the probabilistic model itself — that R and u are independent,
  that the errors are additive — because those are modelling assumptions about a dataset, not theorems.
  The one result in the paper that genuinely needs a distribution, E|e| = 2σ/√π for Gaussian errors, is
  deliberately absent here and is checked numerically instead.

  The expansion of `B e Rhat` is the step where an error could hide, so it is proved rather than assumed.
-/
import Mathlib

namespace Sharedstat

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-- Second-order structure: a symmetric bilinear form standing for covariance. -/
structure Cov (V : Type*) [AddCommGroup V] [Module ℝ V] where
  B : V → V → ℝ
  add_left : ∀ x y z, B (x + y) z = B x z + B y z
  smul_left : ∀ (c : ℝ) (x y), B (c • x) y = c * B x y
  symm : ∀ x y, B x y = B y x

namespace Cov

variable (C : Cov V)

lemma add_right (x y z : V) : C.B x (y + z) = C.B x y + C.B x z := by
  rw [C.symm x (y + z), C.add_left, C.symm y x, C.symm z x]

lemma smul_right (c : ℝ) (x y : V) : C.B x (c • y) = c * C.B x y := by
  rw [C.symm x (c • y), C.smul_left, C.symm y x]

lemma neg_left (x y : V) : C.B (-x) y = -C.B x y := by
  have h : (-x : V) = (-1 : ℝ) • x := by simp
  rw [h, C.smul_left]; ring

lemma sub_left (x y z : V) : C.B (x - y) z = C.B x z - C.B y z := by
  rw [sub_eq_add_neg, C.add_left, C.neg_left]; ring

lemma sub_right (x y z : V) : C.B x (y - z) = C.B x y - C.B x z := by
  rw [C.symm x (y - z), C.sub_left, C.symm y x, C.symm z x]

/-- `var x = cov x x`. -/
def var (x : V) : ℝ := C.B x x

end Cov

/-- The architecture: one measured reference reused as both a component of the contrast and the axis the
contrast is regressed on. `R` is the reference signal, `u` its measurement error, `v` the response arm's
own error, `γ` the gain the analyst is trying to estimate. -/
structure Design (V : Type*) [AddCommGroup V] [Module ℝ V] where
  C : Cov V
  R : V
  u : V
  v : V
  γ : ℝ
  R_u : C.B R u = 0
  R_v : C.B R v = 0
  u_v : C.B u v = 0

namespace Design

variable (D : Design V)

/-- The reference as measured. -/
def Rhat : V := D.R + D.u
/-- The response arm as measured. -/
def A : V := D.γ • D.R + D.v
/-- The contrast the analyst forms. -/
def e : V := D.A - D.Rhat

lemma var_Rhat : D.C.var D.Rhat = D.C.var D.R + D.C.var D.u := by
  unfold Cov.var Rhat
  rw [D.C.add_left, D.C.add_right, D.C.add_right, D.R_u, D.C.symm D.u D.R, D.R_u]
  ring

/-- **The expansion.** This is the step the paper's derivation turns on, so it is proved, not assumed:
the covariance of the contrast with the measured reference is `(γ - 1) var R - var u`. The response
arm's own error `v` does not appear. -/
theorem cov_e_Rhat :
    D.C.B D.e D.Rhat = (D.γ - 1) * D.C.var D.R - D.C.var D.u := by
  unfold e A Rhat Cov.var
  rw [D.C.sub_left, D.C.add_left, D.C.add_left, D.C.smul_left, D.C.add_right, D.C.add_right,
      D.C.add_right, D.R_u, D.C.symm D.v D.R, D.R_v, D.C.symm D.v D.u, D.u_v,
      D.C.symm D.u D.R, D.R_u]
  ring

/-- Reliability of the reference arm: the share of its measured variance that is signal. -/
noncomputable def h : ℝ := D.C.var D.R / D.C.var D.Rhat

/-- The ordinary least-squares slope of the contrast on the measured reference. -/
noncomputable def slope : ℝ := D.C.B D.e D.Rhat / D.C.var D.Rhat

/-- The estimate the analyst reads off, `1 + slope`. -/
noncomputable def naive : ℝ := 1 + D.slope

/-- **The slope in the paper's form**: `(γ - 1)·h - (1 - h)`. -/
theorem slope_eq (hv : D.C.var D.Rhat ≠ 0) :
    D.slope = (D.γ - 1) * D.h - (1 - D.h) := by
  unfold slope h
  rw [D.cov_e_Rhat, D.var_Rhat] at *
  have hR : D.C.var D.R = D.C.var D.Rhat - D.C.var D.u := by rw [D.var_Rhat]; ring
  field_simp [D.var_Rhat] at hv ⊢

/-- **The central result.** The naive estimate is the gain multiplied by the reliability. -/
theorem naive_eq (hv : D.C.var D.Rhat ≠ 0) : D.naive = D.γ * D.h := by
  unfold naive
  rw [D.slope_eq hv]
  unfold h
  field_simp
  ring

/-- **Corollary 1.** With no divergence at all, the estimator returns the reliability. Deterministically:
this is an identity, not a tail event. -/
theorem corollary_one (hv : D.C.var D.Rhat ≠ 0) (hγ : D.γ = 1) : D.naive = D.h := by
  rw [D.naive_eq hv, hγ, one_mul]

/-- **Corollary 3.** The direction of the inference inverts exactly below the threshold `h < 1/γ`. -/
theorem corollary_three (hv : D.C.var D.Rhat ≠ 0) (hγ : 0 < D.γ) :
    D.naive < 1 ↔ D.h < 1 / D.γ := by
  rw [D.naive_eq hv]
  rw [lt_div_iff hγ]
  constructor
  · intro hlt; linarith [hlt]
  · intro hlt; linarith [hlt]

/-- **The coupling is what removes conservatism.** Dilution alone gives `γh + (1 - h)`, which for a true
increase lies strictly between 1 and γ: it can attenuate an effect but can never cross 1, so it can never
reverse the direction of a conclusion. -/
theorem dilution_only_cannot_cross_one
    (hh0 : 0 < D.h) (hh1 : D.h < 1) (hγ : 1 < D.γ) :
    1 < D.γ * D.h + (1 - D.h) ∧ D.γ * D.h + (1 - D.h) < D.γ := by
  constructor
  · nlinarith [mul_pos (sub_pos.mpr hγ) hh0]
  · nlinarith [mul_pos (sub_pos.mpr hh1) (sub_pos.mpr hγ)]

/-- **Identification.** An instrument correlated with the reference's signal but with neither error
recovers the gain, and the reliability does not enter. -/
theorem iv_identifies (S : V) (hSu : D.C.B S D.u = 0) (hSv : D.C.B S D.v = 0)
    (hRS : D.C.B D.R S ≠ 0) :
    1 + D.C.B D.e S / D.C.B D.Rhat S = D.γ := by
  have hden : D.C.B D.Rhat S = D.C.B D.R S := by
    unfold Rhat; rw [D.C.add_left, D.C.symm D.u S, hSu]; ring
  have hnum : D.C.B D.e S = (D.γ - 1) * D.C.B D.R S := by
    unfold e A Rhat
    rw [D.C.sub_left, D.C.add_left, D.C.add_left, D.C.smul_left,
        D.C.symm D.v S, hSv, D.C.symm D.u S, hSu]
    ring
  rw [hden, hnum]
  field_simp

/-- **The independence of the two errors is load-bearing, and its failure is signed.** If `cov(u, v) = c`
the slope gains exactly `c / var(Rhat)`; a positive covariance therefore moves the estimate toward the
truth and makes the artefact look smaller than it is. -/
theorem cov_uv_shifts_slope (C : Cov V) (R u v : V) (γ c : ℝ)
    (hRu : C.B R u = 0) (hRv : C.B R v = 0) (huv : C.B u v = c) (hv : C.var (R + u) ≠ 0) :
    C.B ((γ • R + v) - (R + u)) (R + u) / C.var (R + u)
      = ((γ - 1) * C.var R - C.var u) / C.var (R + u) + c / C.var (R + u) := by
  have hnum : C.B ((γ • R + v) - (R + u)) (R + u)
      = (γ - 1) * C.var R - C.var u + c := by
    unfold Cov.var
    rw [C.sub_left, C.add_left, C.add_left, C.smul_left, C.add_right, C.add_right, C.add_right,
        hRu, C.symm v R, hRv, C.symm v u, huv, C.symm u R, hRu]
    ring
  rw [hnum]
  field_simp

end Design

/-! ### The replicate axis

Averaging `n` independent replicates of the reference divides its measurement variance by `n`. The
reliability then rises as `1 / (1 + τ/n)` with `τ = var u / var R`, so replication helps, but only
through `h`, and no finite `n` reaches `h = 1`. -/

theorem h_of_replicates (vR vu : ℝ) (n : ℝ) (hR : 0 < vR) (hu : 0 ≤ vu) (hn : 0 < n) :
    vR / (vR + vu / n) = 1 / (1 + (vu / vR) / n) := by
  have h1 : vR + vu / n ≠ 0 := by positivity
  have h2 : 1 + (vu / vR) / n ≠ 0 := by positivity
  field_simp

theorem h_lt_one_of_pos_error (vR vu n : ℝ) (hR : 0 < vR) (hu : 0 < vu) (hn : 0 < n) :
    vR / (vR + vu / n) < 1 := by
  have hpos : 0 < vu / n := by positivity
  rw [div_lt_one (by positivity)]
  linarith

/-! ### The same expression in the change-on-baseline parametrisation

Hayes (1988), reviewing Blomqvist (1977), gives the coupled slope as `βc = β·R - (1 - R)` with `R` the
baseline reliability, and the correction that inverts it as `β = (βc + 1 - R)/R`. Setting `β = γ - 1`
and `R = h` makes the first identical to `slope` above, and the second recovers `γ - 1`. Both are
proved here so that the paper's concession of priority is exact rather than approximate. -/

theorem hayes_slope_is_ours (γ h : ℝ) :
    (γ - 1) * h - (1 - h) = γ * h - 1 := by ring

theorem blomqvist_correction_inverts (γ h : ℝ) (hh : h ≠ 0) :
    (((γ - 1) * h - (1 - h)) + 1 - h) / h = γ - 1 := by
  field_simp
  ring

/-- And the structural null they state, `crude slope = R - 1` at no effect, is Corollary 1 in the other
parametrisation. -/
theorem structural_null (h : ℝ) : (1 - 1 : ℝ) * h - (1 - h) = h - 1 := by ring

end Sharedstat
