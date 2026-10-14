# Machine-checked derivations

Lean 4 (v4.14.0) with Mathlib. Fourteen theorems, no `sorry`, and every one depends only on Mathlib's
three standard axioms — `propext`, `Classical.choice`, `Quot.sound` — which `check_axioms.lean` verifies.
The absence of `sorryAx` from that list is the point: a file that merely compiles can still be full of
holes.

## What is proved

| theorem | statement |
|---|---|
| `var_Rhat` | var(R̂) = var(R) + var(u) |
| `cov_e_Rhat` | cov(e, R̂) = (γ − 1)·var(R) − var(u). **The step the paper turns on**, so it is proved rather than assumed |
| `slope_eq` | slope = (γ − 1)·h − (1 − h) |
| `naive_eq` | 1 + slope = γ·h |
| `corollary_one` | γ = 1 ⇒ the estimator returns h |
| `corollary_three` | γ̂ < 1 ⟺ h < 1/γ |
| `dilution_only_cannot_cross_one` | for 1 < γ and 0 < h < 1, dilution alone gives a value strictly between 1 and γ, so it can attenuate but never reverse |
| `iv_identifies` | an instrument orthogonal to both errors recovers γ, and h does not enter |
| `cov_uv_shifts_slope` | cov(u, v) = c shifts the slope by exactly c / var(R̂) |
| `h_of_replicates` | averaging n replicates gives h = 1/(1 + τ/n), τ = var(u)/var(R) |
| `h_lt_one_of_pos_error` | no finite n reaches h = 1 |
| `hayes_slope_is_ours` | Hayes' βc = β·R − (1 − R) is the same expression, with β = γ − 1 and R = h |
| `blomqvist_correction_inverts` | Blomqvist's correction recovers β from the coupled slope |
| `structural_null` | their "crude slope = R − 1 at no effect" is Corollary 1 reparametrised |

## What is deliberately not proved here

Covariance is modelled as a symmetric bilinear form and every orthogonality is an explicit hypothesis.
What the file settles is the **algebra**: given that second-order structure, each identity follows. It does
not settle the probabilistic model — that the errors are additive, or independent of the signal in a
particular dataset — because those are modelling assumptions about data, not theorems.

Corollary 2's "O(1) in the number of features" is also outside the scope. It is a statement about a
sample covariance converging on a population covariance as features are added, which is a limit theorem
about estimation, not an identity in the second-order algebra formalised here.

`E|e| = 2σ/√π` needs a distribution and so is absent. It is verified separately by symbolic integration
(sympy) and against the R implementation by simulation over fifteen (γ, h) cells at n = 4 × 10⁵, where
simulation, theorem and package agree to 3.1 × 10⁻³.

## Reproducing

```bash
lake exe cache get      # Mathlib prebuilt objects
lake build              # expect "Build completed successfully"
lean check_axioms.lean  # expect three axioms per theorem, and no sorryAx
```
