alg_extrapolates(alg::NewmarkBeta) = false
alg_extrapolates(alg::GeneralizedAlpha) = false

alg_order(alg::NewmarkBeta) = alg.γ ≈ 1 / 2 ? 2 : 1
# The step re-derives aₙ from f at the accepted state instead of carrying the
# algorithmic acceleration, so the velocity update is second order iff
# γ(1 - αf) = (1 - αm)/2 rather than on the textbook locus γ = 1/2 - αm + αf
# (measured fixed-step orders 2.000 and ~1 respectively for αm ≠ αf).
alg_order(alg::GeneralizedAlpha) = alg.γ * (1 - alg.αf) ≈ (1 - alg.αm) / 2 ? 2 : 1

is_mass_matrix_alg(alg::NewmarkBeta) = true
is_mass_matrix_alg(alg::GeneralizedAlpha) = true

isadaptive(alg::NewmarkBeta) = true
isadaptive(alg::GeneralizedAlpha) = true
