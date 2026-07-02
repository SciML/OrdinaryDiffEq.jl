alg_extrapolates(alg::NewmarkBeta) = false
alg_extrapolates(alg::GeneralizedAlpha) = false

alg_order(alg::NewmarkBeta) = alg.γ ≈ 1 / 2 ? 2 : 1
alg_order(alg::GeneralizedAlpha) = alg.γ ≈ 1 / 2 - alg.αm + alg.αf ? 2 : 1

is_mass_matrix_alg(alg::NewmarkBeta) = true
is_mass_matrix_alg(alg::GeneralizedAlpha) = true

isadaptive(alg::NewmarkBeta) = true
isadaptive(alg::GeneralizedAlpha) = true
