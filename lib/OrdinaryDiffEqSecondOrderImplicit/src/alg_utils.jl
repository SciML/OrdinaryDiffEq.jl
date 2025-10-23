alg_extrapolates(alg::NewmarkBeta) = false

alg_order(alg::NewmarkBeta) = alg.γ ≈ 1 / 2 ? 2 : 1

is_mass_matrix_alg(alg::NewmarkBeta) = true

isadaptive(alg::NewmarkBeta) = true
