alg_order(alg::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(alg::ExplicitTaylor{P}) where P = P
alg_stability_size(alg::ExplicitTaylor) = 1

alg_order(alg::DAETS) = 2
alg_stability_size(alg::DAETS) = 1
