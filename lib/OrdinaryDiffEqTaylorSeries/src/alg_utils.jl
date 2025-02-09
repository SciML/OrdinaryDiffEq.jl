alg_order(::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(::ExplicitTaylor{P}) where P = P
alg_stability_size(alg::ExplicitTaylor) = 1
