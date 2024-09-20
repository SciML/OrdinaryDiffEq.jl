qmax_default(alg::Union{RadauIIA3, RadauIIA5, RadauIIA9}) = 8

alg_order(alg::RadauIIA3) = 3
alg_order(alg::RadauIIA5) = 5
alg_order(alg::RadauIIA9) = 9
alg_order(alg::AdaptiveRadau) = 5

isfirk(alg::RadauIIA3) = true
isfirk(alg::RadauIIA5) = true
isfirk(alg::RadauIIA9) = true
isfirk(alg::AdaptiveRadau) = true

alg_adaptive_order(alg::RadauIIA3) = 1
alg_adaptive_order(alg::RadauIIA5) = 3
alg_adaptive_order(alg::RadauIIA9) = 5
