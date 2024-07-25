qmax_default(alg::Union{RadauIIA3, RadauIIA5, RadauIIA7}) = 8

alg_order(alg::RadauIIA3) = 3
alg_order(alg::RadauIIA5) = 5
alg_order(alg::RadauIIA7) = 7

isfirk(alg::RadauIIA3) = true
isfirk(alg::RadauIIA5) = true
isfirk(alg::RadauIIA7) = true

alg_adaptive_order(alg::RadauIIA3) = 1
alg_adaptive_order(alg::RadauIIA5) = 3
alg_adaptive_order(alg::RadauIIA7) = 5
