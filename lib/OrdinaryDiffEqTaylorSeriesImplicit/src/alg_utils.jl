alg_extrapolates(::ImplicitTaylor) = true

alg_order(::ImplicitTaylor{P}) where {P} = P

alg_adaptive_order(::ImplicitTaylor) = 0
