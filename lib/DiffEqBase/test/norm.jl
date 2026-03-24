using Test
using ForwardDiff: Dual, gradient, partials

using DiffEqBase: ODE_DEFAULT_NORM
const internalnorm = ODE_DEFAULT_NORM

val = rand(10)
par = rand(10)
u = Dual.(val, par)
reference(val, par) = sqrt((sum(abs2, val) + sum(abs2, par)) / (length(val) + length(par)))
dual_real = internalnorm(u, 1)
dual_dual = internalnorm(u, u[1])
@test reference(val, par) ≈ dual_real
@test reference(val, par) ≈ dual_dual
@test partials(dual_dual, 1) ≈ gradient(x -> internalnorm(x, x[1]), val)'par
