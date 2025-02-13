using OrdinaryDiffEqDefault, OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner,
      OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF, ADTypes
using Test, LinearSolve, LinearAlgebra, SparseArrays, StaticArrays

f_2dlinear = (du, u, p, t) -> (@. du = p * u)

prob_ode_2Dlinear = ODEProblem(f_2dlinear, rand(4, 2), (0.0, 1.0), 1.01)
sol = @inferred solve(prob_ode_2Dlinear)

tsitsol = solve(prob_ode_2Dlinear, Tsit5())
# test that default is the same as Tsit5 (we expect it to use Tsit5 for this).
@test sol.stats.naccept == tsitsol.stats.naccept
@test sol.stats.nf == tsitsol.stats.nf
@test all(isequal(1), sol.alg_choice)
@test sol(0.5) == only(sol([0.5]).u) == tsitsol(0.5)
x = [zeros(4, 2) for _ in 1:5]
@test sol(x, 0:0.1:0.4) == tsitsol(x, 0:0.1:0.4)

sol = solve(prob_ode_2Dlinear, reltol = 1e-10)
vernsol = solve(prob_ode_2Dlinear, Vern7(), reltol = 1e-10)
# test that default is the same as Vern7 (we expect it to use Vern7 for this).
@test sol.stats.naccept == vernsol.stats.naccept
@test sol.stats.nf == vernsol.stats.nf
@test all(isequal(2), sol.alg_choice)
@test sol(0.5) == only(sol([0.5]).u) == vernsol(0.5)

prob_ode_linear_fast = ODEProblem(
    ODEFunction(f_2dlinear, mass_matrix = 2 * I(2)), rand(2), (0.0, 1.0), 1.01)
sol = solve(prob_ode_linear_fast)
@test all(isequal(4), sol.alg_choice)
# for some reason the timestepping here is different from regular Rosenbrock23 (including the initial timestep)

function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    [-k₁ * y₁ + k₃ * y₂ * y₃,
        k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
        k₂ * y₂^2]
end
prob_rober = ODEProblem(rober, [1.0, 0.0, 0.0], (0.0, 1e3), (0.04, 3e7, 1e4))
sol = solve(prob_rober)
rosensol = solve(prob_rober, AutoTsit5(Rosenbrock23(autodiff = AutoFiniteDiff())))
#test that cache is type stable
@test typeof(sol.interp.cache.cache3) == typeof(rosensol.interp.cache.caches[2])
# test that default has the same performance as AutoTsit5(Rosenbrock23()) (which we expect it to use for this).
@test sol.stats.naccept == rosensol.stats.naccept
@test sol.stats.nf == rosensol.stats.nf
@test unique(sol.alg_choice) == [1, 3]
@test sol.alg_choice[1] == 1
@test sol.alg_choice[end] == 3

sol = solve(prob_rober, reltol = 1e-7, abstol = 1e-7)
rosensol = solve(
    prob_rober, AutoVern7(Rodas5P(autodiff = AutoFiniteDiff())), reltol = 1e-7, abstol = 1e-7)
#test that cache is type stable
@test typeof(sol.interp.cache.cache4) == typeof(rosensol.interp.cache.caches[2])
# test that default has the same performance as AutoTsit5(Rosenbrock23()) (which we expect it to use for this).
@test sol.stats.naccept == rosensol.stats.naccept
@test sol.stats.nf == rosensol.stats.nf
@test unique(sol.alg_choice) == [2, 4]
@test sol.alg_choice[1] == 2
@test sol.alg_choice[end] == 4

function exrober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du .= vcat([-k₁ * y₁ + k₃ * y₂ * y₃,
            k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
            k₂ * y₂^2],
        fill(t, length(u) - 3))
end

for n in (100, 600)
    stiffalg = n < 50 ? 4 : n < 500 ? 5 : 6
    linsolve = n < 500 ? nothing : KrylovJL_GMRES()
    jac_prototype = sparse(I(n + 3))
    jac_prototype[1:3, 1:3] .= 1.0

    prob_ex_rober = ODEProblem(ODEFunction(exrober; jac_prototype),
        vcat([1.0, 0.0, 0.0], ones(n)), (0.0, 100.0), (0.04, 3e7, 1e4))
    global sol = solve(prob_ex_rober)
    fsol = solve(prob_ex_rober, AutoTsit5(FBDF(; autodiff = AutoFiniteDiff(), linsolve)))
    # test that default has the same performance as AutoTsit5(Rosenbrock23()) (which we expect it to use for this).
    @test sol.stats.naccept == fsol.stats.naccept
    @test sol.stats.nf == fsol.stats.nf
    @test unique(sol.alg_choice) == [1, stiffalg]
end

function swaplinear(u, p, t)
    [u[2], u[1]] .* p
end
swaplinearf = ODEFunction(swaplinear, mass_matrix = ones(2, 2) - I(2))
prob_swaplinear = ODEProblem(swaplinearf, rand(2), (0.0, 1.0), 1.01)
sol = solve(prob_swaplinear)
@test all(isequal(4), sol.alg_choice)
@test sol(0.5) isa Vector{Float64} # test dense output
# for some reason the timestepping here is different from regular Rodas5P (including the initial timestep)

# test mass matrix DAE where we have to initialize algebraic variables
function rober_mm(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
f = ODEFunction(rober_mm, mass_matrix = [1 0 0; 0 1 0; 0 0 0])
prob_rober_mm = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_rober_mm)
@test all(isequal(4), sol.alg_choice)
@test sol(0.5) isa Vector{Float64} # test dense output

# test callback on ConstantCache (https://github.com/SciML/OrdinaryDiffEq.jl/issues/2287)
using StaticArrays
cb = ContinuousCallback((u, t, integrator) -> t - 1, (integrator) -> nothing)
SA_ode_problem = ODEProblem((u, p, t) -> zero(u), SA[0], 2)
@test solve(SA_ode_problem; callback = cb).retcode == ReturnCode.Success

# test Complex numbers
H(s) = (1 - s) * complex([0 1; 1 0]) + s * complex([1 0; 0 -1])
schrod_eq(state, time, s) = -im * time * H(s) * state

prob_complex = ODEProblem(schrod_eq, complex([1, -1] / sqrt(2)), (0, 1), 100)
complex_sol = solve(prob_complex)
@test complex_sol.retcode == ReturnCode.Success
