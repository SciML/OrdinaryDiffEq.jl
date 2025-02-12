using OrdinaryDiffEq, OrdinaryDiffEqCore, DiffEqBase, Test, ADTypes
using Random, SparseDiffTools
using OrdinaryDiffEqDefault
using ElasticArrays, LinearSolve
Random.seed!(213)
CACHE_TEST_ALGS = [Euler(), Midpoint(), RK4(), SSPRK22(), SSPRK33(), SSPRK43(), SSPRK104(),
    CarpenterKennedy2N54(), SHLDDRK64(), ORK256(), DGLDDRK73_C(),
    CFRLDDRK64(), TSLDDRK74(), CKLLSRK43_2(), ParsaniKetchesonDeconinck3S32(),
    BS3(), BS5(), DP5(), DP8(), Feagin10(), Feagin12(), Feagin14(), TanYam7(),
    Tsit5(), TsitPap8(), Vern6(), Vern7(), Vern8(), Vern9(), OwrenZen3(), OwrenZen4(),
    OwrenZen5(), AutoTsit5(Rosenbrock23()), TRBDF2(), KenCarp4(), ABDF2(),
    OrdinaryDiffEqDefault.DefaultODEAlgorithm()]
broken_CACHE_TEST_ALGS = [
    QNDF(),
    ExtrapolationMidpointHairerWanner(),
    ImplicitEulerExtrapolation(),
    ImplicitDeuflhardExtrapolation()
]
# AitkenNeville(threading=false) fails Elastic but not normal case

using InteractiveUtils

NON_IMPLICIT_ALGS = filter((x) -> isconcretetype(x) && !OrdinaryDiffEqCore.isimplicit(x()),
    union(subtypes(OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm),
        subtypes(OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveAlgorithm)))

f = function (du, u, p, t)
    for i in 1:length(u)
        du[i] = (0.3 / length(u)) * u[i]
    end
end

condition = function (u, t, integrator)
    1 - maximum(u)
end

affect! = function (integrator)
    u = integrator.u
    maxidx = findmax(u)[2]
    resize!(integrator, length(u) + 1)
    Θ = rand() / 5 + 0.25
    u[maxidx] = Θ
    u[end] = 1 - Θ
    nothing
end

callback = ContinuousCallback(condition, affect!)

u0 = [0.2]
tspan = (0.0, 10.0)
prob = ODEProblem(f, u0, tspan)

println("Check for stochastic errors")
for i in 1:10
    @test_nowarn sol = solve(prob, Tsit5(), callback = callback)
end

println("Check some other integrators")
sol = solve(prob, Rosenbrock23(), callback = callback, dt = 1 / 2)
@test length(sol[end]) > 1
sol = solve(prob, Rosenbrock32(), callback = callback, dt = 1 / 2)
@test length(sol[end]) > 1
sol = solve(prob, KenCarp4(), callback = callback, dt = 1 / 2)
@test length(sol[end]) > 1
sol = solve(prob, TRBDF2(), callback = callback, dt = 1 / 2)
@test length(sol[end]) > 1
sol = solve(prob, TRBDF2(linsolve = LinearSolve.KrylovJL_GMRES()),
    callback = callback)
@test length(sol[end]) > 1

for alg in CACHE_TEST_ALGS
    @show alg
    local sol = solve(prob, alg, callback = callback, dt = 1 / 2)
    @test length(sol[end]) > 1
end

for alg in broken_CACHE_TEST_ALGS
    @show alg
    @test_broken length(solve(prob, alg, callback = callback, dt = 1 / 2)[end]) > 1
end

sol = solve(prob, Rodas4(autodiff = AutoForwardDiff(chunksize = 1)),
    callback = callback, dt = 1 / 2)
@test length(sol[end]) > 1
sol = solve(prob, Rodas5(autodiff = AutoForwardDiff(chunksize = 1)),
    callback = callback, dt = 1 / 2)
@test length(sol[end]) > 1

# cache tests resizing multidimensional arrays
println("Check resizing multidimensional arrays")
u0_matrix = ElasticArray(ones(2, 2))
f_matrix = (du, u, p, t) -> du .= u
prob_matrix = ODEProblem(f_matrix, u0_matrix, (0.0, 2.0))
condition_matrix = (u, t, integrator) -> t - 1
affect_matrix! = function (integrator)
    resize!(integrator, (2, 3))
    integrator.u .= 1
    nothing
end
callback_matrix = ContinuousCallback(condition_matrix, affect_matrix!)

for alg in CACHE_TEST_ALGS
    OrdinaryDiffEqCore.isimplicit(alg) && continue # this restriction should be removed in the future
    @show alg
    local sol = solve(prob_matrix, alg, callback = callback_matrix, dt = 1 / 2)
    @test size(sol[end]) == (2, 3)
end

# additional cache tests to find more bugs
println("Additional resize! checks")
u0resize3 = ones(2)
fresize3 = (du, u, p, t) -> du .= u
prob_resize3 = ODEProblem(fresize3, u0resize3, (0.0, 2.0))
condition_resize3 = (u, t, integrator) -> t - 1
affect!_resize3 = function (integrator)
    resize!(integrator, 3)
    integrator.u .= 1
    nothing
end
callback_resize3 = ContinuousCallback(condition_resize3, affect!_resize3)

for alg in CACHE_TEST_ALGS
    (OrdinaryDiffEqCore.isimplicit(alg) || OrdinaryDiffEqCore.alg_order(alg) < 2) &&
        continue
    @show alg
    local sol = solve(prob_resize3, alg, callback = callback_resize3, dt = 0.125)
    @test size(sol[end]) == (3,)
    @test all(sol[end] .== sol[end][1])
    @test sol[end][1]≈exp(1) atol=1.0e-2
end

# additional cache tests to find more bugs
println("Enforced adaptive dt checks")
u0_adapt = zeros(2)
f_adapt = (du, u, p, t) -> du .= t
prob_adapt = ODEProblem(f_adapt, u0_adapt, (0.0, 2.0))
condition_adapt = (u, t, integrator) -> true
affect!_adapt = function (integrator)
    dt = rand((0.0625, 0.125, 0.25))
    set_proposed_dt!(integrator, dt)
    integrator.opts.dtmax = dt
    integrator.dtcache = dt
    u_modified!(integrator, false)
    nothing
end
callback_adapt = DiscreteCallback(condition_adapt, affect!_adapt,
    save_positions = (false, false))

for alg in CACHE_TEST_ALGS
    (OrdinaryDiffEqCore.isimplicit(alg) || OrdinaryDiffEqCore.alg_order(alg) < 2) &&
        continue
    @show alg
    local sol = solve(prob_adapt, alg, callback = callback_adapt, dt = 0.125)
    @test all(idx -> all(isapprox.(sol.u[idx], 0.5 * sol.t[idx]^2, atol = 1.0e-6)),
        eachindex(sol.t))
end

# Force switching

function f2(du, u, p, t)
    @assert length(u)==length(du) "length(u) = $(length(u)), length(du) = $(length(du)) at time $(t)"
    for i in 1:length(u)
        if t > 10
            du[i] = -10000 * u[i]
        else
            du[i] = 0.3 * u[i]
        end
    end
    return du
end

function condition2(u, t, integrator)
    1 - maximum(u)
end

function affect2!(integrator)
    u = integrator.u
    resize!(integrator, length(u) + 1)
    maxidx = findmax(u)[2]
    Θ = rand()
    u[maxidx] = Θ
    u[end] = 1 - Θ
    nothing
end

callback = ContinuousCallback(condition2, affect2!)
u0 = [0.2]
tspan = (0.0, 20.0)
prob = ODEProblem(f2, u0, tspan)
sol = solve(prob, AutoTsit5(Rosenbrock23()), callback = callback)
@test length(sol[end]) > 1
