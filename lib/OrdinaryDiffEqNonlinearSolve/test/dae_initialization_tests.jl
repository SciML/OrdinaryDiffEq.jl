using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqBDF,
    StaticArrays, LinearAlgebra, Test, ADTypes

## Mass Matrix

function rober_oop(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du3 = y₁ + y₂ + y₃ - 1
    return [du1, du2, du3]
end
M = Diagonal([1.0, 1.0, 0.0])
f_oop = ODEFunction(rober_oop, mass_matrix = M)
prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = @inferred solve(
    prob_mm, Rosenbrock23(autodiff = AutoFiniteDiff()), reltol = 1.0e-8, abstol = 1.0e-8
)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = @inferred solve(
    prob_mm, Rosenbrock23(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = ShampineCollocationInit()
)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

integrator = @inferred init(prob_mm, Rodas5(autodiff = AutoForwardDiff(chunksize = 3)))
# It would be nice if this could test that the initialization is fully inferred,
# but since the return is just the integrator, this doesn't meet that goal.
@inferred SciMLBase.initialize_dae!(integrator)

prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = solve(
    prob_mm, Rosenbrock23(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
@test sum(sol.u[1]) ≈ 1
@test sol.u[1] ≈ [1.0, 0.0, 0.0]
for alg in [Rosenbrock23(autodiff = AutoFiniteDiff()), Trapezoid()]
    local sol
    sol = solve(
        prob_mm, alg, reltol = 1.0e-8, abstol = 1.0e-8,
        initializealg = ShampineCollocationInit()
    )
    @test sum(sol.u[1]) ≈ 1
end

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
M = Diagonal([1.0, 1.0, 0.0])
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = @inferred solve(prob_mm, Rodas5(autodiff = AutoFiniteDiff()), reltol = 1.0e-8, abstol = 1.0e-8)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = solve(
    prob_mm, Rodas5(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = ShampineCollocationInit()
)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

integrator = @inferred init(prob_mm, Rodas5(autodiff = AutoForwardDiff(chunksize = 3)))
# It would be nice if this could test that the initialization is fully inferred,
# but since the return is just the integrator, this doesn't meet that goal.
@inferred SciMLBase.initialize_dae!(integrator)

prob_mm = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = solve(
    prob_mm, Rodas5(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
@test sum(sol.u[1]) ≈ 1
@test sol.u[1] ≈ [1.0, 0.0, 0.0]

for alg in [Rodas5(autodiff = AutoFiniteDiff()), Trapezoid()]
    local sol
    sol = solve(
        prob_mm, alg, reltol = 1.0e-8, abstol = 1.0e-8,
        initializealg = ShampineCollocationInit()
    )
    @test sum(sol.u[1]) ≈ 1
end

function rober_no_p(du, u, p, t)
    y₁, y₂, y₃ = u
    (k₁, k₂, k₃) = (0.04, 3.0e7, 1.0e4)
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end

function rober_oop_no_p(du, u, p, t)
    y₁, y₂, y₃ = u
    (k₁, k₂, k₃) = (0.04, 3.0e7, 1.0e4)
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du3 = y₁ + y₂ + y₃ - 1
    return [du1, du2, du3]
end

# test oop and iip ODE initialization with parameters without eltype/length
struct UnusedParam
end
for f in (
        ODEFunction(rober_no_p, mass_matrix = M), ODEFunction(rober_oop_no_p, mass_matrix = M),
    )
    local prob, probp
    prob = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0e5))
    probp = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0e5), UnusedParam)
    for initializealg in (ShampineCollocationInit(), BrownFullBasicInit())
        isapprox(
            init(prob, Rodas5(), abstol = 1.0e-10; initializealg).u,
            init(probp, Rodas5(), abstol = 1.0e-10; initializealg).u
        )
    end
end

# to test that we get the right NL solve we need a broken solver.
struct BrokenNLSolve <: SciMLBase.AbstractNonlinearAlgorithm
    BrokenNLSolve(; kwargs...) = new()
end
function SciMLBase.__solve(
        prob::NonlinearProblem,
        alg::BrokenNLSolve, args...;
        kwargs...
    )
    u = fill(reinterpret(Float64, 0xDEADBEEFDEADBEEF), 3)
    return SciMLBase.build_solution(
        prob, alg, u, copy(u);
        retcode = ReturnCode.Success
    )
end
function f2(u, p, t)
    return u
end
f = ODEFunction(f2, mass_matrix = Diagonal([1.0, 1.0, 0.0]))
prob = ODEProblem(f, ones(3), (0.0, 1.0))
integrator = init(
    prob, Rodas5P(),
    initializealg = ShampineCollocationInit(1.0, BrokenNLSolve())
)
@test all(isequal(reinterpret(Float64, 0xDEADBEEFDEADBEEF)), integrator.u)

@testset "`reinit!` reruns initialization" begin
    initializeprob = NonlinearProblem(1.0, [0.0]) do u, p
        return u^2 - p[1]^2
    end
    initializeprobmap = function (nlsol)
        return [nlsol.prob.p[1], nlsol.u]
    end
    update_initializeprob! = function (iprob, integ)
        iprob.p[1] = integ.u[1]
    end
    initialization_data = SciMLBase.OverrideInitData(
        initializeprob, update_initializeprob!, initializeprobmap, nothing
    )
    fn = ODEFunction(; mass_matrix = [1 0; 0 0], initialization_data) do du, u, p, t
        du[1] = u[1]
        du[2] = u[1]^2 - u[2]^2
    end
    prob = ODEProblem(fn, [2.0, 0.0], (0.0, 1.0))
    integ = init(prob, Rodas5P())
    @test integ.u ≈ [2.0, 2.0] atol = 1.0e-8
    reinit!(integ)
    @test integ.u ≈ [2.0, 2.0] atol = 1.0e-8
    step!(integ, 0.01, true)
    @test SciMLBase.successful_retcode(integ.sol.retcode)
    reinit!(integ, reinit_dae = false)
    @test integ.u ≈ [2.0, 0.0]
    # With reinit_dae=false the algebraic constraint u[1]^2 - u[2]^2 = 0 is violated.
    # Rosenbrock methods (Rodas5P) linearize and don't iterate, so the step succeeds
    # but the constraint remains violated — u[2] stays near 0 instead of tracking u[1].
    step!(integ, 0.01, true)
    @test abs(integ.u[2]) < 1.0e-10  # u[2] stuck near 0, not reinitialized
    @test abs(integ.u[1]) > 1.5    # u[1] still evolving
end

# records the keywords Brown passes to the nonlinear solve
struct SpyNLSolve <: SciMLBase.AbstractNonlinearAlgorithm
    kwargs::Base.RefValue{Any}
end
SpyNLSolve() = SpyNLSolve(Ref{Any}((;)))
function SciMLBase.__solve(prob::NonlinearProblem, alg::SpyNLSolve, args...; kwargs...)
    alg.kwargs[] = NamedTuple(kwargs)
    return SciMLBase.build_solution(prob, alg, prob.u0, prob.u0; retcode = ReturnCode.Success)
end

@testset "BrownFullBasicInit tolerances" begin
    # y₁' = -y₁, 0 = y₂³ + y₂ - y₁
    mm_oop(u, p, t) = [-u[1], u[2]^3 + u[2] - u[1]]
    mm_iip(du, u, p, t) = (du .= mm_oop(u, p, t); nothing)
    dae_oop(du, u, p, t) = [du[1] + u[1], u[2]^3 + u[2] - u[1]]
    dae_iip(r, du, u, p, t) = (r .= dae_oop(du, u, p, t); nothing)
    y₂ = 0.6823278038280193 # root of y₂³ + y₂ = 1

    M = Diagonal([1.0, 0.0])
    dv = [true, false]
    make_probs(u0; kw...) = [
        (ODEProblem(ODEFunction(mm_oop; mass_matrix = M), u0, (0.0, 1.0); kw...), FBDF()),
        (ODEProblem(ODEFunction(mm_iip; mass_matrix = M), u0, (0.0, 1.0); kw...), FBDF()),
        (DAEProblem(dae_oop, [-1.0, 0.0], u0, (0.0, 1.0); differential_vars = dv, kw...), DFBDF()),
        (DAEProblem(dae_iip, [-1.0, 0.0], u0, (0.0, 1.0); differential_vars = dv, kw...), DFBDF()),
    ]

    # per-state tolerances are reduced to scalars for the nonlinear solve
    for (prob, alg) in make_probs([1.0, 2.0])
        for initializealg in (BrownFullBasicInit(), BrownFullBasicInit(abstol = nothing))
            integ = init(prob, alg; initializealg, abstol = [1.0e-8, 1.0e-8], reltol = [1.0e-8, 1.0e-8])
            @test integ.u ≈ [1.0, y₂]
        end
    end

    # the residual is below the solver's abstol but above Brown's default of 1e-10,
    # so only `abstol = nothing` leaves u0 untouched
    u0 = [1.0, y₂ + 1.0e-5]
    for (prob, alg) in make_probs(u0), abstol in (1.0e-3, [1.0e-3, 1.0e-3])
        integ = init(prob, alg; initializealg = BrownFullBasicInit(abstol = nothing), abstol)
        @test integ.u == u0
        integ = init(prob, alg; initializealg = BrownFullBasicInit(), abstol)
        @test integ.u ≈ [1.0, y₂]
        @test integ.u != u0
    end

    # the tolerances reach the nonlinear solve, with vectors reduced to their tightest entry
    for (prob, alg) in make_probs([1.0, 2.0])
        spy = SpyNLSolve()
        init(
            prob, alg; initializealg = BrownFullBasicInit(; abstol = 1.0e-6, nlsolve = spy),
            abstol = 1.0e-3, reltol = [1.0e-4, 1.0e-5]
        )
        @test get(spy.kwargs[], :abstol, nothing) == 1.0e-6
        @test get(spy.kwargs[], :reltol, nothing) == 1.0e-5

        spy = SpyNLSolve()
        init(
            prob, alg; initializealg = BrownFullBasicInit(; abstol = nothing, nlsolve = spy),
            abstol = [1.0e-2, 1.0e-3]
        )
        @test get(spy.kwargs[], :abstol, nothing) == 1.0e-3
    end

    # `abstol = nothing` also works when the algorithm is set on the problem
    spy = SpyNLSolve()
    initializealg = BrownFullBasicInit(; abstol = nothing, nlsolve = spy)
    for (prob, alg) in make_probs([1.0, 2.0]; initializealg)
        spy.kwargs[] = (;)
        init(prob, alg; abstol = 1.0e-3)
        @test get(spy.kwargs[], :abstol, nothing) == 1.0e-3
    end
end
