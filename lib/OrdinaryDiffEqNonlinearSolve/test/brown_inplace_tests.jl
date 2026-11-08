using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve, NonlinearSolve, ADTypes,
    ForwardDiff, LinearAlgebra, SciMLBase, Test

# du₁ = -u₁ + u₂, 0 = u₁ + u₂ - p. The consistent algebraic value is u₂ = p - u₁; then
# u₁' = p - 2u₁, so u₁(t) = p/2 + (u₁(0) - p/2) exp(-2t).
function brown_f!(du, u, p, t)
    du[1] = -u[1] + u[2]
    du[2] = u[1] + u[2] - p
    return nothing
end
const brown_fn = ODEFunction{true, SciMLBase.FullSpecialize}(
    brown_f!; mass_matrix = Diagonal([1.0, 0.0])
)
const brown_nlsolves = (
    nothing, NewtonRaphson(autodiff = AutoFiniteDiff()),
    NewtonRaphson(autodiff = AutoForwardDiff()),
)

@testset "In-place BrownFullBasicInit, nlsolve = $nl" for nl in brown_nlsolves
    ia = BrownFullBasicInit(nlsolve = nl)

    prob = ODEProblem(brown_fn, [1.0, 0.5], (0.0, 1.0), 1.0)
    integ = init(prob, Rodas5P(); initializealg = ia)
    @test integ.u ≈ [1.0, 0.0] atol = 1.0e-10
    @test prob.u0 == [1.0, 0.5]

    prob = ODEProblem(brown_fn, [1.0, 0.5], (0.0, 1.0), 1.0)
    integ = init(
        prob, Rodas5P(); initializealg = ia,
        alias = SciMLBase.ODEAliasSpecifier(alias_u0 = true)
    )
    @test integ.u === prob.u0
    @test prob.u0 ≈ [1.0, 0.0] atol = 1.0e-10

    prob = ODEProblem(brown_fn, [1.0, 0.5], (0.0, 1.0), 1.0)
    sol = solve(prob, Rodas5P(); initializealg = ia, abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    u1 = 0.5 + 0.5exp(-2.0)
    @test sol.u[end] ≈ [u1, 1 - u1] atol = 1.0e-7

    # Reinitialization after a callback that sets p = 2 at t = 0.5 changes only the
    # algebraic value, to 2 - u₁
    cb = DiscreteCallback(
        (u, t, integ) -> t == 0.5, integ -> (integ.p = 2.0); initializealg = ia
    )
    sol = solve(
        prob, Rodas5P(); initializealg = ia, callback = cb, tstops = [0.5],
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol)
    i = findall(==(0.5), sol.t)
    @test length(i) == 2
    @test sol.u[i[2]][1] == sol.u[i[1]][1]
    @test sol.u[i[2]][2] ≈ 2.0 - sol.u[i[1]][1] atol = 1.0e-10
    @test sol.u[i[1]][2] ≈ 1.0 - sol.u[i[1]][1] atol = 1.0e-8
end

@testset "In-place BrownFullBasicInit under ForwardDiff" begin
    # The initialized algebraic value is p - u₁(0) = p - 1
    function initialized_algebraic(p)
        prob = ODEProblem(brown_fn, [one(p), one(p) / 2], (0.0, 1.0), p)
        ia = BrownFullBasicInit(nlsolve = NewtonRaphson(autodiff = AutoForwardDiff()))
        return init(prob, Rodas5P(); initializealg = ia).u[2]
    end
    @test initialized_algebraic(1.5) ≈ 0.5 atol = 1.0e-10
    @test ForwardDiff.derivative(initialized_algebraic, 1.5) ≈ 1.0 atol = 1.0e-8
end

# u = [y, x, z] with algebraic y and z: 0 = y³ - x³, x' = -x, 0 = z³ - exp(3x). The
# consistent algebraic values are y = x and z = exp(x), and x(t) = x(0) exp(-t).
function cubic_f!(du, u, p, t)
    y, x, z = u
    du[1] = y^3 - x^3
    du[2] = -x
    du[3] = z^3 - exp(3x)
    return nothing
end
const cubic_fn = ODEFunction{true, SciMLBase.FullSpecialize}(
    cubic_f!; mass_matrix = Diagonal([0.0, 1.0, 0.0])
)

@testset "Nonlinear interleaved BrownFullBasicInit, nlsolve = $nl" for nl in brown_nlsolves
    ia = BrownFullBasicInit(nlsolve = nl)
    prob = ODEProblem(cubic_fn, [2.0, 1.0, 1.0], (0.0, 1.0))
    integ = init(prob, Rodas5P(); initializealg = ia)
    @test integ.u ≈ [1.0, 1.0, exp(1.0)] atol = 1.0e-8
    @test prob.u0 == [2.0, 1.0, 1.0]
    sol = solve(prob, Rodas5P(); initializealg = ia, abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ [exp(-1.0), exp(-1.0), exp(exp(-1.0))] atol = 1.0e-7
end

# Records its problem and the residual's argument types, then solves with NewtonRaphson
struct RecordingNewton <: SciMLBase.AbstractNonlinearAlgorithm
    probs::Vector{Any}
    residual_types::Vector{Type}
end
RecordingNewton() = RecordingNewton(Any[], Type[])
function SciMLBase.__solve(prob::NonlinearProblem, alg::RecordingNewton, args...; kwargs...)
    push!(alg.probs, prob)
    f! = (du, u, p) -> (push!(alg.residual_types, typeof(u)); prob.f(du, u, p))
    recording = NonlinearProblem(
        NonlinearFunction{true, SciMLBase.FullSpecialize}(f!), prob.u0, prob.p
    )
    return solve(recording, NewtonRaphson(autodiff = AutoFiniteDiff()), args...; kwargs...)
end

@testset "In-place BrownFullBasicInit solves in a dense Vector" begin
    rec = RecordingNewton()
    prob = ODEProblem(brown_fn, [1.0, 0.5], (0.0, 1.0), 1.0)
    integ = init(prob, Rodas5P(); initializealg = BrownFullBasicInit(nlsolve = rec))
    @test integ.u ≈ [1.0, 0.0] atol = 1.0e-10
    nlprob = only(rec.probs)
    @test nlprob.u0 isa Vector{Float64}
    @test nlprob.u0 == [0.5]
    @test unique(rec.residual_types) == [Vector{Float64}]
    for alias_u0 in (false, true)
        alias = SciMLBase.NonlinearAliasSpecifier(; alias_u0)
        nlalg = NewtonRaphson(autodiff = AutoFiniteDiff())
        cache = @inferred init(nlprob, nlalg; alias)
        @test SciMLBase.solve!(cache).u isa Vector{Float64}
    end
end

# Throws once the algebraic value leaves [0.4, 0.6], which the first Newton step (to 0) does
function throwing_brown_f!(du, u, p, t)
    abs(u[2] - 0.5) > 0.1 && throw(DomainError(u[2], "algebraic value left [0.4, 0.6]"))
    return brown_f!(du, u, p, t)
end

@testset "In-place BrownFullBasicInit that throws, nlsolve = $nl" for nl in brown_nlsolves
    fn = ODEFunction{true, SciMLBase.FullSpecialize}(
        throwing_brown_f!; mass_matrix = Diagonal([1.0, 0.0])
    )
    # [0.5, 0.5] is consistent for p = 1, so `init` doesn't solve
    prob = ODEProblem(fn, [0.5, 0.5], (0.0, 1.0), 1.0)
    integ = init(prob, Rodas5P(); initializealg = BrownFullBasicInit())
    integ.u[1] = 1.0
    saved_u = copy(integ.u)
    saved_uprev = copy(integ.uprev)
    @test_throws DomainError OrdinaryDiffEqNonlinearSolve._initialize_dae!(
        integ, integ.sol.prob, BrownFullBasicInit(nlsolve = nl), Val(true)
    )
    @test integ.u == saved_u
    @test integ.uprev == saved_uprev
end
