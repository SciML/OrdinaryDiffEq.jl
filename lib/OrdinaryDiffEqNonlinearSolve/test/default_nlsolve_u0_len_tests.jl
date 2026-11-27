using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: default_nlsolve, _init_autodiff
using NonlinearSolve
using SciMLBase
using SciMLBase: OverrideInitData, OverrideInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using ForwardDiff
using Test

@testset "least-squares initialization solver" begin
    f(u, p) = [u[1]^2 - 4, u[2] - 1, u[1] - 2]
    prob = NonlinearLeastSquaresProblem(NonlinearFunction(f; resid_prototype = zeros(3)), [1.5, 1.5])
    sol = solve(prob, default_nlsolve(nothing, Val(false), prob.u0, prob); abstol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u ≈ [2.0, 1.0] atol = 1.0e-8
    @test f(sol.u, nothing) ≈ zeros(3) atol = 1.0e-8
end

tagof(::AutoForwardDiff{CS, T}) where {CS, T} = T

@testset "ForwardDiff tag matches the init problem's specialization" begin
    resid!(du, u, p) = (du .= u .^ 2 .- 4; nothing)
    auto = NonlinearProblem(resid!, [1.5, 1.5])  # AutoSpecialize by default
    full = NonlinearProblem(NonlinearFunction{true, SciMLBase.FullSpecialize}(resid!), [1.5, 1.5])
    nospec = NonlinearProblem(NonlinearFunction{true, SciMLBase.NoSpecialize}(resid!), [1.5, 1.5])
    oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, [1.5, 1.5])
    u = [1.5, 1.5]

    @test _init_autodiff(u, auto, false, Val(1)) isa AutoFiniteDiff
    @test _init_autodiff(u, full, false, Val(1)) isa AutoFiniteDiff

    for prob in (full, nospec, oop)
        @test tagof(_init_autodiff(u, prob, true, Val(1))) <: ForwardDiff.Tag
    end
    @test tagof(_init_autodiff(u, auto, true, Val(1))) === Nothing
    @test _init_autodiff(u, auto, true, Val(nothing)) isa AutoForwardDiff{1}
    @test _init_autodiff(u, auto, true, Val(0)) isa AutoForwardDiff{1}
    @test _init_autodiff(u, auto, true, Val(2)) isa AutoForwardDiff{2}

    for prob in (auto, full)
        sol = solve(prob, default_nlsolve(nothing, Val(true), u, prob, true); abstol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol.u ≈ [2.0, 2.0]
    end
end

@testset "default OverrideInit reaches the Jacobian step for $specialization" for specialization in (
        SciMLBase.AutoSpecialize, SciMLBase.NoSpecialize, SciMLBase.FullSpecialize,
    )
    saw_dual = Ref(false)
    function linear!(r, u, p, t)
        saw_dual[] |= eltype(u) <: ForwardDiff.Dual
        r[1] = 3u[1] + u[2] - 7
        r[2] = 2u[1] + 5u[2] - 12
        return nothing
    end
    root = [23, 22] ./ 13
    # NoSpecialize and FullSpecialize residuals are not wrapped by NonlinearSolve, so they
    # can call an ODE function wrapper that only accepts OrdinaryDiffEqTag Duals.
    wrapped = ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(linear!, [1.1, 1.2], (0.0, 0.01)).f
    rhs! = specialization === SciMLBase.AutoSpecialize ? linear! : wrapped
    residual!(r, u, p) = rhs!(r, u, p, 0.0)
    iprob = NonlinearProblem(
        NonlinearFunction{true, specialization}(residual!), [1.1, 1.2]; maxiters = 1
    )
    data = OverrideInitData(iprob, nothing, sol -> copy(sol.u), nothing)
    f = ODEFunction(linear!; mass_matrix = zeros(2, 2), initialization_data = data)
    prob = ODEProblem(f, [1.1, 1.2], (0.0, 0.01))
    integ = init(prob, Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test saw_dual[]
    @test integ.sol.retcode != SciMLBase.ReturnCode.InitialFailure
    @test integ.u ≈ root atol = 1.0e-10
    r = zeros(2)
    linear!(r, integ.u, nothing, 0.0)
    @test maximum(abs, r) < 1.0e-10
end

function closure!(f, ang, (l1, l2, x, y))
    f[1] = l1 * cos(ang[1]) + l2 * cos(ang[2]) - x
    f[2] = l1 * sin(ang[1]) + l2 * sin(ang[2]) - y
    return nothing
end

const intended = [-0.7, -1.2]
const guess = intended .+ [0.005, -0.005]
const linkage_p = (
    1.0, 0.7,
    cos(intended[1]) + 0.7cos(intended[2]),
    sin(intended[1]) + 0.7sin(intended[2]),
)

@testset "explicit NewtonRaphson in OverrideInit selects the nearby assembly" begin
    data = OverrideInitData(
        NonlinearProblem(closure!, copy(guess), linkage_p), nothing, sol -> copy(sol.u), nothing
    )
    f = ODEFunction(
        (du, u, p, t) -> closure!(du, u, p); mass_matrix = zeros(2, 2),
        initialization_data = data
    )
    prob = ODEProblem(f, copy(guess), (0.0, 0.01), linkage_p)
    integ = init(prob, Rodas5P(); initializealg = OverrideInit(nlsolve = NewtonRaphson()), abstol = 1.0e-10, reltol = 1.0e-10)
    @test integ.u ≈ intended atol = 1.0e-8
end

@testset "explicit NewtonRaphson in BrownFullBasicInit selects the nearby assembly" begin
    f = ODEFunction((du, u, p, t) -> closure!(du, u, p); mass_matrix = zeros(2, 2))
    prob = ODEProblem(f, copy(guess), (0.0, 0.01), linkage_p)
    integ = init(
        prob, Rodas5P(); initializealg = BrownFullBasicInit(nlsolve = NewtonRaphson()), abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test integ.u ≈ intended atol = 1.0e-8
end

@testset "default initialization preserves root and domain" begin
    for (g, guess, expected) in ((x -> x * (x - 2) * (x - 10), 0.97, 0.0), (x -> log(x / 2), 6.0, 2.0))
        residual!(r, u, p) = (r[1] = g(u[1]); nothing)
        f!(du, u, p, t) = residual!(du, u, p)
        for initializealg in (nothing, BrownFullBasicInit(), ShampineCollocationInit())
            data = isnothing(initializealg) ? OverrideInitData(
                    NonlinearProblem(residual!, [guess]), nothing, sol -> copy(sol.u), nothing
                ) : nothing
            prob = ODEProblem(ODEFunction(f!; mass_matrix = zeros(1, 1), initialization_data = data), [guess], (0.0, 0.01))
            for ad in (AutoForwardDiff(), AutoFiniteDiff())
                kwargs = isnothing(initializealg) ? (;) : (; initializealg)
                integ = init(prob, Rodas5P(autodiff = ad); abstol = 1.0e-10, reltol = 1.0e-10, kwargs...)
                @test only(integ.u) ≈ expected atol = 1.0e-8
                @test abs(g(only(integ.u))) < 1.0e-8
            end
        end
    end
end

@testset "typed supplied initialization residual" begin
    typed_residual!(out::Vector{Float64}, u::Vector{Float64}, p) = (out[1] = u[1] - 2; nothing)
    for specialization in (SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize)
        initprob = NonlinearProblem(NonlinearFunction{true, specialization}(typed_residual!), [1.0])
        data = OverrideInitData(initprob, nothing, sol -> copy(sol.u), nothing)
        f!(du, u, p, t) = (du .= 2 .- u; nothing)
        prob = ODEProblem(ODEFunction(f!; mass_matrix = zeros(1, 1), initialization_data = data), [0.0], (0.0, 0.01))
        integ = init(prob, Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
        @test integ.u ≈ [2.0] atol = 1.0e-8
    end
end

@testset "default solver follows the problem returned by the initializer updater" begin
    g(x) = 0.5 - exp(-1.0e4 * (x - 1)^2)
    root = 1 + sqrt(log(2)) / 100
    residual!(r, u, p) = (r[1] = g(u[1]); nothing)
    f!(du, u, p, t) = residual!(du, u, p)
    supplied = NonlinearProblem(NonlinearFunction{true, SciMLBase.AutoSpecialize}(residual!), [1.01])
    placeholder = NonlinearProblem(NonlinearFunction{true, SciMLBase.FullSpecialize}(residual!), [1.01])
    for updated in (false, true)
        calls = Ref(0)
        data = if updated
            OverrideInitData(
                placeholder, (_, _) -> (calls[] += 1; supplied), sol -> copy(sol.u), nothing;
                is_update_oop = Val(true)
            )
        else
            OverrideInitData(supplied, nothing, sol -> copy(sol.u), nothing)
        end
        fn = ODEFunction{true, SciMLBase.FullSpecialize}(f!; mass_matrix = zeros(1, 1), initialization_data = data)
        prob = ODEProblem(fn, [1.01], (0.0, 0.01))
        integ = init(prob, Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
        @test calls[] == updated
        @test integ.sol.retcode != SciMLBase.ReturnCode.InitialFailure
        @test only(integ.u) ≈ root atol = 1.0e-8
        @test abs(g(only(integ.u))) < 1.0e-10
    end
end
