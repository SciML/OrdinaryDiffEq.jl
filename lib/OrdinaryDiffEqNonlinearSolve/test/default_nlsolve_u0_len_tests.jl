using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: default_nlsolve, _init_autodiff
using NonlinearSolve
using SCCNonlinearSolve
using SciMLBase
using SciMLBase: OverrideInitData, OverrideInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using ForwardDiff
using LinearAlgebra: Diagonal
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

@testset "default solver differentiates a least-squares problem returned by the updater (outer iip = $iip)" for iip in (true, false)
    saw_dual = Ref(false)
    residual!(r, u, p) = (saw_dual[] |= eltype(u) <: ForwardDiff.Dual; r[1] = u[1] - 2; nothing)
    calls = Ref(0)
    updated = NonlinearLeastSquaresProblem(
        NonlinearFunction{true, SciMLBase.AutoSpecialize}(residual!; resid_prototype = zeros(1)), [1.5]
    )
    placeholder = NonlinearLeastSquaresProblem((u, p) -> [1.0], nothing)
    data = OverrideInitData(
        placeholder, (_, _) -> (calls[] += 1; updated), sol -> copy(sol.u), nothing;
        is_update_oop = Val(true)
    )
    rhs = iip ? (du, u, p, t) -> (du .= 2 .- u; nothing) : (u, p, t) -> 2 .- u
    fn = ODEFunction{iip, SciMLBase.FullSpecialize}(rhs; mass_matrix = zeros(1, 1), initialization_data = data)
    integ = init(ODEProblem(fn, [1.5], (0.0, 0.01)), Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test calls[] == 1
    @test saw_dual[]
    @test integ.sol.retcode == SciMLBase.ReturnCode.Default
    @test abs(only(integ.u) - 2) < 1.0e-10
    r = zeros(1)
    residual!(r, integ.u, nothing)
    @test abs(only(r)) < 1.0e-10
end

@testset "default solver accepts an interval problem returned by the initializer updater" begin
    calls = Ref(0)
    placeholder = NonlinearLeastSquaresProblem((u, p) -> [1.0], nothing)
    bracket = IntervalNonlinearProblem((x, p) -> x^2 - 2, (1.0, 2.0))
    data = OverrideInitData(
        placeholder, (_, _) -> (calls[] += 1; bracket), sol -> [sol.u], nothing;
        is_update_oop = Val(true)
    )
    f!(du, u, p, t) = (du[1] = u[1]^2 - 2; nothing)
    fn = ODEFunction{true, SciMLBase.FullSpecialize}(f!; mass_matrix = zeros(1, 1), initialization_data = data)
    prob = ODEProblem(fn, [1.3], (0.0, 0.01))
    integ = init(prob, Rodas5P(); abstol = 1.0e-12, reltol = 1.0e-12)
    @test calls[] == 1
    @test integ.sol.retcode != SciMLBase.ReturnCode.InitialFailure
    @test only(integ.u) ≈ sqrt(2) atol = 1.0e-10
    @test abs(only(integ.u)^2 - 2) < 1.0e-10
end

const explicit = (Returns(nothing), Returns(nothing))
square(u, p) = u .^ 2 .- 4
square!(du, u, p) = (du .= u .^ 2 .- 4; nothing)
square_jac(u, p) = Diagonal(2 .* u)
square_jac!(J, u, p) = (J .= Diagonal(2 .* u); nothing)

@testset "plain and SCC defaults start at Broyden" begin
    oop = NonlinearProblem(square, [1.5, 1.5])
    iip = NonlinearProblem(square!, [1.5, 1.5])
    oop_jac = NonlinearProblem(NonlinearFunction(square; jac = square_jac), [1.5, 1.5])
    iip_jac = NonlinearProblem(NonlinearFunction(square!; jac = square_jac!), [1.5, 1.5])
    plain = ((oop, [1.5, 1.5]), (iip, [1.5, 1.5]), (oop_jac, [1.5, 1.5]), (iip_jac, [1.5, 1.5]))
    sccs = Tuple(
        (SCCNonlinearProblem(blocks, explicit, nothing, Val(true)), zeros(4))
            for blocks in ((oop, iip), (oop_jac, iip), (oop_jac, iip_jac))
    )
    for (prob, u) in (plain..., sccs...), isinplace in (Val(true), Val(false))
        alg = default_nlsolve(nothing, isinplace, u, prob, false)
        @test alg.start_index == 1
        @test first(alg.algs).name === :Broyden
        @test alg.algs[end].autodiff isa AutoFiniteDiff
        @test default_nlsolve(nothing, isinplace, u, prob, true).algs[end].autodiff isa
            AutoForwardDiff
    end
end

@testset "SCC default keeps OrdinaryDiffEqTag and solves each block" begin
    auto = NonlinearProblem(square!, [1.5, 1.5])
    full = NonlinearProblem(NonlinearFunction{true, SciMLBase.FullSpecialize}(square!), [1.5, 1.5])
    for blocks in ((full, full), (full, auto))
        prob = SCCNonlinearProblem(blocks, explicit, nothing, Val(true))
        alg = default_nlsolve(nothing, Val(true), zeros(4), prob, true)
        @test tagof(alg.algs[end].autodiff) <: ForwardDiff.Tag
        for block in prob.probs
            sol = solve(block, alg; abstol = 1.0e-12)
            @test SciMLBase.successful_retcode(sol)
            @test sol.u ≈ [2.0, 2.0]
        end
    end
end

function override_initialization(initprob, u0)
    data = OverrideInitData(initprob, nothing, sol -> copy(sol.u), nothing)
    f = ODEFunction((du, u, p, t) -> (du .= -u; nothing); initialization_data = data)
    return init(ODEProblem(f, u0, (0.0, 1.0)), Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
end

@testset "default initialization solves a cubic from a zero-derivative guess" begin
    cubic! = (r, u, p) -> (r .= u .^ 3 .- 1; nothing)
    block = NonlinearProblem(NonlinearFunction{true, SciMLBase.FullSpecialize}(cubic!), [0.0])
    scc = SCCNonlinearProblem((block, block), explicit, nothing, Val(false))
    for (prob, n) in ((block, 1), (scc, 2))
        integ = override_initialization(prob, zeros(n))
        @test integ.sol.retcode == SciMLBase.ReturnCode.Default
        @test all(abs.(integ.u .- 1) .< 1.0e-10)
        @test all(abs.(integ.u .^ 3 .- 1) .< 1.0e-10)
    end
end

function strict_init_residual!(r::Vector{Float64}, u::Vector{Float64}, p)
    r .= u .- 2
    return nothing
end

@testset "SCC initialization accepts different residual requirements" begin
    raw! = (du, u, p, t) -> (du .= u .- 2; nothing)
    ode = ODEFunction{true, SciMLBase.AutoSpecialize}(raw!)
    wrapped = init(ODEProblem(ode, [1.5], (0.0, 1.0)), Rodas5P()).f
    wblock = NonlinearProblem(
        NonlinearFunction{true, SciMLBase.FullSpecialize}((du, u, p) -> wrapped(du, u, p, 0.0)), [1.5]
    )
    ablock = NonlinearProblem(
        NonlinearFunction{true, SciMLBase.AutoSpecialize}((du, u, p) -> (du .= u .- 3; nothing)), [1.5]
    )
    strict = NonlinearProblem(
        NonlinearFunction{true, SciMLBase.FullSpecialize}(strict_init_residual!), [1.5]
    )
    for (blocks, expected) in (((wblock, ablock), [2.0, 3.0]), ((strict, strict), [2.0, 2.0]))
        prob = SCCNonlinearProblem(blocks, explicit, nothing, Val(false))
        integ = override_initialization(prob, [1.5, 1.5])
        @test integ.sol.retcode == SciMLBase.ReturnCode.Default
        @test maximum(abs, integ.u .- expected) <= 1.0e-10
    end
end

# Smooth at the guess 2.0 and at the root 1.5, but a forward difference from 2.0 leaves the
# domain of the logarithm.
const log_c = 2.0 + 1.0e-9
log_resid!(r, u, p) = (r[1] = log(log_c - 1.5) - log(log_c - u[1]); nothing)
log_jac!(J, u, p) = (J[1, 1] = 1 / (log_c - u[1]); nothing)
# An undamped Newton step from 3.0 lands at 3 - 3log(3) < 0, outside the domain of `log`.
plain_log!(r, u, p) = (r[1] = log(u[1]); nothing)
plain_log_jac!(J, u, p) = (J[1, 1] = 1 / u[1]; nothing)
resid_at(resid!, x) = (r = zeros(1); resid!(r, [x], nothing); r[1])

@testset "default initialization with domain-limited residuals" begin
    cases = (
        (log_resid!, nothing, 2.0, 1.5),
        (log_resid!, log_jac!, 2.0, 1.5),
        (plain_log!, plain_log_jac!, 3.0, 1.0),
    )
    for (resid!, jac, guess, root) in cases
        block = NonlinearProblem(NonlinearFunction{true, SciMLBase.FullSpecialize}(resid!; jac), [guess])
        scc = SCCNonlinearProblem((block, block), explicit, nothing, Val(false))
        for (prob, n) in ((block, 1), (scc, 2))
            integ = override_initialization(prob, fill(guess, n))
            @test integ.sol.retcode == SciMLBase.ReturnCode.Default
            @test all(abs.(integ.u .- root) .<= 1.0e-10)
            @test all(abs(resid_at(resid!, x)) < 1.0e-10 for x in integ.u)
        end
    end
end
