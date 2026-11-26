using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: default_nlsolve, _init_autodiff
using NonlinearSolve
using SCCNonlinearSolve
using SciMLBase
using SciMLBase: OverrideInitData, OverrideInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using ForwardDiff
import ArrayInterface
using Test

start_index(alg) = alg.start_index
newton_first(alg) = first(alg.algs).name === :NewtonRaphson

small_oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, [1.5, 1.5])
small_iip = NonlinearProblem((du, u, p) -> (du .= u .^ 2 .- 4; nothing), [1.5, 1.5])
big_oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, fill(1.5, 26))
edge_oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, fill(1.5, 25))

@testset "polyalg start follows the init problem size" begin
    for ad in (false, true)
        @test newton_first(default_nlsolve(nothing, Val(false), small_oop.u0, small_oop, ad)) == true
        @test newton_first(default_nlsolve(nothing, Val(true), small_iip.u0, small_iip, ad)) == true
        # NonlinearSolve's cutoff is `u0_len <= 25`
        @test newton_first(default_nlsolve(nothing, Val(false), edge_oop.u0, edge_oop, ad)) == true
        @test newton_first(default_nlsolve(nothing, Val(false), big_oop.u0, big_oop, ad)) == false
    end
end

@testset "size comes from the problem, not from `u`" begin
    # BrownFullBasicInit passes the whole state as `u` but only solves the algebraic part,
    # so a big ODE with a handful of algebraic variables should still start with Newton...
    @test newton_first(default_nlsolve(nothing, Val(true), zeros(100), small_iip)) == true
    @test newton_first(default_nlsolve(nothing, Val(false), zeros(100), small_oop)) == true
    # ...and a short `u` must not make a big problem look small.
    @test newton_first(default_nlsolve(nothing, Val(false), [0.0], big_oop)) == false
end

# Stand-in for a GPU array: anything that says scalar indexing is slow. Newton would need a
# dense ForwardDiff Jacobian of it up front, so those keep the Broyden-first order.
struct SlowIndexVector{T} <: AbstractVector{T}
    x::Vector{T}
end
Base.size(v::SlowIndexVector) = size(v.x)
Base.getindex(v::SlowIndexVector, i::Int) = v.x[i]
ArrayInterface.fast_scalar_indexing(::Type{<:SlowIndexVector}) = false

@testset "arrays without fast scalar indexing keep the old order" begin
    slow = NonlinearProblem((u, p) -> u .^ 2 .- 4, SlowIndexVector([1.5, 1.5]))
    @test newton_first(default_nlsolve(nothing, Val(false), slow.u0, slow)) == false
    # a slow `u` alone is enough, e.g. Brown's full state living on the GPU
    @test newton_first(default_nlsolve(nothing, Val(false), SlowIndexVector(zeros(4)), small_oop)) == false
end

const explicit = (Returns(nothing), Returns(nothing))

@testset "SCCNonlinearProblem is sized by its largest block" begin
    # The SCC solve hands the same algorithm to every block, so it starts at Newton only when
    # all of the blocks are small enough for it.
    small = SCCNonlinearProblem((small_oop, small_iip), explicit, nothing, Val(true))
    # two blocks of 25 are 50 unknowns in total, but each solve only sees 25
    edge = SCCNonlinearProblem((edge_oop, edge_oop), explicit, nothing, Val(true))
    big = SCCNonlinearProblem((small_oop, big_oop), explicit, nothing, Val(true))
    slow_block = NonlinearProblem((u, p) -> u .^ 2 .- 4, SlowIndexVector([1.5, 1.5]))
    slow = SCCNonlinearProblem((small_oop, slow_block), explicit, nothing, Val(true))
    for iip in (Val(true), Val(false)), ad in (false, true)
        @test newton_first(default_nlsolve(nothing, iip, zeros(4), small, ad)) == true
        @test newton_first(default_nlsolve(nothing, iip, zeros(50), edge, ad)) == true
        @test newton_first(default_nlsolve(nothing, iip, zeros(28), big, ad)) == false
        # one GPU block keeps the Broyden-first order for the whole problem...
        @test newton_first(default_nlsolve(nothing, iip, zeros(4), slow, ad)) == false
        # ...and so does a GPU `u`
        @test newton_first(default_nlsolve(nothing, iip, SlowIndexVector(zeros(4)), small, ad)) == false
    end
end

@testset "least-squares default is untouched" begin
    nlls = NonlinearLeastSquaresProblem(
        NonlinearFunction((u, p) -> [u[1]^2 - 4, u[2] - 1, u[1] - 2]; resid_prototype = zeros(3)),
        [1.5, 1.5]
    )
    alg = default_nlsolve(nothing, Val(false), nlls.u0, nlls)
    @test start_index(alg) == 1
    @test first(alg.algs).name === :GaussNewton
end

tagof(::AutoForwardDiff{CS, T}) where {CS, T} = T

@testset "ForwardDiff tag matches the init problem's specialization" begin
    resid!(du, u, p) = (du .= u .^ 2 .- 4; nothing)
    auto = NonlinearProblem(resid!, [1.5, 1.5])  # AutoSpecialize by default
    full = NonlinearProblem(NonlinearFunction{true, SciMLBase.FullSpecialize}(resid!), [1.5, 1.5])
    u = [1.5, 1.5]

    @test _init_autodiff(u, auto, false, Val(1)) isa AutoFiniteDiff
    @test _init_autodiff(u, full, false, Val(1)) isa AutoFiniteDiff

    # FullSpecialize residuals (what BrownFullBasicInit / ShampineCollocationInit build) call
    # the ODE's FunctionWrappers, so they keep the OrdinaryDiffEqTag.
    @test tagof(_init_autodiff(u, full, true, Val(1))) <: ForwardDiff.Tag
    # AutoSpecialize init problems get wrapped by NonlinearSolve for its own tag, which it
    # only stamps onto an untagged backend.
    @test tagof(_init_autodiff(u, auto, true, Val(1))) === Nothing
    # chunk size normalization is the same either way
    @test _init_autodiff(u, auto, true, Val(nothing)) isa AutoForwardDiff{1}
    @test _init_autodiff(u, auto, true, Val(0)) isa AutoForwardDiff{1}
    @test _init_autodiff(u, auto, true, Val(2)) isa AutoForwardDiff{2}

    # and the default solver actually goes through Newton's Jacobian on the wrapped problem
    for prob in (auto, full)
        sol = solve(prob, default_nlsolve(nothing, Val(true), u, prob, true); abstol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol.u ≈ [2.0, 2.0]
    end

    full_scc = SCCNonlinearProblem((full, full), explicit, nothing, Val(true))
    mixed_scc = SCCNonlinearProblem((full, auto), explicit, nothing, Val(true))
    for prob in (full_scc, mixed_scc)
        alg = default_nlsolve(nothing, Val(true), zeros(4), prob, true)
        @test newton_first(alg)
        @test first(alg.algs).autodiff isa AutoFiniteDiff
        @test alg.algs[2].start_index == 1
        for block in prob.probs
            sol = solve(block, alg; abstol = 1.0e-12)
            @test SciMLBase.successful_retcode(sol)
            @test sol.u ≈ [2.0, 2.0]
        end
    end
end

# The reproducer from the issue: a closed two-rod linkage between fixed anchors, unknowns
# are the absolute rod angles. The guess is 0.005 rad away from the intended assembly, and
# the closure equations have other (far away) solutions that Broyden used to land on.
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

@testset "OverrideInit converges to the nearby assembly (#5005)" begin
    data = OverrideInitData(
        NonlinearProblem(closure!, copy(guess), linkage_p), nothing, sol -> copy(sol.u), nothing
    )
    f = ODEFunction(
        (du, u, p, t) -> closure!(du, u, p); mass_matrix = zeros(2, 2),
        initialization_data = data
    )
    prob = ODEProblem(f, copy(guess), (0.0, 0.01), linkage_p)
    integ = init(prob, Rodas5P(); initializealg = OverrideInit(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test integ.u ≈ intended atol = 1.0e-8
end

@testset "BrownFullBasicInit converges to the nearby assembly (#5005)" begin
    f = ODEFunction((du, u, p, t) -> closure!(du, u, p); mass_matrix = zeros(2, 2))
    prob = ODEProblem(f, copy(guess), (0.0, 0.01), linkage_p)
    integ = init(
        prob, Rodas5P(); initializealg = BrownFullBasicInit(), abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test integ.u ≈ intended atol = 1.0e-8
end

function override_initialization(nlprob, u0; nlsolve = nothing)
    data = OverrideInitData(nlprob, nothing, sol -> copy(sol.u), nothing)
    f = ODEFunction((du, u, p, t) -> (du .= -u; nothing); initialization_data = data)
    prob = ODEProblem(f, u0, (0.0, 1.0))
    return init(
        prob, Rodas5P(); initializealg = OverrideInit(; nlsolve),
        abstol = 1.0e-10, reltol = 1.0e-10
    )
end

@testset "Newton initialization retains Jacobian-free fallback" begin
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
        @test integ.u ≈ expected atol = 1.0e-10
        @test maximum(abs, integ.u .- expected) <= 1.0e-10
    end
end
