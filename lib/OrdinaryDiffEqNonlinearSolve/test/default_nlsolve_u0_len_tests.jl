using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: default_nlsolve, _init_autodiff
using NonlinearSolve
using SciMLBase
using SciMLBase: OverrideInitData, OverrideInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using ForwardDiff
using Test

# SciML/OrdinaryDiffEq.jl#5005: the default initialization solver has to tell
# `FastShortcutNonlinearPolyalg` how big the system is, the same way `solve(prob, nothing)`
# does in NonlinearSolve. Small systems then start at NewtonRaphson instead of Broyden,
# which could wander off to a far-away root and still report Success.
#
# Starting at Newton means the init solve now takes a ForwardDiff Jacobian right away, so the
# backend's tag has to match how NonlinearSolve wrapped the init function.

start_index(alg) = alg.start_index

small_oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, [1.5, 1.5])
small_iip = NonlinearProblem((du, u, p) -> (du .= u .^ 2 .- 4; nothing), [1.5, 1.5])
big_oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, fill(1.5, 26))
edge_oop = NonlinearProblem((u, p) -> u .^ 2 .- 4, fill(1.5, 25))

@testset "polyalg start follows the init problem size" begin
    for ad in (false, true)
        @test start_index(default_nlsolve(nothing, Val(false), small_oop.u0, small_oop, ad)) == 3
        @test start_index(default_nlsolve(nothing, Val(true), small_iip.u0, small_iip, ad)) == 3
        # NonlinearSolve's cutoff is `u0_len <= 25`
        @test start_index(default_nlsolve(nothing, Val(false), edge_oop.u0, edge_oop, ad)) == 3
        @test start_index(default_nlsolve(nothing, Val(false), big_oop.u0, big_oop, ad)) == 1
    end
end

@testset "size comes from the problem, not from `u`" begin
    # BrownFullBasicInit passes the whole state as `u` but only solves the algebraic part,
    # so a big ODE with a handful of algebraic variables should still start with Newton...
    @test start_index(default_nlsolve(nothing, Val(true), zeros(100), small_iip)) == 3
    @test start_index(default_nlsolve(nothing, Val(false), zeros(100), small_oop)) == 3
    # ...and a short `u` must not make a big problem look small.
    @test start_index(default_nlsolve(nothing, Val(false), [0.0], big_oop)) == 1
end

@testset "SCCNonlinearProblem default is not changed here" begin
    # SCC init problems keep their current Broyden-first default for now; see the note next
    # to `_init_u0_len` in src/initialize_dae.jl.
    explicit = (Returns(nothing), Returns(nothing))
    scc = SCCNonlinearProblem((small_oop, small_oop), explicit, nothing, Val(true))
    for iip in (Val(true), Val(false))
        @test start_index(default_nlsolve(nothing, iip, zeros(4), scc)) == 1
    end
end

# Stand-in for a GPU array: anything that says scalar indexing is slow. Newton would need a
# dense ForwardDiff Jacobian of it up front, so those keep the Broyden-first order.
struct SlowIndexVector{T} <: AbstractVector{T}
    x::Vector{T}
end
Base.size(v::SlowIndexVector) = size(v.x)
Base.getindex(v::SlowIndexVector, i::Int) = v.x[i]
OrdinaryDiffEqNonlinearSolve.ArrayInterface.fast_scalar_indexing(::Type{<:SlowIndexVector}) = false

@testset "arrays without fast scalar indexing keep the old order" begin
    slow = NonlinearProblem((u, p) -> u .^ 2 .- 4, SlowIndexVector([1.5, 1.5]))
    @test start_index(default_nlsolve(nothing, Val(false), slow.u0, slow)) == 1
    # a slow `u` alone is enough, e.g. Brown's full state living on the GPU
    @test start_index(default_nlsolve(nothing, Val(false), SlowIndexVector(zeros(4)), small_oop)) == 1
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
