using OrdinaryDiffEqDifferentiation, OrdinaryDiffEqCore, Test
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF, OrdinaryDiffEqDefault
using ADTypes: AutoFiniteDiff, AutoForwardDiff

# The state Jacobian differentiates with respect to `u`, so its finite-difference
# stencil has nothing to do with which way `t` is running. When `jacobian`/
# `jacobian!` took their `dir` from `diffdir(integrator)`, ∂f/∂u came out as a
# forward difference going forwards and a backward difference going backwards -
# two matrices that agree only to O(sqrt(eps)) - and a solve was no longer
# reproducible under time reversal. `solve(prob)` is affected because
# `prepare_alg(::Nothing, ...)` supplies
# `DefaultODEAlgorithm(autodiff = AutoFiniteDiff())`.
#
# On a tspan symmetric about zero the mirror map `t -> -t` is exact in IEEE
# arithmetic and round-to-nearest is symmetric under negation, so the mirrored
# problem `g(u, p, t) = -f(u, p, -t)` must reproduce the forward `|dt|` sequence
# bitwise.

"""Accepted `|dt|` sequence of `prob` under `alg`."""
function accepted_dts(prob, alg; kwargs...)
    integ = init(prob, alg; save_everystep = false, kwargs...)
    return [abs(Float64(integ.dt)) for _ in integ]
end

"""Forward and mirrored-backward `|dt|` sequences on tspan `(-T, T)`."""
function reversal_dts(f, f!, u0, p, T, alg; inplace = false, kwargs...)
    if inplace
        fwd = ODEProblem(f!, copy(u0), (-T, T), p)
        g! = (du, u, q, t) -> (f!(du, u, q, -t); du .= .-du; nothing)
        bwd = ODEProblem(g!, copy(u0), (T, -T), p)
    else
        fwd = ODEProblem(f, copy(u0), (-T, T), p)
        g = (u, q, t) -> -f(u, q, -t)
        bwd = ODEProblem(g, copy(u0), (T, -T), p)
    end
    return accepted_dts(fwd, alg; kwargs...), accepted_dts(bwd, alg; kwargs...)
end

function rober(u, p, t)
    k1, k2, k3 = p
    y1, y2, y3 = u
    return [-k1 * y1 + k3 * y2 * y3, k1 * y1 - k2 * y2^2 - k3 * y2 * y3, k2 * y2^2]
end
function rober!(du, u, p, t)
    k1, k2, k3 = p
    y1, y2, y3 = u
    du[1] = -k1 * y1 + k3 * y2 * y3
    du[2] = k1 * y1 - k2 * y2^2 - k3 * y2 * y3
    du[3] = k2 * y2^2
    return nothing
end

# Non-autonomous, so the time-derivative path is exercised too.
nonaut(u, p, t) = [cos(p[1] * t) * u[1] - 0.3 * u[2], sin(p[1] * t) * u[1]]
function nonaut!(du, u, p, t)
    du[1] = cos(p[1] * t) * u[1] - 0.3 * u[2]
    du[2] = sin(p[1] * t) * u[1]
    return nothing
end

const CASES = (
    ("rober", rober, rober!, [1.0, 0.0, 0.0], [0.04, 3.0e7, 1.0e4], 50.0),
    ("nonautonomous", nonaut, nonaut!, [1.0, 0.5], [3.0], 2.5),
)

@testset "finite-difference Jacobian direction" begin
    algs = (
        "Rodas5P/FiniteDiff" => Rodas5P(autodiff = AutoFiniteDiff()),
        "Rosenbrock23/FiniteDiff" => Rosenbrock23(autodiff = AutoFiniteDiff()),
        "FBDF/FiniteDiff" => FBDF(autodiff = AutoFiniteDiff()),
        # complex-step does not use `dir` at all, but the in-place path used to
        # pick between two prepared configs by `diffdir`
        "Rodas5P/FiniteDiff complex" => Rodas5P(
            autodiff = AutoFiniteDiff(fdtype = Val(:complex))
        ),
        # guards: forward-mode must stay reversible
        "Rodas5P/ForwardDiff" => Rodas5P(autodiff = AutoForwardDiff()),
    )
    for (algname, alg) in algs,
            (name, f, f!, u0, p, T) in CASES,
            inplace in (false, true)

        fdts, bdts = reversal_dts(
            f, f!, u0, p, T, alg;
            inplace, abstol = 1.0e-8, reltol = 1.0e-8
        )
        @testset "$algname $name $(inplace ? "iip" : "oop")" begin
            @test fdts == bdts
        end
    end

    # `solve(prob)` with no algorithm routes through
    # `DefaultODEAlgorithm(autodiff = AutoFiniteDiff())`, so it must be reversible too.
    @testset "default algorithm" begin
        for (name, f, f!, u0, p, T) in CASES, inplace in (false, true)
            fdts, bdts = reversal_dts(
                f, f!, u0, p, T, nothing;
                inplace, abstol = 1.0e-8, reltol = 1.0e-8
            )
            @testset "$name $(inplace ? "iip" : "oop")" begin
                @test fdts == bdts
            end
        end
    end
end
