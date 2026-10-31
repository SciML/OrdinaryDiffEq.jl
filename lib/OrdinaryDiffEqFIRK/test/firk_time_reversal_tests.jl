using OrdinaryDiffEqFIRK, OrdinaryDiffEqCore, Test
using LinearAlgebra: mul!

# An ODE integrated forwards and its time-mirrored counterpart integrated
# backwards are the same problem, so an adaptive solver must take the same steps
# either way. With a tspan symmetric about zero the mirror map `t -> -t` is exact
# in IEEE arithmetic and round-to-nearest is symmetric under negation, so every
# stage value, Newton iterate and error estimate mirrors *bitwise*.
#
# `RadauIIA3Cache` used `needfactor = iter == 1`, so it handed `A = W1` to the
# linear solver on the first Newton iteration of every step, including steps
# where `new_W` was false and `W1` therefore still held the previous step's
# factorization. Re-factorizing those factors corrupts the solve; because an LU
# has an implicit unit diagonal the corruption is not sign-symmetric either.
# RadauIIA5/RadauIIA9/AdaptiveRadau all guard with `&& new_W`.

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

const STIFFA = [-2000.0 1.0 0.0; 1.0 -3.0 1.0; 0.0 1.0 -0.5]
sys(u, p, t) = STIFFA * u
sys!(du, u, p, t) = (mul!(du, STIFFA, u); nothing)

const CASES = (
    ("rober", rober, rober!, [1.0, 0.0, 0.0], [0.04, 3.0e7, 1.0e4], 50.0),
    ("stiff linear system", sys, sys!, [1.0, 1.0, 1.0], [nothing], 1.0),
)

@testset "FIRK time-reversal symmetry" begin
    algs = (
        "RadauIIA3" => RadauIIA3(), "RadauIIA5" => RadauIIA5(),
        "RadauIIA9" => RadauIIA9(), "AdaptiveRadau" => AdaptiveRadau(),
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
end
