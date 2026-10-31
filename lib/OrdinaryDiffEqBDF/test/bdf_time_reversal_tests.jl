using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, Test
using LinearAlgebra: mul!

# An ODE integrated forwards and its time-mirrored counterpart integrated
# backwards are the same problem, so an adaptive solver must take the same steps
# either way.  With a tspan symmetric about zero the mirror map `t -> -t` is
# exact in IEEE arithmetic and round-to-nearest is symmetric under negation, so
# every stage value, error estimate and Newton iterate mirrors *bitwise*: the
# `|dt|` sequences must agree to the last bit, not merely approximately.
#
# `step_accept_controller!` for QNDF/QBDF ranked the order-(k-1)/k/(k+1)
# candidate step sizes with `>` against a `0.0` sentinel meaning "this order is
# not viable", using the signed `integrator.dt`.  For a backward solve every
# candidate is negative, so the sentinel won every comparison: the order was
# driven to `max_order` immediately and `|dt|` never grew.  On `u' = 1.01u` that
# was 43 forward steps against 1075 backward.

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

lin(u, p, t) = p[1] .* u
lin!(du, u, p, t) = (du .= p[1] .* u; nothing)

const STIFFA = [-2000.0 1.0 0.0; 1.0 -3.0 1.0; 0.0 1.0 -0.5]
sys(u, p, t) = STIFFA * u
sys!(du, u, p, t) = (mul!(du, STIFFA, u); nothing)

const CASES = (
    ("linear", lin, lin!, [0.5], [1.01], 0.5),
    ("stiff linear system", sys, sys!, [1.0, 1.0, 1.0], [nothing], 1.0),
)

@testset "BDF time-reversal symmetry" begin
    # QBDF/QBDF2 are aliases of QNDF/QNDF2, so label the algorithms explicitly
    # rather than by type name.
    algs = (
        "QNDF" => QNDF(), "QNDF2" => QNDF2(),
        "QBDF" => QBDF(), "QBDF2" => QBDF2(), "FBDF" => FBDF(),
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
