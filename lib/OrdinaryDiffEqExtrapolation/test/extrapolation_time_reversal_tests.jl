using OrdinaryDiffEqExtrapolation, OrdinaryDiffEqCore, Test

# An ODE integrated forwards and its time-mirrored counterpart integrated
# backwards are the same problem, so an adaptive solver must take the same steps
# either way. With a tspan symmetric about zero the mirror map `t -> -t` is exact
# in IEEE arithmetic and round-to-nearest is symmetric under negation, so every
# stage value and error estimate mirrors *bitwise*: the `|dt|` sequences must
# agree to the last bit, not merely approximately.
#
# The Hairer-Wanner order selector formed `dt_new[i] = integrator.dt / Q[i]` with
# the sign kept, so `work[i] = s[i] / dt_new[i]` was negative for a backward
# solve and every `work[a] < sigma * work[b]` comparison (sigma = 9//10 > 0) came
# out inverted. The Deuflhard selector takes `abs` first and was unaffected.

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

function lv(u, p, t)
    a, b, c, d = p
    return [a * u[1] - b * u[1] * u[2], -c * u[2] + d * u[1] * u[2]]
end
function lv!(du, u, p, t)
    a, b, c, d = p
    du[1] = a * u[1] - b * u[1] * u[2]
    du[2] = -c * u[2] + d * u[1] * u[2]
    return nothing
end

const CASES = (
    ("linear", lin, lin!, [0.5], [1.01], 0.5),
    ("lotka", lv, lv!, [1.0, 1.0], [1.5, 1.0, 3.0, 1.0], 5.0),
)

@testset "Extrapolation time-reversal symmetry" begin
    algs = (
        "ExtrapolationMidpointHairerWanner" => ExtrapolationMidpointHairerWanner(),
        "ImplicitHairerWannerExtrapolation" => ImplicitHairerWannerExtrapolation(),
        "ImplicitEulerExtrapolation" => ImplicitEulerExtrapolation(),
        "ImplicitEulerBarycentricExtrapolation" => ImplicitEulerBarycentricExtrapolation(),
        # Deuflhard family already selected on magnitudes; kept as a guard.
        "ExtrapolationMidpointDeuflhard" => ExtrapolationMidpointDeuflhard(),
        "ImplicitDeuflhardExtrapolation" => ImplicitDeuflhardExtrapolation(),
    )
    for (algname, alg) in algs, (name, f, f!, u0, p, T) in CASES
        fdts, bdts = reversal_dts(
            f, f!, u0, p, T, alg;
            inplace = false, abstol = 1.0e-8, reltol = 1.0e-8
        )
        @testset "$algname $name oop" begin
            @test fdts == bdts
        end
    end
end
