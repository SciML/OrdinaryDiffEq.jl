using OrdinaryDiffEqExtrapolation, OrdinaryDiffEqCore, Test

# The mutable and immutable code paths implement the same method, so they must
# agree on the solution and on how much work the solve takes.
#
# The order-extension block of the in-place ImplicitDeuflhardExtrapolation step
# rebuilt W with `jacobian2W!` but never handed the new matrix to the linear
# solver, so LinearSolve kept a factorization whose storage `jacobian2W!` had
# just overwritten; the midpoint residual in the same loop was also scaled by an
# extra factor of `dt_int`. Both were confined to T[n_curr+1], the only entry
# that block computes.

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

"""Accepted `|dt|` sequence of `prob` under `alg`."""
function accepted_dts(prob, alg; kwargs...)
    integ = init(prob, alg; save_everystep = false, kwargs...)
    return [abs(Float64(integ.dt)) for _ in integ]
end

@testset "Extrapolation in-place/out-of-place consistency" begin
    algs = (
        "ImplicitDeuflhardExtrapolation" => ImplicitDeuflhardExtrapolation(),
        "ImplicitHairerWannerExtrapolation" => ImplicitHairerWannerExtrapolation(),
    )
    kw = (abstol = 1.0e-8, reltol = 1.0e-8, save_everystep = false)

    @testset "same effort as out-of-place" begin
        for (algname, alg) in algs, (name, f, f!, u0, p, T) in CASES
            soop = solve(ODEProblem(f, copy(u0), (-T, T), p), alg; kw...)
            siip = solve(ODEProblem(f!, copy(u0), (-T, T), p), alg; kw...)
            @testset "$algname $name" begin
                @test siip.stats.naccept <= 2 * soop.stats.naccept
                @test siip.stats.nreject <= 2 * soop.stats.nreject + 4
                @test siip.u ≈ soop.u rtol = 1.0e-5
            end
        end
    end

    # A corrupted solve is also not sign-symmetric (an LU has an implicit unit
    # diagonal, so LU(-X) != -LU(X)), which cost the in-place method its
    # time-reversal symmetry: on a tspan symmetric about zero the mirrored
    # problem `g(u, p, t) = -f(u, p, -t)` must reproduce the forward |dt|
    # sequence bitwise.
    # Only the Deuflhard family here: the Hairer-Wanner order selector has its own
    # direction bug in `step_accept_controller!`, fixed separately.
    @testset "in-place time-reversal symmetry" begin
        for (algname, alg) in ("ImplicitDeuflhardExtrapolation" => ImplicitDeuflhardExtrapolation(),),
                (name, f, f!, u0, p, T) in CASES
            fwd = ODEProblem(f!, copy(u0), (-T, T), p)
            g! = (du, u, q, t) -> (f!(du, u, q, -t); du .= .-du; nothing)
            bwd = ODEProblem(g!, copy(u0), (T, -T), p)
            @testset "$algname $name" begin
                @test accepted_dts(fwd, alg; kw...) == accepted_dts(bwd, alg; kw...)
            end
        end
    end
end
