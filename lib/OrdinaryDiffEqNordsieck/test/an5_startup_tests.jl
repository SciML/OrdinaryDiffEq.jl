using OrdinaryDiffEqNordsieck, OrdinaryDiffEqCore, Test

# AN5 builds its Nordsieck history in one shot from a Tsit5 step, reading the
# higher derivatives out of the Tsit5 interpolant. `ode_interpolant`'s first
# argument is the normalised parameter Θ ∈ [0, 1], but the startup passed the
# absolute time `t`. The vector is assembled at the *start* of the step
# (`z[1] = uprev`, `z[2] = k[1] * dt`), so the correct value is `Θ = 0`.
#
# Every standard test problem starts at t = 0, where `t` happens to equal the
# right answer, which is why this went unnoticed. Away from the origin the
# startup derivatives are wrong, and since Θ then depends on where the interval
# sits, the method loses both time-translation invariance and time-reversal
# symmetry.

lin(u, p, t) = p[1] .* u
lin!(du, u, p, t) = (du .= p[1] .* u; nothing)

"""Accepted `|dt|` sequence of `prob` under `alg`."""
function accepted_dts(prob, alg; kwargs...)
    integ = init(prob, alg; save_everystep = false, kwargs...)
    return [abs(Float64(integ.dt)) for _ in integ]
end

@testset "AN5 startup" begin
    kw = (abstol = 1.0e-8, reltol = 1.0e-8, save_everystep = false)

    # An autonomous ODE does not care where the interval sits on the time axis.
    @testset "time-translation invariance" begin
        for f in (lin, lin!)
            u0 = [0.5]
            a = solve(ODEProblem(f, copy(u0), (0.0, 1.0), [1.01]), AN5(); kw...)
            b = solve(ODEProblem(f, copy(u0), (100.0, 101.0), [1.01]), AN5(); kw...)
            @test a.stats.naccept == b.stats.naccept
            @test a.u[end] ≈ b.u[end] rtol = 1.0e-6
        end
    end

    # On a tspan symmetric about zero the mirror map `t -> -t` is exact in IEEE
    # arithmetic and round-to-nearest is symmetric under negation, so the
    # mirrored problem `g(u, p, t) = -f(u, p, -t)` must reproduce the forward
    # `|dt|` sequence bitwise.
    @testset "time-reversal symmetry" begin
        for inplace in (false, true)
            if inplace
                fwd = ODEProblem(lin!, [0.5], (-0.5, 0.5), [1.01])
                g! = (du, u, p, t) -> (lin!(du, u, p, -t); du .= .-du; nothing)
                bwd = ODEProblem(g!, [0.5], (0.5, -0.5), [1.01])
            else
                fwd = ODEProblem(lin, [0.5], (-0.5, 0.5), [1.01])
                g = (u, p, t) -> -lin(u, p, -t)
                bwd = ODEProblem(g, [0.5], (0.5, -0.5), [1.01])
            end
            @testset "$(inplace ? "iip" : "oop")" begin
                @test accepted_dts(fwd, AN5(); kw...) == accepted_dts(bwd, AN5(); kw...)
            end
        end
    end
end
