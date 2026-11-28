using OrdinaryDiffEqNordsieck, OrdinaryDiffEqCore, Test
import SciMLBase

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

@testset "AN5 rejection and callback accuracy" begin
    for inplace in (false, true)
        f = inplace ? ((du, u, p, t) -> (du .= u)) : ((u, p, t) -> u)
        for modified in (false, true)
            done = Ref(false)
            cb = DiscreteCallback(
                (u, t, it) -> !done[] && t >= 0.2,
                it -> (done[] = true; u_modified!(it, modified));
                save_positions = (false, false)
            )
            sol = solve(
                ODEProblem(f, [1.0], (0.0, 0.4)), AN5();
                dt = 0.01, callback = cb, tstops = [0.2], reltol = 1.0e-10, abstol = 1.0e-12
            )
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 0.4
            @test isapprox(sol.u[end][1], exp(0.4); rtol = 1.0e-6, atol = 1.0e-9)
        end
        f = inplace ? ((du, u, p, t) -> (du .= (1 + t) .* u)) : ((u, p, t) -> (1 + t) .* u)
        for dt in (1.0e-4, 0.01, 0.1)
            sol = solve(
                ODEProblem(f, [1.0], (0.0, 0.4)), AN5();
                dt, reltol = 1.0e-10, abstol = 1.0e-12
            )
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 0.4
            if dt == 0.1
                @test sol.stats.nreject > 0
            end
            @test isapprox(sol.u[end][1], exp(0.4 + 0.4^2 / 2); rtol = 1.0e-6, atol = 1.0e-9)
        end
    end
end

@testset "AN5 forced minimum timestep" begin
    for inplace in (false, true)
        f = inplace ? ((du, u, p, t) -> (du .= (1 + 20t^4) .* u)) :
            ((u, p, t) -> (1 + 20t^4) .* u)
        sol = solve(
            ODEProblem(f, [1.0], (0.0, 0.4)), AN5(); dt = 0.05, dtmin = 0.05,
            force_dtmin = true, reltol = 1.0e-8, abstol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(sol)
        @test isapprox(sol.u[end][1], exp(0.4 + 4 * 0.4^5); rtol = 1.0e-4)
    end
end
