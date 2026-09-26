using StochasticDiffEq
using Test
using Random

# StochasticDiffEqCore resolves abstol/reltol to 1//10^2 before calling
# OrdinaryDiffEqCore._ode_init. Prove those SDE defaults are unchanged by the
# Core tolerance-resolution wrapper (which must pass concrete values through).

function expected_sde_default_tolerances(u0)
    uBottomEltype = eltype(u0)
    tol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
    return tol, tol
end

function linear_sde!(du, u, p, t)
    return du[1] = 1.01u[1]
end
function linear_noise!(du, u, p, t)
    return du[1] = 0.87u[1]
end

@testset "SDE default tolerances equal explicit 1e-2 (SOSRI, EM, SRIW1)" begin
    u0 = [1 / 2]
    tspan = (0.0, 1.0)
    prob = SDEProblem(linear_sde!, linear_noise!, u0, tspan)
    exp_a, exp_r = expected_sde_default_tolerances(u0)
    @test exp_a == 1.0e-2
    @test exp_r == 1.0e-2

    cases = (
        (SOSRI(), NamedTuple()),
        (EM(), (; dt = 0.01)),
        (SRIW1(), NamedTuple()),
    )
    for (alg, base_kw) in cases
        Random.seed!(12345)
        sol_default = solve(prob, alg; base_kw...)
        Random.seed!(12345)
        sol_explicit = solve(prob, alg; abstol = exp_a, reltol = exp_r, base_kw...)
        @test sol_default.stats.naccept == sol_explicit.stats.naccept
        @test sol_default.stats.nreject == sol_explicit.stats.nreject
        @test sol_default.t == sol_explicit.t
        @test sol_default.u == sol_explicit.u
        integ = init(prob, alg; base_kw...)
        @test integ.opts.abstol == exp_a
        @test integ.opts.reltol == exp_r
    end
end

@testset "complex user SDE tolerances are stored as real" begin
    # SDE entry applies real.(...) on user tolerances (master `_ode_init` did the
    # same) so adaptive residual calculations stay real.
    drift(u, p, t) = -u
    diffusion(u, p, t) = 0.1u
    for (u0, abstol, reltol) in (
            (1.0 + 0.5im, 1.0e-3 + 0im, 1.0e-3 + 0im),
            (1.0 + 0.5im, 1.0e-3 + 0.1im, 1.0e-3 + 0.2im),
            (1.0, 1.0e-3 + 0im, 1.0e-3 + 0im),
        )
        prob = SDEProblem(drift, diffusion, u0, (0.0, 0.1))
        for (alg, base_kw) in ((SOSRI(), (; dt = 0.01, seed = 1234)), (EM(), (; dt = 0.01, seed = 1234)))
            integ = init(prob, alg; abstol, reltol, base_kw...)
            @test integ.opts.abstol isa Real
            @test integ.opts.reltol isa Real
            @test integ.opts.abstol == real(abstol)
            @test integ.opts.reltol == real(reltol)
            if alg isa typeof(SOSRI())
                sol = solve!(integ)
                @test SciMLBase.successful_retcode(sol)
            end
        end
    end
end
