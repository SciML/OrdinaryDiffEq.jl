using OrdinaryDiffEqDifferentiation, OrdinaryDiffEqBDF, DiffEqBase, SciMLBase, ADTypes, Test
using OrdinaryDiffEqDifferentiation: prepare_ADType

# A despecialized (wrapped) residual must be differentiated with chunk size 1 for every
# state length, for `DAEFunction` exactly as for `ODEFunction`; otherwise the chunk chosen
# from `length(u)` ends up in the nonlinear solver's Jacobian config and the whole DFBDF
# step recompiles for every new problem size.
function resid!(out, du, u, p, t)
    n = length(u)
    dx = 1 / (n - 1)
    out[1] = u[1]
    out[n] = u[n]
    @inbounds for i in 2:(n - 1)
        out[i] = du[i] - p.α * (u[i + 1] - 2u[i] + u[i - 1]) / dx^2
    end
    return nothing
end
function heat_dae(n)
    dx = 1 / (n - 1)
    u0 = [sin(π * (i - 1) * dx) for i in 1:n]
    du0 = [0.0; [(u0[i + 1] - 2u0[i] + u0[i - 1]) / dx^2 for i in 2:(n - 1)]; 0.0]
    return DAEProblem{true, SciMLBase.AutoDespecialize}(
        resid!, du0, u0, (0.0, 0.1), (; α = 1.0);
        differential_vars = [false; trues(n - 2); false]
    )
end

@testset "wrapped DAEFunction is differentiated with chunk size 1" begin
    for n in (8, 12)
        prob = heat_dae(n)
        cprob = DiffEqBase.get_concrete_problem(prob, true; alg = DFBDF())
        @test cprob.f.f isa DiffEqBase.FunctionWrappersWrappers.FunctionWrappersWrapper
        ad = prepare_ADType(AutoForwardDiff(), cprob, cprob.u0, cprob.p, true)
        @test ad isa AutoForwardDiff{1}
    end
end

@testset "DFBDF integrator type is independent of the state length" begin
    integs = map((8, 12)) do n
        init(heat_dae(n), DFBDF(); initializealg = BrownFullBasicInit())
    end
    @test typeof(integs[1]) === typeof(integs[2])
    @test typeof(integs[1].cache) === typeof(integs[2].cache)
    for integ in integs
        sol = solve!(integ)
        @test SciMLBase.successful_retcode(sol)
    end
end
