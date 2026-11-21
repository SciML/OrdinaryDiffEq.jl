using OrdinaryDiffEqBDF, LinearAlgebra, Test
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, ShampineCollocationInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using SciMLBase: AutoDespecialize, DAEFunction, DAEProblem, successful_retcode
import DifferentiationInterface as DI

afd_cs3 = AutoForwardDiff(chunksize = 3)

struct DespecializedDAEResidual{ID} end

function (::DespecializedDAEResidual)(resid, du, u, p, t)
    resid[1] = du[1] + p.rate * u[1]
    return nothing
end

function despecialized_dae_problem(model, p)
    f = DAEFunction{true, AutoDespecialize}(model)
    return DAEProblem(f, [-p.rate], [1.0], (0.0, 0.1), p)
end

@testset "AutoDespecialize reuses DAE solver caches" begin
    problems = (
        despecialized_dae_problem(
            DespecializedDAEResidual{1}(), (rate = 0.5,)
        ),
        despecialized_dae_problem(
            DespecializedDAEResidual{2}(), (rate = 0.5, unused = 1)
        ),
    )

    for Alg in (DFBDF, DNordsieckBDF)
        for autodiff in (AutoForwardDiff(), AutoFiniteDiff())
            alg = Alg(; autodiff)
            integrators = map(prob -> init(prob, alg), problems)
            @test typeof(integrators[1].sol.prob) === typeof(integrators[2].sol.prob)
            @test typeof(integrators[1].cache) === typeof(integrators[2].cache)

            for prob in problems
                sol = solve(prob, alg)
                @test successful_retcode(sol)
                @test sol.u[end][1] ≈ exp(-0.05) rtol = 1.0e-3
            end
        end
    end
end

function f(out, du, u, p, t)
    out[1] = -p[1] * u[1] + p[3] * u[2] * u[3] - du[1]
    out[2] = +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
    return
end
function f(du, u, p, t)
    return [
        -p[1] * u[1] + p[3] * u[2] * u[3] - du[1],
        +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3] - du[2],
        u[1] + u[2] + u[3] - 1.0,
    ]
end
function f_ode(du, u, p, t)
    return du .= [
        -p[1] * u[1] + p[3] * u[2] * u[3],
        +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3],
        u[1] + u[2] + u[3] - 1.0,
    ]
end
function f_ode(u, p, t)
    return [
        -p[1] * u[1] + p[3] * u[2] * u[3],
        +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3],
        u[1] + u[2] + u[3] - 1.0,
    ]
end
p = [0.04, 3.0e7, 1.0e4]
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
M = Diagonal([1.0, 1.0, 0.0])
prob = DAEProblem(f, du₀, u₀, tspan, p; differential_vars)
prob_oop = DAEProblem{false}(f, du₀, u₀, tspan, p; differential_vars)
f_mm = ODEFunction{true}(f_ode, mass_matrix = M)
prob_mm = ODEProblem(f_mm, u₀, tspan, p)
f_mm_oop = ODEFunction{false}(f_ode, mass_matrix = M)
prob_mm_oop = ODEProblem(f_mm_oop, u₀, tspan, p)
if VERSION >= v"1.12"
    sol1 = @inferred solve(
        prob, DFBDF(autodiff = afd_cs3), dt = 1.0e-5, abstol = 1.0e-8, reltol = 1.0e-8
    )
    sol2 = @inferred solve(
        prob_oop, DFBDF(autodiff = afd_cs3), dt = 1.0e-5, abstol = 1.0e-8, reltol = 1.0e-8
    )
    sol3 = @inferred solve(
        prob_mm, FBDF(autodiff = afd_cs3), dt = 1.0e-5, abstol = 1.0e-8, reltol = 1.0e-8
    )
end

# These tests flex differentiation of the solver and through the initialization
# To only test the solver part and isolate potential issues, set the initialization to consistent
@testset "Inplace: $(isinplace(_prob)), DAEProblem: $(_prob isa DAEProblem), BrownBasic: $(initalg isa BrownFullBasicInit), Autodiff: $autodiff" for _prob in [
            prob, prob_oop, prob_mm, prob_mm_oop,
        ],
        initalg in [BrownFullBasicInit(), ShampineCollocationInit()],
        autodiff in [afd_cs3, AutoFiniteDiff()]

    alg = (_prob isa DAEProblem) ? DFBDF(; autodiff) : FBDF(; autodiff)
    function f_loss(p)
        sol = solve(
            remake(_prob; p), alg, abstol = 1.0e-14,
            reltol = 1.0e-14, initializealg = initalg
        )
        sum(sum, sol.u)
    end
    @test DI.gradient(f_loss, AutoForwardDiff(), [0.04, 3.0e7, 1.0e4]) ≈ [0, 0, 0] atol = 1.0e-8
end
