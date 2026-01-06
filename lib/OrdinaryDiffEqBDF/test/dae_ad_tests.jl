using OrdinaryDiffEqBDF, LinearAlgebra, ForwardDiff, Test
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, ShampineCollocationInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff

afd_cs3 = AutoForwardDiff(chunksize = 3)

function f(out, du, u, p, t)
    out[1] = -p[1] * u[1] + p[3] * u[2] * u[3] - du[1]
    out[2] = +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3] - du[2]
    return out[3] = u[1] + u[2] + u[3] - 1.0
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
prob = DAEProblem(f, du₀, u₀, tspan, p, differential_vars = differential_vars)
prob_oop = DAEProblem{false}(f, du₀, u₀, tspan, p, differential_vars = differential_vars)
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
@testset "Inplace: $(isinplace(_prob)), DAEProblem: $(_prob isa DAEProblem), BrownBasic: $(initalg isa BrownFullBasicInit), Autodiff: $autodiff" for _prob in
        [
            prob, prob_oop, prob_mm, prob_mm_oop,
        ],
        initalg in [BrownFullBasicInit(), ShampineCollocationInit()],
        autodiff in [afd_cs3, AutoFiniteDiff()]

    alg = (_prob isa DAEProblem) ? DFBDF(; autodiff) : FBDF(; autodiff)
    function f(p)
        sol = solve(
            remake(_prob, p = p), alg, abstol = 1.0e-14,
            reltol = 1.0e-14, initializealg = initalg
        )
        sum(sol)
    end
    @test ForwardDiff.gradient(f, [0.04, 3.0e7, 1.0e4]) ≈ [0, 0, 0] atol = 1.0e-8
end
