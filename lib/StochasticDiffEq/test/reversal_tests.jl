using StochasticDiffEq, DiffEqNoiseProcess, Test, Random
using SDEProblemLibrary
# automatically construct SDE transformation for Ito reversal
using ModelingToolkit
import SciMLBase

# tested solvers
additive_noise_solver = [
    SRA(),
    SRA1(),
    SRA2(),
    SRA3(),
    SOSRA(),
    SOSRA2(),
    SKenCarp(),
]

Stratonovich_solver = [
    # non-stiff methods
    EulerHeun(),
    LambaEulerHeun(),
    RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    RKMilCommute(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    # S-Rock methods
    SROCK1(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    # stiff methods
    ImplicitEulerHeun(),
    ImplicitRKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    ISSEulerHeun(),
]

Ito_solver = [
    # non-stiff methods
    EM(),
    LambaEM(),
    WangLi3SMil_A(),
    WangLi3SMil_B(),
    WangLi3SMil_C(),
    WangLi3SMil_D(),
    WangLi3SMil_E(),
    WangLi3SMil_F(),
    RKMil(),
    SRI(),
    SRIW1(),
    SRIW2(),
    SOSRI(),
    SOSRI2(),
    # S-Rock methods
    SROCK1(),
    SROCKEM(),
    SROCK2(),
    SKSROCK(),
    TangXiaoSROCK2(),
    # stiff methods
    ImplicitEM(),
    ImplicitRKMil(),
    ISSEM(),
]

seed = 62356236796
Random.seed!(seed)

dt = 1.0e-6
n = round(Int, 1 / dt) # number of steps, probs have time span (0,1)
W = [0.0; cumsum(sqrt(dt) * randn(n + 1))] # n + 1 to avoid past step in noise grid
Z = [0.0; cumsum(sqrt(dt) * randn(n + 1))]
ts = collect(0:dt:(1 + dt))

W_forward = NoiseGrid(ts, W, Z)
W_reverse = reverse(W_forward)

@testset "Additive Noise Solver Reversal Tests ($(["out-of-place", "in-place"][i]))" for i in
    1:2

    prob = (
        SDEProblemLibrary.prob_sde_additive,
        SDEProblemLibrary.prob_sde_additivesystem,
    )[i]

    prob_forward = remake(prob, noise = W_forward)

    for solver in additive_noise_solver
        println("solver: ", solver)
        sol_forward = solve(prob_forward, solver, dt = dt, adaptive = false)

        prob_reverse = remake(
            prob_forward, noise = W_reverse, tspan = reverse(prob.tspan), u0 = sol_forward.u[end]
        )
        sol_reverse = solve(prob_reverse, solver, dt = dt, adaptive = false)

        @test sol_forward(ts).u ≈ sol_reverse(ts).u rtol = 1.0e-3
        @test length(sol_forward.t) == length(sol_reverse.t)
        GC.gc()
    end
end

@testset "Stratonovich Solver Reversal Tests ($(["out-of-place", "in-place"][i]))" for i in
    1:2

    prob = (
        SDEProblemLibrary.prob_sde_linear_stratonovich,
        SDEProblemLibrary.prob_sde_2Dlinear_stratonovich,
    )[i]

    prob_forward = remake(prob, noise = W_forward)

    for solver in Stratonovich_solver
        println("solver: ", solver)
        sol_forward = solve(prob_forward, solver, dt = dt, adaptive = false)

        prob_reverse = remake(
            prob_forward, noise = W_reverse, tspan = reverse(prob.tspan), u0 = sol_forward.u[end]
        )
        sol_reverse = solve(prob_reverse, solver, dt = dt, adaptive = false)

        @test sol_forward(ts).u ≈ sol_reverse(ts).u rtol = 1.0e-2
        @test length(sol_forward.t) == length(sol_reverse.t)
        GC.gc()
    end
end

@testset "Ito Solver Reversal Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        SDEProblemLibrary.prob_sde_linear,
        SDEProblemLibrary.prob_sde_2Dlinear,
    )[i]
    if i == 1
        prob_forward = remake(prob, noise = W_forward)
    else
        prob_forward = remake(prob, noise = W_forward, u0 = vec(prob.u0))
    end
    sys = complete(modelingtoolkitize(prob_forward))
    mtkps = MTKParameters(sys, [])
    sys2 = stochastic_integral_transform(sys, -1 // 1)
    fdrift = generate_rhs(sys2; expression = Val{false})[i]
    fdif = generate_diffusion_function(sys2; expression = Val{false})[i]

    for solver in Ito_solver
        println("solver: ", solver)
        sol_forward = solve(prob_forward, solver, dt = dt, adaptive = false)
        if i == 1
            _u0 = [sol_forward.u[end]]
        else
            _u0 = sol_forward.u[end]
        end
        prob_reverse = remake(
            prob_forward, f = SDEFunction(fdrift, fdif), noise = W_reverse,
            tspan = reverse(prob.tspan), u0 = _u0, p = mtkps
        )
        sol_reverse = solve(prob_reverse, solver, dt = dt, adaptive = false)

        if i == 1
            @test sol_forward(ts).u ≈ vcat(sol_reverse(ts).u...) rtol = 1.0e-2
        else
            @test sol_forward(ts).u ≈ sol_reverse(ts).u rtol = 1.0e-2
        end
        @test length(sol_forward.t) == length(sol_reverse.t)
        GC.gc()
    end
end
