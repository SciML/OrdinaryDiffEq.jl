using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using ADTypes, SciMLBase
using Test

# `nlstep_data` is a field of `ODEFunction`/`SplitFunction` only, so `build_nlsolver` and
# the `NonlinearSolveAlg` step have to look it up by field: reading it outright threw on
# every other `AbstractSciMLFunction` before the first step.

@testset "DynamicalODEFunction carries no nlstep_data" begin
    # A `SecondOrderODEProblem` builds a `DynamicalODEFunction`; `u'' = -u` has
    # `u = cos(t)`, `v = -sin(t)`.
    ff(dv, v, u, p, t) = (dv .= -u)
    prob = SecondOrderODEProblem(ff, [0.0], [1.0], (0.0, 1.0))

    alg = ImplicitEuler(
        nlsolve = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoFiniteDiff()))
    )
    sol = solve(prob, alg; reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end].x[1][1] ≈ -sin(1.0) rtol = 1.0e-3
    @test sol.u[end].x[2][1] ≈ cos(1.0) rtol = 1.0e-3
end

# The `nlstep_data` stage systems below are hand-rolled stand-ins for what ModelingToolkit
# generates, so this sublibrary can test the contract without depending on it. Each holds
# its stage coefficients in a mutable parameter object that the six hooks write into.
mutable struct StageParams
    γ₁::Float64
    γ₂::Float64
    γ₃::Float64
    c::Float64
    inner::Vector{Float64}
    outer::Vector{Float64}
end

set_γ_c!(prob, γc) = ((prob.p.γ₁, prob.p.γ₂, prob.p.γ₃, prob.p.c) = γc; nothing)
set_outer!(prob, tmp) = (copyto!(prob.p.outer, tmp); nothing)
set_inner!(prob, tmp) = (copyto!(prob.p.inner, tmp); nothing)

# `0 = du[1] + u[1] - u[2]`, `0 = u[1] + u[2] - 1` is index-1 with the closed form
# `u[1] = 1/2 + exp(-2t)/2`, `u[2] = 1 - u[1]`.
function dae_f!(out, du, u, p, t)
    out[1] = du[1] + u[1] - u[2]
    out[2] = u[1] + u[2] - 1
    return nothing
end

# The fully implicit stage system: `F(γ₁ z + outer_tmp, γ₂ z + inner_tmp, p, c)`.
function dae_stage!(res, z, p)
    du = p.γ₁ .* z .+ p.outer
    u = p.γ₂ .* z .+ p.inner
    dae_f!(res, du, u, nothing, p.c)
    return nothing
end

@testset "DAEProblem with nlstep_data solves the fully implicit stage system" begin
    # `DAEFunction` only gained the field in SciMLBase 3.44.
    if :nlstep_data in fieldnames(DAEFunction)
        nlstep = SciMLBase.ODENLStepData(
            NonlinearProblem(
                NonlinearFunction{true}(dae_stage!), zeros(2),
                StageParams(1.0, 1.0, 1.0, 0.0, zeros(2), zeros(2))
            ),
            1:2, set_γ_c!, set_outer!, set_inner!, (ztmp, sol) -> (ztmp .= sol.u)
        )
        u0, du0, tspan = [1.0, 0.0], [-1.0, 0.0], (0.0, 1.0)
        dvars = [true, false]
        prob = DAEProblem(
            DAEFunction(dae_f!; nlstep_data = nlstep), du0, u0, tspan;
            differential_vars = dvars
        )
        alg = DFBDF(nlsolve = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoFiniteDiff())))
        sol = solve(prob, alg; reltol = 1.0e-10, abstol = 1.0e-10)

        # Before the stage assignment existed, both `initialize!` arms stated the stage
        # system in mass-matrix ODE conventions, and this returned `Success` sitting at the
        # initial condition.
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ 0.5 + exp(-2.0) / 2 rtol = 1.0e-6
        @test sol.u[end][2] ≈ 0.5 - exp(-2.0) / 2 rtol = 1.0e-6
        @test sol.u[end] ≉ u0
    end
end

odef!(du, u, p, t) = (du .= .-u .^ 3 .- u)

# The mass-matrix stage system: `outer_tmp + γ₁ f(γ₂ z + inner_tmp, p, c) - γ₃ M z`.
function ode_stage!(res, z, p)
    u = p.γ₂ .* z .+ p.inner
    k = similar(u)
    odef!(k, u, nothing, p.c)
    @. res = p.outer + p.γ₁ * k - p.γ₃ * z
    return nothing
end

@testset "COEFFICIENT_MULTISTEP nlstep stage residual is on the solution's scale" begin
    # `FBDF` drives the stage system through the `COEFFICIENT_MULTISTEP` arm, whose raw
    # residual carries an `α * invγdt` factor. The `nlstep_data` branch measures that
    # residual for convergence, so leaving it unscaled read as divergence at small `dt`
    # and walked the solver to `dtmin`: `Unstable`, sitting at `u0`.
    u0, tspan = [1.0], (0.0, 1.0)
    ref = solve(ODEProblem(odef!, u0, tspan), FBDF(); reltol = 1.0e-10, abstol = 1.0e-10)

    nlstep = SciMLBase.ODENLStepData(
        NonlinearProblem(
            NonlinearFunction{true}(ode_stage!), zeros(1),
            StageParams(1.0, 1.0, 1.0, 0.0, zeros(1), zeros(1))
        ),
        1:1, set_γ_c!, set_outer!, set_inner!, (ztmp, sol) -> (ztmp .= sol.u)
    )
    prob = ODEProblem(ODEFunction(odef!; nlstep_data = nlstep), u0, tspan)
    alg = FBDF(nlsolve = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoFiniteDiff())))
    sol = solve(prob, alg; reltol = 1.0e-10, abstol = 1.0e-10)

    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-6
    @test sol.u[end] ≉ u0
end
