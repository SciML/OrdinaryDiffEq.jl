using OrdinaryDiffEqSDIRK
using Test

vdp!(du, u, p, t) = (du[1] = u[2]; du[2] = 1.0e3 * (1 - u[1]^2) * u[2] - u[1]; nothing)
prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3))
ref = solve(prob, ESDIRK659L2SA(predictor = Predictor.Trivial), reltol = 1.0e-12, abstol = 1.0e-12)

methods = [
    ESDIRK325L2SA, ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2,
    ESDIRK659L2SA,
]
predictors = [
    Predictor.Trivial, Predictor.Linear, Predictor.CopyPrev, Predictor.StageExtrap,
    Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder, Predictor.Tableau,
]

@testset "Stage predictors" begin
    for M in methods
        @testset "$(nameof(M))" begin
            triv = solve(prob, M(predictor = Predictor.Trivial), reltol = 1.0e-8, abstol = 1.0e-8)
            for p in predictors
                s = solve(prob, M(predictor = p), reltol = 1.0e-8, abstol = 1.0e-8)
                @test all(isfinite, s.u[end])
                @test maximum(abs.(s.u[end] .- ref.u[end])) < 1.0e-5
            end
            se = solve(prob, M(predictor = Predictor.StageExtrap), reltol = 1.0e-8, abstol = 1.0e-8)
            @test se.stats.nf <= triv.stats.nf
        end
    end
end

@testset "ImplicitEuler first-stage predictor seed" begin
    # The direct ImplicitEuler solve reaches the implicit-first-stage branch of
    # `generic_imex_perform_step.jl`. The predictor selects the nlsolver's stage-1
    # initial guess, so recording the argument of every `f` call pins the seed:
    # `Trivial` guesses z = 0 (first residual at uprev), `Linear` guesses
    # z = dt*f(uprev) (first residual at uprev + z). The fsal evaluation at uprev
    # is always the first recorded call.
    for isinplace in (false, true)
        calls = Float64[]
        if isinplace
            function f!(du, u, p, t)
                eltype(u) === Float64 && push!(calls, u[1])
                du[1] = 2u[1] - u[1]^2 + 2
                return nothing
            end
            local_prob = ODEProblem(f!, [0.0], (0.0, 1.0))
        else
            function f(u, p, t)
                u isa Float64 && push!(calls, u)
                return 2u - u^2 + 2
            end
            local_prob = ODEProblem(f, 0.0, (0.0, 1.0))
        end

        solve(
            local_prob, ImplicitEuler(predictor = Predictor.Trivial);
            adaptive = false, dt = 1.0
        )
        @test calls[1:2] == [0.0, 0.0]

        empty!(calls)
        solve(
            local_prob, ImplicitEuler(predictor = Predictor.Linear);
            adaptive = false, dt = 1.0
        )
        @test calls[1:2] == [0.0, 2.0]
    end
end

@testset "ImplicitEuler MaxOrder stage-1 seed" begin
    # MaxOrder seeds stage 1 from the previous step's interpolant (from step 2 on).
    # Out-of-place scalar `u` is immutable, so the seed must be built without the
    # in-place `current_extrapolant!`. The predictor only changes the nlsolver's
    # initial guess, so the converged trajectory must match `Trivial`.
    for u0 in (1.0, [1.0])  # scalar out-of-place, then in-place
        prob = u0 isa Number ?
            ODEProblem((u, p, t) -> -u, u0, (0.0, 1.0)) :
            ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), u0, (0.0, 1.0))
        sol_max = solve(
            prob, ImplicitEuler(predictor = Predictor.MaxOrder);
            adaptive = false, dt = 0.01
        )
        sol_triv = solve(
            prob, ImplicitEuler(predictor = Predictor.Trivial);
            adaptive = false, dt = 0.01
        )
        @test all(isfinite, Array(sol_max))
        @test Array(sol_max) ≈ Array(sol_triv) rtol = 1.0e-6
        @test Array(sol_max)[end] ≈ exp(-1.0) rtol = 0.02
    end
end

# deprecated `extrapolant` keyword still maps onto a predictor
@testset "extrapolant deprecation" begin
    @test ImplicitEuler(extrapolant = :linear).predictor == Predictor.Linear
    @test ImplicitEuler(extrapolant = :constant).predictor == Predictor.Trivial
    @test KenCarp4(extrapolant = :interpolant).predictor == Predictor.MaxOrder
end
