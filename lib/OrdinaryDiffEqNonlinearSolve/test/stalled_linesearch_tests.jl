using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqSDIRK, DiffEqBase, LineSearches, LinearAlgebra, Test
using OrdinaryDiffEqBDF: DImplicitEuler, QNDF

@testset "Stalled Newton line search" begin
    ϵ = 1.0e-5
    Bmat = [1.0 2 / ϵ; 0.0 1.0]
    Amat = Bmat * [1.0 0.0; -ϵ 1 - (1 - ϵ)]
    cvec = [1.0, 0.0]
    hswitch(x) = x <= 0.25 ? 0.0 : x >= 0.75 ? 1.0 : (y = 2 * (x - 0.25); y * y * (3 - 2y))
    Fstage(u) = Amat * u + cvec + hswitch(u[2]) * (Bmat - Amat) * (u .- 1)

    # The stale Jacobian at uprev gives an ascent direction at the zero predictor.
    for magnitude in (0.0, 1.0e7, 1.0e9, 1.0e12), iip in (true, false),
            always_new in (false, true)
        function rhs!(du, u, p, t)
            du[1:2] .= u[1:2] .- 1 .- Fstage(u[1:2])
            if !iszero(p)
                du[3] = p
            end
            return du
        end
        rhs(u, p, t) = rhs!(similar(u), u, p, t)
        function predictor!(u, uprev, p, t, dt)
            fill!(u, 0.0)
            if !iszero(p)
                u[3] = uprev[3] + p * dt
            end
            return u
        end
        predictor(uprev, p, t, dt) = predictor!(similar(uprev), uprev, p, t, dt)
        u0 = iszero(magnitude) ? ones(2) : [1.0, 1.0, 0.0]
        prob = ODEProblem{iip}(iip ? rhs! : rhs, u0, (0.0, 1.0), magnitude)
        alg = ImplicitEuler(
            nlsolve = NLNewton(
                relax = BackTracking(), predictor = iip ? predictor! : predictor,
                always_new = always_new, max_iter = 50
            )
        )
        @testset "padding=$magnitude iip=$iip fresh=$always_new" begin
            sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 0.0)
            expected = always_new ? SciMLBase.ReturnCode.Success : SciMLBase.ReturnCode.ConvergenceFailure
            @test sol.retcode == expected
            @test sol.retcode != SciMLBase.ReturnCode.Success ||
                norm(Fstage(sol.u[end][1:2]), Inf) <= 1.0e-3
        end
    end
end

@testset "Out-of-place line-search residuals" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 0.1))
    for alg in (ImplicitEuler, QNDF)
        reference = solve(prob, alg(); dt = 0.1, adaptive = false)
        sol = solve(
            prob, alg(nlsolve = NLNewton(relax = BackTracking()));
            dt = 0.1, adaptive = false
        )
        @test reference.retcode == sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-12
    end
    dae = DAEProblem{false}((du, u, p, t) -> du .+ u, [-1.0], [1.0], (0.0, 0.1))
    sol = solve(
        dae, DImplicitEuler(nlsolve = NLNewton(relax = BackTracking()));
        dt = 0.1, adaptive = false
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end][1] ≈ 1 / 1.1 rtol = 1.0e-12
end
