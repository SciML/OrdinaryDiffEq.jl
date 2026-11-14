using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqSDIRK, DiffEqBase, LineSearches, LinearAlgebra, Test
using OrdinaryDiffEqBDF: DImplicitEuler, QNDF, QBDF

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

@testset "Coefficient-multistep correction tolerance" begin
    magnitude = 1.0e16
    f!(du, u, p, t) = (du[1] = 5 * (u[1] - magnitude) + magnitude + 2)
    jac!(J, u, p, t) = (J[1, 1] = 5)
    predictor!(u, uprev, p, t, dt) = fill!(u, magnitude)
    prob = ODEProblem(ODEFunction(f!; jac = jac!), [1.0], (0.0, 1.0))
    # The correction is below one ULP of the predictor and within its relative tolerance.
    for relax in (nothing, BackTracking(), LineSearches.Static()), always_new in (false, true)
        alg = QBDF(nlsolve = NLNewton(; relax, predictor = predictor!, always_new))
        sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 1.0e-6)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end] == [magnitude]
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

@testset "Strictly decreasing stalled Newton line search" begin
    h(x) = x <= 0.25 ? 0.0 : x >= 0.75 ? 1.0 : (y = 2 * (x - 0.25); y * y * (3 - 2y))
    function rhs!(du, u, p, t)
        du[1] = -2
        du[2] = -1 + 1.0e6 * (u[1] - h(u[2]) * (u[1] - 1))
        return du
    end
    rhs(u, p, t) = rhs!(similar(u), u, p, t)
    predictor!(u, uprev, p, t, dt) = fill!(u, 0)
    predictor(uprev, p, t, dt) = zero(uprev)
    for iip in (true, false), always_new in (false, true)
        @testset "iip=$iip fresh=$always_new" begin
            accepted = Tuple{Float64, Float64, Float64}[]
            function linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
                α, merit = BackTracking()(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
                push!(accepted, (α, ϕ0, merit))
                return α, merit
            end
            prob = ODEProblem{iip}(iip ? rhs! : rhs, ones(2), (0.0, 1.0))
            alg = ImplicitEuler(
                nlsolve = NLNewton(;
                    relax = linesearch, predictor = iip ? predictor! : predictor,
                    always_new, max_iter = 50
                )
            )
            sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 0.0)
            expected = always_new ? SciMLBase.ReturnCode.Success : SciMLBase.ReturnCode.ConvergenceFailure
            @test sol.retcode == expected
            if always_new
                @test sol.u[end][1] ≈ -1 atol = 1.0e-3 rtol = 0
                @test norm(sol.u[end] - ones(2) - rhs(sol.u[end], nothing, 1.0), Inf) <= 1.0e-3
            else
                α, initial, merit = first(accepted)
                @test α < 1.0e-10
                @test merit < initial
            end
        end
    end
end

@testset "Multidimensional in-place line search" begin
    rhs!(du, u, p, t) = (du .= -u; du)
    rhs(u, p, t) = -u
    for dims in ((2, 2), (2, 1, 2)), iip in (true, false)
        u0 = reshape(collect(1.0:prod(dims)), dims)
        prob = ODEProblem{iip}(iip ? rhs! : rhs, u0, (0.0, 0.1))
        sol = solve(
            prob, ImplicitEuler(nlsolve = NLNewton(relax = BackTracking()));
            dt = 0.1, adaptive = false
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end] ≈ u0 ./ 1.1 rtol = 1.0e-12
    end
end

@testset "Damped Newton progress with growing full corrections" begin
    Fmono(x) = x + 0.99 * sin(x) - 2π
    rhs!(du, u, p, t) = (du[1] = u[1] - 1 - Fmono(u[1]); du)
    rhs(u, p, t) = rhs!(similar(u), u, p, t)
    predictor!(u, uprev, p, t, dt) = fill!(u, 0.0)
    predictor(uprev, p, t, dt) = zero(uprev)
    for iip in (true, false)
        accepted = Tuple{Float64, Float64, Float64}[]
        function linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
            α, merit = BackTracking()(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
            push!(accepted, (α, ϕ0, merit))
            return α, merit
        end
        prob = ODEProblem{iip}(iip ? rhs! : rhs, [1.0], (0.0, 1.0))
        alg = ImplicitEuler(
            nlsolve = NLNewton(;
                relax = linesearch, predictor = iip ? predictor! : predictor,
                always_new = true, max_iter = 50
            )
        )
        sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-9, reltol = 0.0)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ 2π rtol = 1.0e-12
        @test length(accepted) >= 3
        @test accepted[2][1] ≈ 0.01
        @test accepted[2][3] < accepted[2][2]
    end
end

@testset "Float32 StrongWolfe line search" begin
    for g in (u -> -u, u -> -u .^ 3), iip in (true, false)
        rhs!(du, u, p, t) = (du .= g(u); du)
        rhs(u, p, t) = g(u)
        prob = ODEProblem{iip}(iip ? rhs! : rhs, Float32[1, 2], (0.0f0, 1.0f0))
        alg = ImplicitEuler(nlsolve = NLNewton(relax = StrongWolfe()))
        sol = solve(prob, alg; dt = 0.1f0, adaptive = false)
        plain = solve(
            prob, ImplicitEuler(nlsolve = NLNewton()); dt = 0.1f0, adaptive = false
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test eltype(sol.u[end]) == Float32
        @test maximum(abs.(sol.u[end] .- plain.u[end])) <= 1.0f-3
    end
end

@testset "Coefficient-multistep line-search weights" begin
    rhs!(du, u, p, t) = (du .= -u; du)
    rhs(u, p, t) = -u
    predictor!(u, uprev, p, t, dt) = (u .= (1000.0, 0.0); u)
    predictor(uprev, p, t, dt) = [1000.0, 0.0]
    first_merits = Float64[]
    for iip in (true, false)
        first_merit = Ref{Union{Nothing, Float64}}(nothing)
        function linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
            first_merit[] === nothing && (first_merit[] = ϕ0)
            return BackTracking()(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        end
        prob = ODEProblem{iip}(iip ? rhs! : rhs, [1.0, 1.0], (0.0, 1.0))
        alg = QBDF(
            nlsolve = NLNewton(;
                relax = linesearch, predictor = iip ? predictor! : predictor
            )
        )
        sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-6, reltol = 1.0e-3)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end] ≈ [0.5, 0.5]
        push!(first_merits, something(first_merit[]))
    end
    @test first_merits[1] ≈ first_merits[2] rtol = 1.0e-12
    @test first_merits[1] ≈ 1580.1892298099906 rtol = 1.0e-12
end

@testset "Tiny undamped Newton corrections" begin
    rhs!(du, u, p, t) = (du[1] = 2 * (u[1] - 1) + p; du)
    rhs(u, p, t) = rhs!(similar(u), u, p, t)
    predictor!(u, uprev, p, t, dt) = (u .= prevfloat.(uprev))
    predictor(uprev, p, t, dt) = prevfloat.(uprev)
    for iip in (true, false), c in (0.0, 1.0e-10), always_new in (false, true)
        prob = ODEProblem{iip}(iip ? rhs! : rhs, [1.0], (0.0, 1.0), c)
        alg = ImplicitEuler(
            nlsolve = NLNewton(;
                relax = BackTracking(), predictor = iip ? predictor! : predictor, always_new
            )
        )
        sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-12, reltol = 0.0)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ 1 - c atol = 1.0e-12 rtol = 0
    end
end
