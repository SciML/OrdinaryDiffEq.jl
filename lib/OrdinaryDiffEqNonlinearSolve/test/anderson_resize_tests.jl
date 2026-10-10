using OrdinaryDiffEqNonlinearSolve: NLAnderson, NLFunctional, NLNewton
using OrdinaryDiffEqSDIRK
using SciMLBase
using Test

fdecay!(du, u, p, t) = (du .= -p[1] .* u; nothing)
feq!(du, u, p, t) = (du .= 1 .- u; nothing)

# Until SciML/OrdinaryDiffEq.jl#4722, callback grow must init new uprev/uprev2 from u.
@testset "NLAnderson resize" begin
    @testset "grow/shrink matches NLNewton" begin
        prob = ODEProblem(feq!, fill(0.2, 3), (0.0, 2.0))

        function make_cb()
            phase = Ref(0)
            condition = function (u, t, integ)
                return (phase[] == 0 && t >= 0.5) || (phase[] == 1 && t >= 1.0)
            end
            affect! = function (integ)
                if phase[] == 0
                    resize!(integ, 5)
                    integ.u[4:5] .= 0.5
                    integ.uprev[4:5] .= integ.u[4:5]
                    integ.uprev2[4:5] .= integ.u[4:5]
                    phase[] = 1
                else
                    resize!(integ, 3)
                    phase[] = 2
                end
                u_modified!(integ, true)
                return nothing
            end
            return DiscreteCallback(condition, affect!)
        end

        kwargs = (;
            tstops = [0.5, 1.0],
            reltol = 1.0e-12, abstol = 1.0e-12,
            dt = 0.025, adaptive = false, saveat = 0.025,
        )
        sol_ref = solve(
            prob, ImplicitEuler(nlsolve = NLNewton()); callback = make_cb(), kwargs...
        )
        sol_and = solve(
            prob, ImplicitEuler(nlsolve = NLAnderson()); callback = make_cb(), kwargs...
        )
        sol_fun = solve(
            prob, ImplicitEuler(nlsolve = NLFunctional()); callback = make_cb(), kwargs...
        )

        @test SciMLBase.successful_retcode(sol_ref)
        @test SciMLBase.successful_retcode(sol_and)
        @test SciMLBase.successful_retcode(sol_fun)
        @test 5 in length.(sol_ref.u)
        @test 5 in length.(sol_and.u)
        @test 5 in length.(sol_fun.u)
        @test length(sol_and.u[end]) == 3
        @test length(sol_fun.u[end]) == 3
        @test sol_and.t == sol_ref.t
        @test sol_fun.t == sol_ref.t
        for i in eachindex(sol_ref.t)
            @test length(sol_and.u[i]) == length(sol_ref.u[i])
            @test length(sol_fun.u[i]) == length(sol_ref.u[i])
            @test sol_and.u[i] ≈ sol_ref.u[i] rtol = 1.0e-8 atol = 1.0e-8
            @test sol_fun.u[i] ≈ sol_ref.u[i] rtol = 1.0e-8 atol = 1.0e-8
        end
    end

    @testset "rebuilds Q/R and Δz₊s when max_history is unchanged" begin
        prob = ODEProblem(fdecay!, ones(10), (0.0, 1.0), 1.0)
        integ = init(prob, TRBDF2(nlsolve = NLAnderson()); dt = 0.1, adaptive = false)
        step!(integ)
        nlc = integ.cache.nlsolver.cache
        @test size(nlc.Q) == (10, 5)
        @test all(v -> length(v) == 10, nlc.Δz₊s)

        resize!(integ, 15)
        nlc = integ.cache.nlsolver.cache
        @test size(nlc.Q) == (15, 5)
        @test size(nlc.R) == (5, 5)
        @test length(nlc.Δz₊s) == 5
        @test all(v -> length(v) == 15, nlc.Δz₊s)
        @test nlc.history == 0
    end

    @testset "callback grow with fixed max_history" begin
        grew = Ref(false)
        cb = DiscreteCallback(
            (u, t, integ) -> !grew[] && t >= 0.5,
            function (integ)
                resize!(integ, 15)
                integ.u[11:15] .= 1.0
                integ.uprev[11:15] .= integ.u[11:15]
                integ.uprev2[11:15] .= integ.u[11:15]
                grew[] = true
                u_modified!(integ, true)
                return nothing
            end,
        )
        sol = solve(
            ODEProblem(fdecay!, ones(10), (0.0, 1.0), 50.0),
            ImplicitEuler(nlsolve = NLAnderson());
            callback = cb, tstops = [0.5], reltol = 1.0e-4, abstol = 1.0e-6,
        )
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.u[end]) == 15
        @test 15 in length.(sol.u)
    end
end
