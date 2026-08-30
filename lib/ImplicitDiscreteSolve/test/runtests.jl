#runtests
using Test
using ImplicitDiscreteSolve
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using NonlinearSolve
using SimpleNonlinearSolve
using SciMLBase
using SafeTestsets
using SciMLTesting

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run functional tests
if TEST_GROUP != "QA"
    # Test implicit Euler using ImplicitDiscreteProblem
    @testset "Implicit Euler" begin
        function lotkavolterra(u, p, t)
            [1.5 * u[1] - u[1] * u[2], -3.0 * u[2] + u[1] * u[2]]
        end

        function f!(resid, u_next, u, p, t)
            lv = lotkavolterra(u_next, p, t)
            resid[1] = u_next[1] - u[1] - 0.01 * lv[1]
            resid[2] = u_next[2] - u[2] - 0.01 * lv[2]
            nothing
        end
        u0 = [1.0, 1.0]
        tspan = (0.0, 0.5)

        idprob = ImplicitDiscreteProblem(f!, u0, tspan, []; dt = 0.01)
        idsol = solve(idprob, IDSolve(); adaptive = false)

        oprob = ODEProblem(lotkavolterra, u0, tspan)
        osol = solve(oprob, ImplicitEuler())

        @test isapprox(idsol.u[end], osol.u[end], atol = 0.1)

        ### free-fall
        # y, dy
        function ff(u, p, t)
            [u[2], -9.8]
        end

        function g!(resid, u_next, u, p, t)
            f = ff(u_next, p, t)
            resid[1] = u_next[1] - u[1] - 0.01 * f[1]
            resid[2] = u_next[2] - u[2] - 0.01 * f[2]
            nothing
        end
        u0 = [10.0, 0.0]
        tspan = (0, 0.2)

        idprob = ImplicitDiscreteProblem(g!, u0, tspan, []; dt = 0.01)
        idsol = solve(idprob, IDSolve(); adaptive = false)

        oprob = ODEProblem(ff, u0, tspan)
        osol = solve(oprob, ImplicitEuler())

        @test isapprox(idsol.u[end], osol.u[end], atol = 0.1)
    end

    @testset "Solver initializes" begin
        function periodic!(resid, u_next, u, p, t)
            resid[1] = u_next[1] - u[1] - sin(t * π / 4)
            resid[2] = 16 - u_next[2]^2 - u_next[1]^2
        end

        tsteps = 15
        u0 = [1.0, 3.0]
        idprob = ImplicitDiscreteProblem(periodic!, u0, (0, tsteps), [])
        integ = init(idprob, IDSolve())
        @test integ.u[1]^2 + integ.u[2]^2 ≈ 16

        for ts in 1:tsteps
            step!(integ)
            @test integ.u[1]^2 + integ.u[2]^2 ≈ 16
        end
    end

    @testset "Reinitialize integrator" begin
        function recurrence!(resid, u_next, u, p, t)
            resid[1] = u_next[1] - u[1] - 1
            return nothing
        end

        idprob = ImplicitDiscreteProblem(
            recurrence!, [0.0], (0.0, 2.0), nothing; dt = 1.0
        )
        integ = init(idprob, IDSolve(); adaptive = false)
        nlcache = integ.cache.nlcache
        fresh = solve!(integ)
        fresh_t = copy(fresh.t)
        fresh_u = deepcopy(fresh.u)

        reinit!(integ, [0.0]; t0 = 0.0, tf = 2.0)
        @test integ.u == [1.0]
        @test integ.cache.nlcache === nlcache
        reused = solve!(integ)
        @test reused.t == fresh_t
        @test reused.u == fresh_u

        integ = init(
            idprob, IDSolve(); adaptive = false, save_start = false
        )
        fresh = solve!(integ)
        fresh_t = copy(fresh.t)
        fresh_u = deepcopy(fresh.u)
        reinit!(integ, [0.0]; t0 = 0.0, tf = 2.0)
        @test integ.u == [1.0]
        reused = solve!(integ)
        @test reused.t == fresh_t
        @test reused.u == fresh_u

        function redundant_recurrence!(resid, u_next, u, p, t)
            step_residual = u_next[1] - u[1] - 1
            resid[1] = step_residual
            resid[2] = 2step_residual
            return nothing
        end

        nllsprob = ImplicitDiscreteProblem(
            ImplicitDiscreteFunction(
                redundant_recurrence!; resid_prototype = zeros(2)
            ),
            [0.0], (0.0, 2.0), nothing; dt = 1.0
        )
        integ = init(nllsprob, IDSolve(); adaptive = false)
        nlcache = integ.cache.nlcache
        fresh = solve!(integ)
        fresh_t = copy(fresh.t)
        fresh_u = deepcopy(fresh.u)
        reinit!(integ, [0.0]; t0 = 0.0, tf = 2.0)
        @test integ.u == [1.0]
        @test integ.cache.nlcache === nlcache
        reused = solve!(integ)
        @test reused.t == fresh_t
        @test reused.u == fresh_u

        inconsistent = [false]
        function conditional_recurrence!(resid, u_next, u, p, t)
            step_residual = u_next[1] - u[1] - 1
            resid[1] = step_residual
            resid[2] = p[1] ? step_residual + 1 : 2step_residual
            return nothing
        end

        failureprob = ImplicitDiscreteProblem(
            ImplicitDiscreteFunction(
                conditional_recurrence!; resid_prototype = zeros(2)
            ),
            [0.0], (0.0, 1.0), inconsistent; dt = 1.0
        )
        integ = init(failureprob, IDSolve(); adaptive = false)
        @test check_error(integ) == ReturnCode.Success
        inconsistent[1] = true
        reinit!(integ, [0.0]; t0 = 0.0, tf = 1.0)
        @test check_error(integ) == ReturnCode.InitialFailure
        inconsistent[1] = false
        reinit!(integ, [0.0]; t0 = 0.0, tf = 1.0)
        @test check_error(integ) == ReturnCode.Success
        @test integ.u == [1.0]
    end

    @testset "Constant extrapolant resets nonlinear iterate" begin
        first_guess = Ref(NaN)
        recording = Ref(false)
        function recording_recurrence!(resid, u_next, u, p, t)
            if recording[] && eltype(u_next) === Float64 && isnan(first_guess[])
                first_guess[] = u_next[1]
            end
            resid[1] = u_next[1] - u[1] - 1
            return nothing
        end

        idprob = ImplicitDiscreteProblem(
            recording_recurrence!, [0.0], (0.0, 2.0), nothing; dt = 1.0
        )
        integ = init(idprob, IDSolve(); adaptive = false)
        accepted_u = copy(integ.u)
        SciMLBase.reinit!(integ.cache.nlcache, [100.0])
        recording[] = true
        step!(integ)
        @test first_guess[] == only(accepted_u)
        @test integ.u == [2.0]
    end

    if VERSION >= v"1.11"
        @testset "Cached operations do not allocate" begin
            function nonlinear_recurrence!(resid, u_next, u, p, t)
                resid[1] = u_next[1]^2 - u[1]^2 - 1
                return nothing
            end
            function step_allocations!(integ)
                return @allocated step!(integ)
            end
            function reinit_allocations!(integ, u0)
                return @allocated reinit!(integ, u0; t0 = 0.0, tf = 100.0)
            end

            u0 = [1.0]
            idprob = ImplicitDiscreteProblem(
                nonlinear_recurrence!, u0, (0.0, 100.0), nothing; dt = 1.0
            )
            integ = init(idprob, IDSolve(); adaptive = false, save_on = false)

            step_allocations!(integ)
            @test step_allocations!(integ) == 0
            reinit_allocations!(integ, u0)
            @test reinit_allocations!(integ, u0) == 0
        end
    end

    @testset "Hard problem" begin
        function hard!(resid, u, u_prev, p, t)
            resid[1] = tanh((u[1] - 10t)^2) / 2
        end

        u0 = [0.0]
        idprob = ImplicitDiscreteProblem(hard!, u0, (0.0, 1.0), [])
        integrator = init(idprob, IDSolve())
        idsol = solve!(integrator)
        @test idsol.retcode == ReturnCode.Success
    end

    @testset "Kantorovich strict acceptance" begin
        strict = ImplicitDiscreteSolve.KantorovichTypeController(;
            Θmin = 1 // 8, p = 1
        )
        @test ImplicitDiscreteSolve._accept_kantorovich_step(strict, Float64[])
        @test ImplicitDiscreteSolve._accept_kantorovich_step(strict, [0.2, 0.95])
        @test !ImplicitDiscreteSolve._accept_kantorovich_step(strict, [0.2, 0.951])

        nonstrict = ImplicitDiscreteSolve.KantorovichTypeController(;
            Θmin = 1 // 8, p = 1, strict = false
        )
        @test ImplicitDiscreteSolve._accept_kantorovich_step(nonstrict, [Inf])
    end

    @testset "Handle nothing in u0" begin
        emptyiip(residual, u_next, u, p, t) = nothing
        emptyoop(u_next, u, p, t) = nothing

        tsteps = 5
        u0 = nothing

        idprob = ImplicitDiscreteProblem(emptyiip, u0, (0, tsteps), [])
        sol = solve(idprob, IDSolve())
        @test sol.retcode == ReturnCode.Success

        idprob2 = ImplicitDiscreteProblem(emptyoop, u0, (0, tsteps), [])
        sol2 = solve(idprob2, IDSolve())
        @test sol2.retcode == ReturnCode.Success
    end

    @testset "Create NonlinearLeastSquaresProblem" begin
        function over(u_next, u, p, t)
            [u_next[1] - 1, u_next[2] - 1, u_next[1] - u_next[2]]
        end

        tsteps = 5
        u0 = [1.0, 1.0]
        idprob = ImplicitDiscreteProblem(
            ImplicitDiscreteFunction(over, resid_prototype = zeros(3)), u0, (0, tsteps), []
        )
        integ = init(idprob, IDSolve())
        @test integ.cache.nlcache.prob isa NonlinearLeastSquaresProblem

        function under(u_next, u, p, t)
            [u_next[1] - u_next[2] - 1]
        end
        idprob = ImplicitDiscreteProblem(
            ImplicitDiscreteFunction(under; resid_prototype = zeros(1)), u0, (0, tsteps), []
        )
        integ = init(idprob, IDSolve())
        @test integ.cache.nlcache.prob isa NonlinearLeastSquaresProblem

        function full(u_next, u, p, t)
            [u_next[1]^2 - 3, u_next[2] - u[1]]
        end
        idprob = ImplicitDiscreteProblem(
            ImplicitDiscreteFunction(full; resid_prototype = zeros(2)), u0, (0, tsteps), []
        )
        integ = init(idprob, IDSolve())
        @test integ.cache.nlcache.prob isa NonlinearProblem
    end

    @testset "Null u0 bypasses NonlinearSolve" begin
        # IDSolve overrides allows_null_u0 so `nothing` reaches alg_cache unchanged
        # and the `::Nothing` dispatch returns a cache with nlcache=nothing. Without
        # the opt-in, u0 is coerced to Float64[] and NonlinearSolve fails building
        # a 0x1 Jacobian for the least-squares branch.
        fn = ImplicitDiscreteFunction(; resid_prototype = [0.0]) do du, unext, u, p, t
            du[1] = 1
        end
        pr = ImplicitDiscreteProblem(fn, nothing, (0, 1))
        integ = init(pr, IDSolve())
        @test integ.cache.u === nothing
        @test integ.cache.nlcache === nothing
        sol = solve(pr, IDSolve())
        @test SciMLBase.successful_retcode(sol)
    end

    @testset "nlsolve accepts polyalgorithms and SimpleNonlinearSolve" begin
        function halving!(resid, u_next, u, p, t)
            resid[1] = u_next[1] - u[1] / 2
            return nothing
        end
        idprob = ImplicitDiscreteProblem(halving!, [1.0], (0, 5), nothing)

        reference = solve(idprob, IDSolve())
        @test SciMLBase.successful_retcode(reference)

        for nlsolve in (
                RobustMultiNewton(), FastShortcutNonlinearPolyalg(),
                SimpleNewtonRaphson(), SimpleTrustRegion(),
            )
            sol = solve(idprob, IDSolve(nlsolve = nlsolve))
            @test SciMLBase.successful_retcode(sol)
            @test sol.t == reference.t
            @test sol.u[end] ≈ reference.u[end]
        end
    end

    @testset "InitialFailure thrown" begin
        function bad(u_next, u, p, t)
            [u_next[1] - u_next[2], u_next[1] - 3, u_next[2] - 4]
        end

        u0 = [3.0, 4.0]
        idprob = ImplicitDiscreteProblem(bad, u0, (0, 0), [])
        integ = init(idprob, IDSolve())
        @test check_error(integ) == ReturnCode.InitialFailure
        sol = solve(idprob, IDSolve())
        @test length(sol.u) == 1
        @test !SciMLBase.successful_retcode(sol)
    end
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "Core" && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
