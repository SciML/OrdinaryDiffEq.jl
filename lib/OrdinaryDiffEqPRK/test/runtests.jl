using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

@time @safetestset "SciMLBase reexport" begin
    using OrdinaryDiffEqPRK, Test
    exported = (
        :ODEProblem, :ODEFunction, :SplitODEProblem, :solve, :init, :step!,
        :remake, :ReturnCode, :CallbackSet, :ContinuousCallback, :terminate!,
        :u_modified!, :add_tstop!, :get_du, :EnsembleProblem,
    )
    @test all(Base.isexported.(Ref(OrdinaryDiffEqPRK), exported))
    internal = (
        :build_solution, :isinplace, :has_jac, :AbstractODEProblem,
        :StandardODEProblem, :UJacobianWrapper, :LinearProblem,
        :ConvexOptimizationProblem,
    )
    @test !any(Base.isexported.(Ref(OrdinaryDiffEqPRK), internal))
end

@time @safetestset "Threaded stages inside another threaded loop" begin
    using OrdinaryDiffEqPRK, Test
    prob_iip = ODEProblem((du, u, p, t) -> (du .= -u), [1.0, 2.0], (0.0, 1.0))
    prob_oop = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    for prob in (prob_iip, prob_oop)
        serial = solve(prob, KuttaPRK2p5(threading = false); dt = 0.1)
        sols = Vector{Any}(undef, 2)
        Threads.@threads for i in 1:2
            sols[i] = solve(prob, KuttaPRK2p5(threading = true); dt = 0.1)
        end
        @test all(sol -> sol.u == serial.u, sols)
    end
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
