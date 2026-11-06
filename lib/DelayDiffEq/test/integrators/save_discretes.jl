using DelayDiffEq, OrdinaryDiffEqTsit5, SciMLBase, Test
using RecursiveArrayTools: DiffEqArray
using SymbolicIndexingInterface
using SymbolicIndexingInterface:
    ParameterTimeseriesCollection, ParameterTimeseriesIndex, SymbolCache

struct SaveDiscretesDDESys
    sc::SymbolCache
end
SymbolicIndexingInterface.symbolic_container(s::SaveDiscretesDDESys) = s.sc

function SciMLBase.create_parameter_timeseries_collection(::SaveDiscretesDDESys, ps, tspan)
    dea = DiffEqArray(Vector{Float64}[], Float64[])
    return ParameterTimeseriesCollection((dea,), deepcopy(ps))
end
function SciMLBase.get_saveable_values(::SaveDiscretesDDESys, p::Vector{Float64}, tsidx)
    return [p[1]]
end

function save_discretes_dde_prob()
    sc = SymbolCache(
        [:x], Dict(:h => 1), :t;
        timeseries_parameters = Dict(:h => ParameterTimeseriesIndex(1, 1))
    )
    sys = SaveDiscretesDDESys(sc)
    f!(du, u, h, p, t) = (du[1] = -u[1] + p[1]; nothing)
    hist(p, t) = [1.0]
    return DDEProblem(
        DDEFunction(f!; sys = sys), [1.0], hist, (0.0, 0.2), [0.0];
        constant_lags = [0.1]
    )
end

@testset "save_discretes=false skips DDE initialization discrete save" begin
    prob = save_discretes_dde_prob()
    alg = MethodOfSteps(Tsit5())
    cb = DiscreteCallback(
        (u, t, integ) -> false,
        integ -> nothing;
        save_positions = (false, true),
        saved_clock_partitions = (1,),
    )

    sol_default = solve(deepcopy(prob), alg; callback = cb)
    @test !isempty(sol_default.discretes.collection[1].t)

    sol_false = solve(deepcopy(prob), alg; callback = cb, save_discretes = false)
    @test isempty(sol_false.discretes.collection[1].t)
    @test isempty(sol_false.discretes.collection[1].u)

    integrator = init(deepcopy(prob), alg; callback = cb, save_discretes = false)
    solve!(integrator)
    @test isempty(integrator.sol.discretes.collection[1].t)
    reinit!(integrator)
    solve!(integrator)
    @test isempty(integrator.sol.discretes.collection[1].t)
    @test isempty(integrator.sol.discretes.collection[1].u)
end
