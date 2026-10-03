using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, Test
using SymbolicIndexingInterface
using SymbolicIndexingInterface:
    ParameterTimeseriesCollection, ParameterTimeseriesIndex, SymbolCache

const DiffEqArray = OrdinaryDiffEqCore.DiffEqArray

# Minimal index provider with one discrete-parameter partition so initialization
# can save through `saved_clock_partitions` without ModelingToolkit.
struct SaveDiscretesSys
    sc::SymbolCache
end
SymbolicIndexingInterface.symbolic_container(s::SaveDiscretesSys) = s.sc

function SciMLBase.create_parameter_timeseries_collection(::SaveDiscretesSys, ps, tspan)
    dea = DiffEqArray(Vector{Float64}[], Float64[])
    return ParameterTimeseriesCollection((dea,), deepcopy(ps))
end
function SciMLBase.get_saveable_values(::SaveDiscretesSys, p::Vector{Float64}, tsidx)
    return [p[1]]
end

@testset "save_discretes=false skips initialization discrete save" begin
    sc = SymbolCache([:x], Dict(:h => 1), :t;
        timeseries_parameters = Dict(:h => ParameterTimeseriesIndex(1, 1)))
    sys = SaveDiscretesSys(sc)
    f!(du, u, p, t) = (du[1] = -u[1] + p[1]; nothing)
    prob = ODEProblem(ODEFunction(f!; sys = sys), [1.0], (0.0, 0.2), [0.0])
    # Never-firing callback: only the initialize_callbacks! discrete save can populate storage.
    cb = DiscreteCallback(
        (u, t, integ) -> false,
        integ -> nothing;
        save_positions = (false, true),
        saved_clock_partitions = (1,),
    )

    sol_true = solve(deepcopy(prob), Tsit5(); callback = cb, save_discretes = true)
    @test sol_true.discretes.collection[1].t == [0.0]
    @test sol_true.discretes.collection[1].u == [[0.0]]

    sol_false = solve(deepcopy(prob), Tsit5(); callback = cb, save_discretes = false)
    @test isempty(sol_false.discretes.collection[1].t)
    @test isempty(sol_false.discretes.collection[1].u)
end
