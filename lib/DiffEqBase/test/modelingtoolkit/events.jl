using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D, SymbolicContinuousCallback

@testset "All simultaenous events are saved" begin
    # https://github.com/SciML/ModelingToolkit.jl/issues/4870 Case 2
    @variables x(t) = 0.0
    c_up = only(@discretes c_up(t) = 0.0)
    c_dn = only(@discretes c_dn(t) = 0.0)
    cb_up = SymbolicContinuousCallback([x ~ 0.5], [c_up => c_up + 1]; affect_neg = nothing)
    cb_dn = SymbolicContinuousCallback([x ~ 0.5], nothing; affect_neg = [c_dn => c_dn + 1])
    @mtkcompile sys2 = System([D(x) ~ cos(t)], t, [x], [c_up, c_dn]; continuous_events = [cb_up, cb_dn])
    prob2 = ODEProblem(sys2, [], (0.0, 13.0))
    # Ensure we have two identical events
    @test prob2.kwargs[:callback].len == 2
    buffer = [NaN, Inf]
    prob2.kwargs[:callback].condition(buffer, prob2.u0, prob2.p, init(prob2, Tsit5()))
    @test allequal(buffer)

    sol2 = solve(prob2, Tsit5())
    @test length(sol2[c_up]) == length(sol2[c_dn])
    @test sol2[c_up][end] == sol2[c_dn][end] == 2
end
