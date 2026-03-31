using Unitful, OrdinaryDiffEq, ForwardDiff, Test
f(du, u, p, t) = du .= 3 * u"1/s" * u
prob = ODEProblem(f, [2.0u"m"], (0.0u"s", Inf * u"s"))
intg = init(prob, Tsit5())
@test_nowarn step!(intg, 0.02u"s", true)

@test SciMLBase.unitfulvalue(u"1/s") == u"1/s"
@test SciMLBase.value(ForwardDiff.Dual(1) * u"1/s") == 1
@test SciMLBase.unitfulvalue(ForwardDiff.Dual(1) * u"1/s") == u"1/s"
