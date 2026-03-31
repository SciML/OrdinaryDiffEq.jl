using FlexUnits, FlexUnits.UnitRegistry, OrdinaryDiffEq, Test

f(du, u, p, t) = du .= 3 * u"1/s" * u
prob = ODEProblem(f, [2.0u"m"], (0.0u"s", 1.0u"s"))
intg = init(prob, Tsit5(), dt = 0.01u"s")
@test_nowarn step!(intg, 0.02u"s", true)

@test SciMLBase.unitfulvalue(1.0u"1/s") == 1.0u"1/s"
@test SciMLBase.value(1.0u"1/s") isa Real
