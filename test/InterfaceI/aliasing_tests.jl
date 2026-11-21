using OrdinaryDiffEq, Test

import ODEProblemLibrary: prob_ode_linear

# Test that ODEAliasSpecifier works for alias_u0
u0_alias_sol = solve(
    prob_ode_linear, Tsit5(), alias = ODEAliasSpecifier(alias_u0 = true)
)
u0_no_alias_sol = solve(prob_ode_linear, Tsit5())

@test u0_alias_sol.u == u0_no_alias_sol.u
