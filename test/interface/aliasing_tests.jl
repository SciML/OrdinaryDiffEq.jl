using OrdinaryDiffEq, Test

import ODEProblemLibrary: prob_ode_linear

# Test that the old keyword works, and that the new AliasSpecier works.
u0_old_alias_kwarg_sol = solve(prob_ode_linear, Tsit5(), alias_u0 = true)
u0_new_alias_kwarg_sol = solve(
    prob_ode_linear, Tsit5(), alias = ODEAliasSpecifier(alias_u0 = true))

@test u0_old_alias_kwarg_sol == u0_new_alias_kwarg_sol
