using Test
using StochasticDiffEq
using SDEProblemLibrary: prob_sde_linear
using SciMLBase

# Test that the old keyword works, and that the new AliasSpecier works.
@test_warn x -> contains("`alias_u0`", x) solve(prob_sde_linear, EM(), dt = 0.1, alias_u0 = true)

@test_nowarn solve(prob_sde_linear, EM(), dt = 0.1)

@test_nowarn solve(prob_sde_linear, EM(), dt = 0.1, alias = SciMLBase.SDEAliasSpecifier(alias_u0 = true))

function f3(u, p, t, W)
    return 2u * sin(W)
end
u0 = 1.0
tspan = (0.0, 5.0)
prob_rode = RODEProblem(f3, u0, tspan)

@test_warn x -> contains("`alias_u0`", x) sol = solve(prob_rode, RandomEM(), dt = 1 / 100, alias_u0 = true)

@test_nowarn solve(prob_rode, RandomEM(), dt = 0.1)

@test_nowarn solve(prob_rode, RandomEM(), dt = 0.1, alias = SciMLBase.RODEAliasSpecifier(alias_u0 = true))
