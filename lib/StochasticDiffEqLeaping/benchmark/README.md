Compare native mass-action solves with rate/update adapters using a fresh Julia
process for each invocation. Use Julia 1.12 or later to report compilation time.
From the repository root, after instantiating the
leaping package environment:

```sh
julia --project=lib/StochasticDiffEqLeaping lib/StochasticDiffEqLeaping/benchmark/massaction.jl adapter implicit 1000
julia --project=lib/StochasticDiffEqLeaping lib/StochasticDiffEqLeaping/benchmark/massaction.jl native implicit 1000
```

Arguments select `adapter` or `native`, `explicit` or `implicit`, and the reaction
count. Each row reports mode, method, reactions, first-solve seconds, compilation
seconds, median warmed-up seconds, and median allocated bytes. Warmed-up medians
use nine solves. Repeat fresh processes to assess compilation variability.

Both paths use the same mass-action operations; the adapter supplies them through
`RegularJump` callbacks.
