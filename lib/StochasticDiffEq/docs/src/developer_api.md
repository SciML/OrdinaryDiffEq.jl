# Developer API

These names are exported for the solver subpackages in the OrdinaryDiffEq
monorepo. They are documented and versioned so downstream solver developers can
extend the stochastic solver infrastructure, but they are not general
user-facing API. User code should prefer the documented solver constructors and
the high-level `solve` interface.

## Extension contract

Solver subpackages extend this layer through a small set of generic dispatches:

- An iterated-integral evaluator is selected with [`StochasticDiffEqCore.get_Jalg`](@ref).
  Implementations must provide the matching out-of-place
  [`StochasticDiffEqCore.get_iterated_I`](@ref) method, or the in-place
  [`StochasticDiffEqCore.get_iterated_I!`](@ref) method that writes into the
  evaluator's preallocated `J` field. Use `AbstractJDiagonal` only for diagonal
  noise and `AbstractJCommute` only when the noise is commutative.
- An SDE algorithm's cache is created through
  [`StochasticDiffEqCore.alg_cache`](@ref). A solver package supplies one method
  for its algorithm and returns the cache consumed by `perform_step!`; the
  `Val{iip}` argument must select the matching in-place or out-of-place layout.
  Cache accessors must preserve the resize and random-number-buffer conventions
  used by the generic cache utilities.
- Algorithm-selection methods such as
  [`StochasticDiffEqCore.alg_compatible`](@ref),
  [`StochasticDiffEqCore.alg_needs_extra_process`](@ref), and
  [`StochasticDiffEqCore.alg_stability_size`](@ref) describe properties of the
  algorithm, not the problem's user-facing solve interface. Add methods for a
  new algorithm only when the generic solver calls them for that property.

These methods are extension points between solver packages. They are versioned
developer API, but applications should use the solver constructors and
`solve(prob, alg)` rather than calling them directly.

```@docs
StochasticDiffEqCore.AbstractJ
StochasticDiffEqCore.AbstractJCommute
StochasticDiffEqCore.AbstractJDiagonal
StochasticDiffEqCore.AutoAlgSwitch
StochasticDiffEqCore.AutoSwitch
StochasticDiffEqCore.DiffCache
StochasticDiffEqCore.DiffEqNLSolveTag
StochasticDiffEqCore.IICommutative
StochasticDiffEqCore.IIFNLSolveFunc
StochasticDiffEqCore.IILevyArea
StochasticDiffEqCore.Ihat2
StochasticDiffEqCore.IteratedIntegralAlgorithm_iip
StochasticDiffEqCore.IteratedIntegralApprox
StochasticDiffEqCore.JCommute_iip
StochasticDiffEqCore.JCommute_oop
StochasticDiffEqCore.JDiagonal_iip
StochasticDiffEqCore.JDiagonal_oop
StochasticDiffEqCore.NLSOLVEJL_SETUP
StochasticDiffEqCore.SDEAlgTypes
StochasticDiffEqCore.SDEIntegrator
StochasticDiffEqCore.SDEOptions
StochasticDiffEqCore.StochasticCompositeCache
StochasticDiffEqCore.TauLeapingDrift
StochasticDiffEqCore._resolve_rng
StochasticDiffEqCore._sde_init
StochasticDiffEqCore._z_prototype
StochasticDiffEqCore.addat_noise!
StochasticDiffEqCore.alg_cache
StochasticDiffEqCore.alg_can_repeat_jac
StochasticDiffEqCore.alg_compatible
StochasticDiffEqCore.alg_control_rate
StochasticDiffEqCore.alg_mass_matrix_compatible
StochasticDiffEqCore.alg_needs_extra_process
StochasticDiffEqCore.alg_stability_size
StochasticDiffEqCore.calc_threepoint_random
StochasticDiffEqCore.calc_threepoint_random!
StochasticDiffEqCore.calc_twopoint_random
StochasticDiffEqCore.calc_twopoint_random!
StochasticDiffEqCore.concrete_prob
StochasticDiffEqCore.deleteat_noise!
StochasticDiffEqCore.delta_default
StochasticDiffEqCore.determine_chunksize
StochasticDiffEqCore.fill_new_noise_caches!
StochasticDiffEqCore.get_Jalg
StochasticDiffEqCore.get_chunksize
StochasticDiffEqCore.get_current_alg_order
StochasticDiffEqCore.get_iterated_I
StochasticDiffEqCore.get_iterated_I!
StochasticDiffEqCore.is_split_step
StochasticDiffEqCore.resize_noise!
StochasticDiffEqCore.supports_regular_jumps
StochasticDiffEqCore.unwrap_alg
StochasticDiffEqHighOrder.RosslerSRA
StochasticDiffEqHighOrder.RosslerSRI
StochasticDiffEqHighOrder.checkSRAOrder
StochasticDiffEqHighOrder.checkSRIOrder
StochasticDiffEqHighOrder.constructExplicitSKenCarp
StochasticDiffEqHighOrder.constructSKenCarp
StochasticDiffEqHighOrder.constructSOSRA
StochasticDiffEqHighOrder.constructSOSRA2
StochasticDiffEqHighOrder.constructSRA1
StochasticDiffEqHighOrder.constructSRA2
StochasticDiffEqHighOrder.constructSRA3
StochasticDiffEqHighOrder.constructSRIOpt1
StochasticDiffEqHighOrder.constructSRIOpt2
StochasticDiffEqHighOrder.constructSRIW1
StochasticDiffEqHighOrder.constructSRIW2
StochasticDiffEqHighOrder.du_cache
StochasticDiffEqHighOrder.u_cache
StochasticDiffEqHighOrder.user_cache
StochasticDiffEqLeaping.ImplicitTauLeaping
StochasticDiffEqLeaping.ThetaTrapezoidalTauLeaping
StochasticDiffEqWeak.IRI1
StochasticDiffEqWeak.KomoriNON
StochasticDiffEqWeak.KomoriNON2
StochasticDiffEqWeak.RDI1WM
StochasticDiffEqWeak.RoesslerRI
StochasticDiffEqWeak.RoesslerRS
StochasticDiffEqWeak.checkNONOrder
StochasticDiffEqWeak.checkRIOrder
StochasticDiffEqWeak.checkRSOrder
StochasticDiffEqWeak.constructDRI1
StochasticDiffEqWeak.constructNON
StochasticDiffEqWeak.constructNON2
StochasticDiffEqWeak.constructRDI1WM
StochasticDiffEqWeak.constructRDI2WM
StochasticDiffEqWeak.constructRDI3WM
StochasticDiffEqWeak.constructRDI4WM
StochasticDiffEqWeak.constructRI1
StochasticDiffEqWeak.constructRI3
StochasticDiffEqWeak.constructRI5
StochasticDiffEqWeak.constructRI6
StochasticDiffEqWeak.constructRS1
StochasticDiffEqWeak.constructRS2
```
