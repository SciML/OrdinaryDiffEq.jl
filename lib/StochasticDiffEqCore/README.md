# StochasticDiffEqCore.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqCore.jl is the core infrastructure package for stochastic differential equation (SDE) solvers within the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo. It provides the shared infrastructure used by all SDE solver subpackages, including:

- `__init` and `__solve` entry points for `SDEProblem` and `RODEProblem`
- Abstract algorithm types (`StochasticDiffEqAlgorithm`, `StochasticDiffEqAdaptiveAlgorithm`, etc.)
- Algorithm traits (`alg_order`, `isadaptive`, `alg_compatible`, etc.)
- `SDEIntegrator` type alias and SDE-specific integrator utilities
- Iterated integral approximation interface
- Composite algorithm support (`StochasticCompositeAlgorithm`)

This package is not intended to be used directly by end users. Instead, use one of the solver subpackages (e.g., `StochasticDiffEqLowOrder`, `StochasticDiffEqHighOrder`) or the umbrella package [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solver Subpackages

| Package | Solvers |
|---------|---------|
| StochasticDiffEqLowOrder | EM, EulerHeun, LambaEM, LambaEulerHeun, SimplifiedEM, SplitEM, RKMil, RKMilCommute, PCEuler |
| StochasticDiffEqHighOrder | SRI, SRIW1, SRIW2, SOSRI, SOSRI2, SRA, SRA1, SRA2, SRA3, SOSRA, SOSRA2 |
| StochasticDiffEqMilstein | RKMilGeneral, WangLi3SMil_A-F |
| StochasticDiffEqROCK | SROCK1, SROCK2, KomBurSROCK2, SROCKC2, SROCKEM, SKSROCK, TangXiaoSROCK2 |
| StochasticDiffEqImplicit | ImplicitEM, ImplicitEulerHeun, ImplicitRKMil, STrapezoid, SImplicitMidpoint, ISSEM, ISSEulerHeun, SKenCarp |
| StochasticDiffEqWeak | DRI1, RI1, RI3, RI5, RI6, RDI1WM-RDI4WM, W2Ito1, RS1, RS2, PL1WM, NON, COM, IRI1, and more |
| StochasticDiffEqIIF | IIF1M, IIF2M, IIF1Mil |
| StochasticDiffEqLeaping | TauLeaping, CaoTauLeaping, ImplicitTauLeaping, ThetaTrapezoidalTauLeaping |
| StochasticDiffEqRODE | RandomEM, RandomHeun, RandomTamedEM, BAOAB |
