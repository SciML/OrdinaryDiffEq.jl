# Adding Algorithms

New algorithms can either be added by extending one of the current solver (or
add-on packages), or by contributing a new package to the organization. If
it's a new problem (a new PDE, a new type of differential equation, a new
subclass of problems for which special methods exist, etc.) then the problem
and solution types should be added to `DiffEqBase` first.

After the problem and solutions are defined, the `__solve` method should be implemented.
It should take in keyword arguments which match the common interface (implement
"as many as possible"). One should note and document the amount of compatibility
with the common interface and Julia-defined types. After that, testing should be
done using `DiffEqDevTools`. Convergence tests and benchmarks should be included
to show the effectiveness of the algorithm and the correctness. Do not worry if
the algorithm is not "effective": the implementation can improve over time and
some algorithms are useful just for the comparison they give!

After some development, one may want to document the algorithm in DiffEqBenchmarks
and DiffEqTutorials.

## Adding new algorithms to OrdinaryDiffEq

This recipe has been used to add the strong stability preserving Runge-Kutta methods
`SSPRK22`, `SSPRK33`, and `SSPRK104` to `OrdinaryDiffEq`. `SSPRK22` will be used
as an example.

  - To create a new solver, two (three) types have to be created.
    The first is the algorithm `SSPRK22` used for dispatch, the other ones are
    the corresponding caches `SSPRK22Cache` (for inplace updates) and
    `SSPRK22ConstantCache`.
  - The algorithm is defined in `algorithms.jl` as
    `struct SSPRK22 <: OrdinaryDiffEqAlgorithm end`.
    Although it does not have the FSAL property, this is set to true since the derivative
    at the start and the end of the interval are used for the Hermite interpolation,
    and so this is FSAL'd so that way only a single extra function evaluation occurs
    over the whole integration. This is done in `alg_utils.jl` via
    `isfsal(alg::SSPRK22) = true`. Additionally, the order is set in the same
    file via `alg_order(alg::SSPRK22) = 2`.
  - The algorithm `SSPRK22` is exported in `OrdinaryDiffEq.jl`.
  - In `caches.jl`, the two cache types `SSPRK22Cache` (for inplace updates) and
    `SSPRK22ConstantCache` are defined, similarly to the other ones.
    Note: `u_cache(c::SSPRK22Cache) = ()` and
    `du_cache(c::SSPRK22Cache) = (c.k,c.du,c.fsalfirst)` return the parts of the
    modifiable cache that are changed if the size of the ODE changes.
  - A new file `perform_step/ssprk_perform_step.jl` has been used for the new
    implementations. For both types of caches, the functions `initialize!`
    and `perform_step!` are defined there.
  - Finally, tests are added. A new file `test/ode/ode_ssprk_tests.jl` is created
    and included in `tests/runtests.jl` via
    `@time @testset "SSPRK Tests" begin include("ode/ode_ssprk_tests.jl") end`.
  - Additionally, regression tests for the dense output are added in
    `test/ode/ode_dense_tests.jl`.

For more details, refer to https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/pull/40

### Self-Contained Example

```julia
using OrdinaryDiffEq
import OrdinaryDiffEq:
                       OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache,
                       OrdinaryDiffEqConstantCache,
                       alg_order, alg_cache, initialize!, perform_step!, trivial_limiter!,
                       constvalue,
                       @muladd, @unpack, @cache, @..

struct RK_ALG{StageLimiter, StepLimiter} <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
end
RK_ALG(stage_limiter! = trivial_limiter!) = RK_ALG(stage_limiter!, trivial_limiter!)
export RK_ALG
alg_order(alg::RK_ALG) = 3

@cache struct RK_ALGCache{uType, rateType, StageLimiter, StepLimiter, TabType} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    tmp::uType
    u₂::uType
    fsalfirst::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    tab::TabType
end

struct RK_ALGConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    α40::T
    α41::T
    α43::T
    α62::T
    α65::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    β65::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
end

function RK_ALGConstantCache(T, T2)
    α40 = T(0.476769811285196)
    α41 = T(0.098511733286064)
    α43 = T(0.424718455428740)
    α62 = T(0.155221702560091)
    α65 = T(0.844778297439909)
    β10 = T(0.284220721334261)
    β21 = T(0.284220721334261)
    β32 = T(0.284220721334261)
    β43 = T(0.120713785765930)
    β54 = T(0.284220721334261)
    β65 = T(0.240103497065900)
    c1 = T2(0.284220721334261)
    c2 = T2(0.568441442668522)
    c3 = T2(0.852662164002783)
    c4 = T2(0.510854218958172)
    c5 = T2(0.795074940292433)

    RK_ALGConstantCache(
        α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4, c5)
end

function alg_cache(alg::RK_ALG, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{true})
    tmp = similar(u)
    u₂ = similar(u)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    tab = RK_ALGConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    RK_ALGCache(u, uprev, k, tmp, u₂, fsalfirst, alg.stage_limiter!, alg.step_limiter!, tab)
end

function alg_cache(alg::RK_ALG, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{false})
    RK_ALGConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function initialize!(integrator, cache::RK_ALGConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.destats.nf += 1
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::RK_ALGConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4,
    c5 = cache

    # u1 -> stored as u
    u = uprev + β10 * dt * integrator.fsalfirst
    k = f(u, p, t + c1 * dt)
    # u2
    u₂ = u + β21 * dt * k
    k = f(u₂, p, t + c2 * dt)
    # u3
    tmp = u₂ + β32 * dt * k
    k = f(tmp, p, t + c3 * dt)
    # u4
    tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
    k = f(tmp, p, t + c4 * dt)
    # u5
    tmp = tmp + β54 * dt * k
    k = f(tmp, p, t + c5 * dt)
    # u
    u = α62 * u₂ + α65 * tmp + β65 * dt * k

    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    integrator.destats.nf += 6
    integrator.k[1] = integrator.fsalfirst
    integrator.u = u
end

function initialize!(integrator, cache::RK_ALGCache)
    @unpack k, fsalfirst = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::RK_ALGCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k, tmp, u₂, fsalfirst, stage_limiter!, step_limiter! = cache
    @unpack α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4,
    c5 = cache.tab

    # u1 -> stored as u
    @.. u = uprev + β10 * dt * integrator.fsalfirst
    stage_limiter!(u, f, p, t + c1 * dt)
    f(k, u, p, t + c1 * dt)
    # u2
    @.. u₂ = u + β21 * dt * k
    stage_limiter!(u₂, f, p, t + c2 * dt)
    f(k, u₂, p, t + c2 * dt)
    # u3
    @.. tmp = u₂ + β32 * dt * k
    stage_limiter!(tmp, f, p, t + c3 * dt)
    f(k, tmp, p, t + c3 * dt)
    # u4
    @.. tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
    stage_limiter!(tmp, f, p, t + c4 * dt)
    f(k, tmp, p, t + c4 * dt)
    # u5
    @.. tmp = tmp + β54 * dt * k
    stage_limiter!(tmp, f, p, t + c5 * dt)
    f(k, tmp, p, t + c5 * dt)
    # u
    @.. u = α62 * u₂ + α65 * tmp + β65 * dt * k
    stage_limiter!(u, f, p, t + dt)
    step_limiter!(u, f, p, t + dt)
    integrator.destats.nf += 6
    f(k, u, p, t + dt)
end

#oop test
f = ODEFunction((u, p, t) -> 1.01u,
    analytic = (u0, p, t) -> u0 * exp(1.01t))
prob = ODEProblem(f, 1.01, (0.0, 1.0))
sol = solve(prob, RK_ALG(), dt = 0.1)

using Plots
plot(sol)
plot(sol, denseplot = false, plot_analytic = true)

using DiffEqDevTools
dts = (1 / 2) .^ (8:-1:1)
sim = test_convergence(dts, prob, RK_ALG())
sim.𝒪est[:final]
plot(sim)

# Example of a good one!
sim = test_convergence(dts, prob, BS3())
sim.𝒪est[:final]
plot(sim)

#iip test
f = ODEFunction((du, u, p, t) -> (du .= 1.01 .* u),
    analytic = (u0, p, t) -> u0 * exp(1.01t))
prob = ODEProblem(f, [1.01], (0.0, 1.0))
sol = solve(prob, RK_ALG(), dt = 0.1)

plot(sol)
plot(sol, denseplot = false, plot_analytic = true)

dts = (1 / 2) .^ (8:-1:1)
sim = test_convergence(dts, prob, RK_ALG())
sim.𝒪est[:final]
plot(sim)

# Example of a good one!
sim = test_convergence(dts, prob, BS3())
sim.𝒪est[:final]
plot(sim)
```

### Adding new exponential algorithms

The exponential algorithms follow the same recipe as the general algorithms, but there
are automation utilities that make this easier. It is recommended that you refer to one
of the model algorithms for reference:

  - For traditional exponential Runge-Kutta type methods (that come with a corresponding
    Butcher table), refer to `ETDRK2`.
  - For adaptive exponential Rosenbrock type methods, refer to `Exprb32`.
  - For exponential propagation iterative Runge-Kutta methods (EPIRK), refer to `EPIRK5P1`.

The first two classes support two modes of operation: operator caching and Krylov
approximation. The `perform_step!` method in `perform_step/exponential_rk_perform_step.jl`,
as a result, is split into two branches depending on whether `alg.krylov` is true. The
caching branch utilizes precomputed operators, which are calculated by the `expRK_operators`
method in `caches/linear_nonlinear_caches.jl`. Both `expRK_operators` and the `arnoldi`/`phiv`
methods in `perform_step!` comes from the
[ExponentialUtilities](https://github.com/JuliaDiffEq/ExponentialUtilities.jl) package.

The EPIRK methods can only use Krylov approximation, and unlike the previous two they use
the timestepping variant `phiv_timestep`. The timestepping method follows the convention
of Neisen & Wright, and can be toggled to use adaptation by `alg.adaptive_krylov`.

Although the exponential integrators (especially the in-place version) can seem complex, they
share similar structures. The infrastructure for the existing exponential methods utilize the
fact to reduce boilerplate code. In particular, the cache construction code in
`caches/linear_nonlinear_caches.jl` and the `initialize!` method in
`perform_step/exponential_rk_perform_step.jl` can be mostly automated and only `perform_step!`
needs implementing.

Finally, to construct tests for the new exponential algorithm, append the new algorithm to
the corresponding algorithm class in `test/linear_nonlinear_convergence_tests.jl` and
`test/linear_nonlinear_krylov_tests.jl`.
