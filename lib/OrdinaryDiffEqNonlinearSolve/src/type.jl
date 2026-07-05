# algorithms
"""
    NLFunctional(; κ = 1//100, max_iter = 10, fast_convergence_cutoff = 1//5)

Functional (fixed-point) iteration solver for the implicit stage equations,
`z ← g(z)`. No Jacobian/`W` is formed, so it is cheap per iteration but only
converges for mildly stiff problems.

# Keyword arguments

  - `κ`: relative tolerance on the increment used in the convergence test.
  - `max_iter`: maximum number of fixed-point iterations per solve.
  - `fast_convergence_cutoff`: convergence-rate threshold for fast convergence.
"""
struct NLFunctional{K, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
end

function NLFunctional(; κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5)
    return NLFunctional(κ, fast_convergence_cutoff, max_iter)
end

"""
    NLAnderson(; κ = 1//100, max_iter = 10, max_history = 5, aa_start = 1,
               droptol = nothing, fast_convergence_cutoff = 1//5)

Anderson-accelerated fixed-point iteration for the implicit stage equations. Like
[`NLFunctional`](@ref) but mixes in `max_history` previous residuals via a
least-squares update to accelerate convergence.

# Keyword arguments

  - `κ`, `max_iter`, `fast_convergence_cutoff`: as in [`NLFunctional`](@ref).
  - `max_history`: number of past iterates kept for the acceleration.
  - `aa_start`: iteration at which acceleration starts.
  - `droptol`: optional condition-number threshold for dropping history columns.
"""
struct NLAnderson{K, D, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
    max_history::Int
    aa_start::Int
    droptol::D
end

function NLAnderson(;
        κ = 1 // 100, max_iter = 10, max_history::Int = 5, aa_start::Int = 1,
        droptol = nothing, fast_convergence_cutoff = 1 // 5
    )
    return NLAnderson(κ, fast_convergence_cutoff, max_iter, max_history, aa_start, droptol)
end

"""
    NLNewton(; κ = 1//100, max_iter = 10, fast_convergence_cutoff = 1//5,
             new_W_dt_cutoff = 1//5, always_new = false, check_div = true, relax = nothing)

Quasi-Newton nonlinear solver for the implicit stage equations. Uses the
`W = M/(γΔt) - J` matrix (reused/refactorized across steps and stages when
possible) to solve `g(z) = 0`.

# Keyword arguments

  - `κ`: relative tolerance on the Newton increment used in the convergence test.
  - `max_iter`: maximum number of Newton iterations per solve.
  - `fast_convergence_cutoff`: convergence-rate threshold below which convergence
    is deemed fast (allowing `W` reuse).
  - `new_W_dt_cutoff`: relative change in `γΔt` above which `W` is refactorized.
  - `always_new`: force recomputation of `W` on every solve.
  - `check_div`: enable early divergence detection.
  - `relax`: optional relaxation parameter in `[0, 1)` damping the Newton update.
"""
struct NLNewton{K, C1, C2, R} <: AbstractNLSolverAlgorithm
    κ::K
    max_iter::Int
    fast_convergence_cutoff::C1
    new_W_dt_cutoff::C2
    always_new::Bool
    check_div::Bool
    relax::R
end

function NLNewton(;
        κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true,
        relax = nothing
    )
    if relax isa Number && !(0 <= relax < 1)
        throw(ArgumentError("The relaxation parameter must be in [0, 1), got `relax = $relax`"))
    end

    return NLNewton(
        κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff, always_new, check_div,
        relax
    )
end

struct NonlinearSolveAlg{K, C1, C2, A} <: AbstractNLSolverAlgorithm
    κ::K
    max_iter::Int
    fast_convergence_cutoff::C1
    new_W_dt_cutoff::C2
    always_new::Bool
    check_div::Bool
    alg::A
end

function NonlinearSolveAlg(
        alg = NewtonRaphson(autodiff = AutoFiniteDiff());
        κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true
    )
    return NonlinearSolveAlg(
        κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff, always_new, check_div,
        alg
    )
end

"""
    HomotopyNonlinearSolveAlg(alg = HomotopySweep(inner = NewtonRaphson(autodiff = AutoFiniteDiff()));
        κ = 1//100, max_iter = 10, fast_convergence_cutoff = 1//5,
        abstol = nothing, reltol = nothing)

Solve the implicit stage equations by homotopy continuation in the step size instead of
a plain Newton iteration. The stage equation

```math
0 = dt⋅f(\\mathrm{tmp} + γ⋅z, p, t + c⋅dt) - M z
```

is embedded into the one-parameter family

```math
H(z, λ) = λ⋅dt⋅f(\\mathrm{tmp} + γ⋅z, p, t + c⋅dt) - M z, \\quad λ ∈ (0, 1)
```

(with the analogous embedding `tmp + λ⋅f(z) - α/(γ dt)⋅M z` for multistep-form
methods), built as a `SciMLBase.HomotopyProblem` and handed to a NonlinearSolve.jl
continuation solver. At `λ = 0` the solution is known (`z = 0` for DIRK-form methods, a
single linear solve for multistep form), and the continuation tracks the solution branch
of the implicit method as the effective step size grows from `0` to `dt` — the
"principal branch" of Green, Patrick & Spiteri (*On theoretical upper limits for valid
timesteps of implicit ODE methods*, AIMS Mathematics 4(6), 2019). This is much more
expensive than [`NLNewton`](@ref) per stage, but it converges from the exact anchor at
`λ = 0` regardless of predictor quality, which makes it useful on problems where the
Newton iteration fails even after step-size reduction.

When the continuation cannot reach `λ = 1` (for example at a fold, where the connected
solution branch of the stage equation turns back and no consistent solution at the full
`dt` exists), the solve reports divergence and the step is rejected, so the integrator
retries with a smaller `dt` — exactly the semantically correct response to a fold in the
step-size homotopy.

# Positional arguments

  - `alg`: the continuation algorithm used for the `HomotopyProblem`. Defaults to
    `HomotopySweep(inner = NewtonRaphson(autodiff = AutoFiniteDiff()))`
    (natural-parameter continuation). Any NonlinearSolve.jl homotopy solver works, e.g.
    `ArcLengthContinuation(inner = NewtonRaphson(autodiff = AutoFiniteDiff()))` to track
    the branch around folds. Note that the inner solver must not use ForwardDiff-based
    autodiff: the stage residual closes over preallocated `Float64` buffers (the same
    restriction as `NonlinearSolveAlg`).

# Keyword arguments

  - `κ`, `max_iter`, `fast_convergence_cutoff`: kept for interface compatibility with the
    other nonlinear-solver algorithms; the continuation solve does not run the shared
    Newton convergence loop, so only `max_iter` (as a bound on inner iterations reported
    to the integrator statistics) has an effect.
  - `abstol`, `reltol`: tolerances passed to the continuation `solve`. The residual is
    kept in `u` units, and `abstol = nothing` (the default) resolves at solve time to
    `κ * abstol_integrator`, mirroring the `κ ⋅ tol` convergence criterion of the Newton
    solvers. `reltol = nothing` uses the NonlinearSolve.jl default.

!!! warning

    DAE problems (`DAEFunction`s and singular mass matrices) are not supported: the
    λ-embedding scales the whole right-hand side, which degenerates the algebraic
    equations at `λ = 0`.
"""
struct HomotopyNonlinearSolveAlg{K, C1, A, AT, RT} <: AbstractNLSolverAlgorithm
    κ::K
    max_iter::Int
    fast_convergence_cutoff::C1
    abstol::AT
    reltol::RT
    alg::A
end

function HomotopyNonlinearSolveAlg(
        alg = HomotopySweep(inner = NewtonRaphson(autodiff = AutoFiniteDiff()));
        κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        abstol = nothing, reltol = nothing
    )
    return HomotopyNonlinearSolveAlg(
        κ, max_iter, fast_convergence_cutoff, abstol, reltol, alg
    )
end

# solver

mutable struct NLSolver{
        algType, iip, uType, gamType, tmpType, tType,
        C <: AbstractNLSolverCache, E,
    } <: AbstractNLSolver{algType, iip}
    z::uType
    tmp::uType # DIRK and multistep methods only use tmp
    tmp2::tmpType # for GLM if necessary
    ztmp::uType
    γ::gamType
    c::tType
    α::tType
    alg::algType
    κ::tType
    fast_convergence_cutoff::tType
    ηold::tType
    iter::Int
    maxiters::Int
    status::NLStatus
    cache::C
    method::MethodType
    nfails::Int
    prev_θ::E
end

# default to DIRK
function NLSolver{iip, tType}(
        z, tmp, ztmp, γ, c, α, alg, κ, fast_convergence_cutoff, ηold,
        iter, maxiters, status, cache, method = DIRK, tmp2 = nothing,
        nfails::Int = 0
    ) where {iip, tType}
    RT = real(eltype(z))
    return NLSolver{typeof(alg), iip, typeof(z), typeof(γ), typeof(tmp2), tType, typeof(cache), RT}(
        z,
        tmp,
        tmp2,
        ztmp,
        float(γ),
        convert(tType, c),
        convert(tType, α),
        alg,
        convert(tType, κ),
        convert(tType, fast_convergence_cutoff),
        convert(tType, ηold),
        iter,
        maxiters,
        status,
        cache,
        method,
        nfails,
        one(RT)
    )
end

# caches

mutable struct NLNewtonCache{
        uType,
        tType,
        tType2,
        rateType,
        J,
        W,
        ufType,
        jcType,
        lsType,
        DC,
    } <: AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    dz::uType
    J::J
    W::W
    new_W::Bool
    firststage::Bool
    firstcall::Bool
    W_γdt::tType
    du1::uType
    uf::ufType
    jac_config::jcType
    linsolve::lsType
    weight::uType
    invγdt::tType2
    new_W_γdt_cutoff::tType
    J_t::tType
    dae_jacobians::DC
end

mutable struct NLNewtonConstantCache{tType, tType2, J, W, ufType, DC} <: AbstractNLSolverCache
    tstep::tType
    J::J
    W::W
    new_W::Bool
    firststage::Bool
    firstcall::Bool
    W_γdt::tType
    uf::ufType
    invγdt::tType2
    new_W_γdt_cutoff::tType
    J_t::tType
    dae_jacobians::DC
end

mutable struct NLFunctionalCache{uType, tType, rateType} <: AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    dz::uType
end

mutable struct NLFunctionalConstantCache{tType} <: AbstractNLSolverCache
    tstep::tType
end

mutable struct NLAndersonCache{uType, tType, rateType, uEltypeNoUnits} <:
    AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    dz::uType
    """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
    dzold::uType
    """value `g(zprev)` of previous fixed-point iteration"""
    z₊old::uType
    Δz₊s::Vector{uType}
    Q::Matrix{uEltypeNoUnits}
    R::Matrix{uEltypeNoUnits}
    γs::Vector{uEltypeNoUnits}
    history::Int
    aa_start::Int
    droptol::Union{Nothing, tType}
end

mutable struct NLAndersonConstantCache{uType, tType, uEltypeNoUnits} <:
    AbstractNLSolverCache
    tstep::tType
    dz::uType
    """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
    dzold::uType
    """value `g(zprev)` of previous fixed-point iteration"""
    z₊old::uType
    Δz₊s::Vector{uType}
    Q::Matrix{uEltypeNoUnits}
    R::Matrix{uEltypeNoUnits}
    γs::Vector{uEltypeNoUnits}
    history::Int
    aa_start::Int
    droptol::Union{Nothing, tType}
end

mutable struct HomotopyNonlinearSolveCache{uType, tType, rateType, tType2, F, R} <:
    AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    invγdt::tType2
    nlfunc::F
    # residual-evaluation counter (a `Ref(0)`): the continuation solutions do not carry
    # stats, so the residual itself counts its `f` calls for the integrator statistics
    nf::R
end

# Analytic-jac callback handing the ODE-side W to the inner NonlinearSolve under
# W-reuse. Deliberately a struct holding a Ref rather than a closure: its type depends
# only on W's type, not on a captured binding, so when `resize!` replaces the W array
# the same NonlinearFunction/NonlinearProblem types remain valid and `initialize!` just
# swaps the Ref target.
struct WReuseJac{W} <: Function
    W::Base.RefValue{W}
end
(j::WReuseJac)(J_out, z, p) = (copyto!(J_out, j.W[]); J_out)

mutable struct NonlinearSolveCache{uType, tType, rateType, tType2, P, C, JType, WType, ufType, jcType, du1Type, weightType} <:
    AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    invγdt::tType2
    prob::P
    cache::C
    J::JType
    W::WType
    uf::ufType
    jac_config::jcType
    du1::du1Type
    weight::weightType
    W_γdt::tType
    new_W::Bool
end
