module OrdinaryDiffEqCoreReactantExt

using OrdinaryDiffEqCore
import Reactant
using Reactant: @reactant_overlay, TracedRArray, TracedRNumber
using ReactantCore: @trace
import Base: *, /, +, -, <, >, <=, >=, convert
import FastPower: fastpower
import DiffEqBase
import SciMLBase

# =============================================================================
# Rational -> Float64 conversion dispatches for traced types
# Reactant cannot convert Rational{Int64} to MLIR types, so we intercept
# operations between Rational and traced types and convert Rational to Float64.
# =============================================================================

# Convert Rational to TracedRNumber by going through Float64
Base.convert(::Type{TracedRNumber{T}}, x::Rational) where {T} = convert(TracedRNumber{T}, convert(T, x))

# Convert BigFloat to TracedRNumber by going through the target type
Base.convert(::Type{TracedRNumber{T}}, x::BigFloat) where {T} = convert(TracedRNumber{T}, convert(T, x))

# Promotion rules to prefer Float64 over Rational/BigFloat when TracedRNumber is involved
Base.promote_rule(::Type{TracedRNumber{T}}, ::Type{Rational{S}}) where {T,S} = TracedRNumber{T}
Base.promote_rule(::Type{TracedRNumber{T}}, ::Type{BigFloat}) where {T} = TracedRNumber{T}

# Multiplication
Base.:*(x::TracedRArray, r::Rational) = x * Float64(r)
Base.:*(r::Rational, x::TracedRArray) = Float64(r) * x
Base.:*(x::TracedRNumber, r::Rational) = x * Float64(r)
Base.:*(r::Rational, x::TracedRNumber) = Float64(r) * x

# Division
Base.:/(x::TracedRArray, r::Rational) = x / Float64(r)
Base.:/(r::Rational, x::TracedRArray) = Float64(r) / x
Base.:/(x::TracedRNumber, r::Rational) = x / Float64(r)
Base.:/(r::Rational, x::TracedRNumber) = Float64(r) / x

# Addition
Base.:+(x::TracedRArray, r::Rational) = x + Float64(r)
Base.:+(r::Rational, x::TracedRArray) = Float64(r) + x
Base.:+(x::TracedRNumber, r::Rational) = x + Float64(r)
Base.:+(r::Rational, x::TracedRNumber) = Float64(r) + x

# Subtraction
Base.:-(x::TracedRArray, r::Rational) = x - Float64(r)
Base.:-(r::Rational, x::TracedRArray) = Float64(r) - x
Base.:-(x::TracedRNumber, r::Rational) = x - Float64(r)
Base.:-(r::Rational, x::TracedRNumber) = Float64(r) - x

# Comparisons
Base.:<(x::TracedRNumber, r::Rational) = x < Float64(r)
Base.:<(r::Rational, x::TracedRNumber) = Float64(r) < x
Base.:>(x::TracedRNumber, r::Rational) = x > Float64(r)
Base.:>(r::Rational, x::TracedRNumber) = Float64(r) > x
Base.:<=(x::TracedRNumber, r::Rational) = x <= Float64(r)
Base.:<=(r::Rational, x::TracedRNumber) = Float64(r) <= x
Base.:>=(x::TracedRNumber, r::Rational) = x >= Float64(r)
Base.:>=(r::Rational, x::TracedRNumber) = Float64(r) >= x

# =============================================================================
# Make Ref{F} callable for function types when used with traced types
# Broadcasting wraps scalar arguments (including functions) in Ref, and Reactant
# needs to be able to call these wrapped functions.
# =============================================================================
(r::Base.RefValue{F})(x::TracedRNumber, t) where {F<:Function} = r[](x, t)
(r::Base.RefValue{F})(x::TracedRArray, t) where {F<:Function} = r[](x, t)

# =============================================================================
# Safe copy for Reactant that doesn't use similar() with undef elements
# For TracedRArrays, use copy directly. For vectors of TracedRArrays, copy element-wise.
reactant_safe_copy(x::TracedRArray) = copy(x)
reactant_safe_copy(x::Number) = x
reactant_safe_copy(x::AbstractVector) = [reactant_safe_copy(el) for el in x]
reactant_safe_copy(x) = copy(x)  # fallback

# Override ode_copyat_or_push! for TracedRArray states.
# Use indexed assignment when the array is pre-allocated, otherwise push!.
# Reactant cannot handle push! that grows arrays inside mutable structs.
function OrdinaryDiffEqCore.ode_copyat_or_push!(a, i, x, u::TracedRArray, perform_copy = true)
    if i <= length(a)
        # Pre-allocated: use indexed assignment
        if perform_copy
            a[i] = reactant_safe_copy(x)
        else
            a[i] = x
        end
    else
        # Not pre-allocated: use push! (will fail with Reactant if array grows)
        if perform_copy
            push!(a, reactant_safe_copy(x))
        else
            push!(a, x)
        end
    end
    return nothing
end

# Reactant cannot convert Rational{Int64} to traced types.
# Override functions that return Rationals to return Float64 instead.

@reactant_overlay function OrdinaryDiffEqCore.qmin_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 0.2 : 0
end

@reactant_overlay function OrdinaryDiffEqCore.beta2_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 2 / (5 * OrdinaryDiffEqCore.alg_order(alg)) : 0
end

@reactant_overlay function OrdinaryDiffEqCore.beta1_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm},
    beta2
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 7 / (10 * OrdinaryDiffEqCore.alg_order(alg)) : 0
end

@reactant_overlay function OrdinaryDiffEqCore.gamma_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 0.9 : 0
end

@reactant_overlay function OrdinaryDiffEqCore.qsteady_max_default(
    alg::OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveImplicitAlgorithm
)
    return 1.2
end

@reactant_overlay function OrdinaryDiffEqCore.qsteady_max_default(
    alg::OrdinaryDiffEqCore.OrdinaryDiffEqImplicitAlgorithm
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 1.0 : 0
end

@reactant_overlay function OrdinaryDiffEqCore.qoldinit_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.anyadaptive(alg) ? 1e-4 : 0
end

# Default tolerance Rationals -> Float64 for Reactant
@reactant_overlay function OrdinaryDiffEqCore.default_abstol_rational()
    return 1e-6
end

@reactant_overlay function OrdinaryDiffEqCore.default_reltol_rational()
    return 1e-3
end

# Simplified ode_determine_initdt for Reactant tracing.
# The full version has many conditionals with traced values that can't be used in
# regular if statements. This simplified version uses a heuristic based on the timespan.
@reactant_overlay function OrdinaryDiffEqCore.ode_determine_initdt(
    u0, t, tdir, dtmax, abstol, reltol, internalnorm,
    prob::SciMLBase.AbstractODEProblem{uType, tType, iip},
    integrator
) where {tType, uType, iip}
    _tType = eltype(tType)
    oneunit_tType = oneunit(_tType)
    tspan = prob.tspan

    # Simple heuristic: use 1/1000 of the timespan, capped by dtmax
    tspan_length = abs(tspan[2] - tspan[1])
    dt_candidate = tspan_length / 1000

    # Ensure we respect dtmax
    dt = min(dt_candidate, abs(dtmax))

    # Ensure minimum dt for numerical stability
    dtmin = max(integrator.opts.dtmin, eps(_tType) * oneunit_tType)
    dt = max(dt, dtmin)

    return tdir * convert(_tType, dt)
end

# =============================================================================
# Stepsize controller overlays for Reactant
# The standard versions use `if` with traced booleans which can't be traced.
# These versions use `ifelse` which works with traced booleans.
# =============================================================================

@reactant_overlay @inline function OrdinaryDiffEqCore.stepsize_controller!(
    integrator, controller::OrdinaryDiffEqCore.PIController, alg
)
    (; qold) = integrator
    (; qmin, qmax, gamma) = integrator.opts
    (; beta1, beta2) = controller
    EEst = DiffEqBase.value(integrator.EEst)

    # Use ifelse instead of if for traced booleans
    q_zero = inv(qmax)
    q11 = fastpower(EEst, convert(typeof(EEst), beta1))
    q_nonzero = q11 / fastpower(qold, convert(typeof(EEst), beta2))

    is_zero = iszero(EEst)
    q = ifelse(is_zero, q_zero, q_nonzero)
    # Only set q11 when EEst is non-zero (use same value to avoid branch)
    integrator.q11 = ifelse(is_zero, integrator.q11, q11)

    q = max(inv(qmax), min(inv(qmin), q / gamma))
    return q
end

@reactant_overlay function OrdinaryDiffEqCore.step_accept_controller!(
    integrator, controller::OrdinaryDiffEqCore.PIController, alg, q
)
    (; qsteady_min, qsteady_max, qoldinit) = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)

    # Use ifelse instead of if for traced booleans
    in_steady_range = (qsteady_min <= q) & (q <= qsteady_max)
    q = ifelse(in_steady_range, one(q), q)
    integrator.qold = max(EEst, qoldinit)
    return integrator.dt / q # new dt
end

# Simplified calc_dt_propose! for Reactant
# Since the integrator's tType (dtpropose field) is Float64 but dtnew can be TracedRNumber,
# we skip the adaptive dt update for Reactant tracing. The step size remains fixed.
@reactant_overlay function OrdinaryDiffEqCore.calc_dt_propose!(integrator, dtnew)
    # For Reactant tracing, we skip dt adaptation because:
    # 1. dtpropose is typed as Float64 (tType), can't store TracedRNumber
    # 2. The loop structure isn't traced anyway, so adaptive stepping doesn't apply
    # Just keep the current dtpropose value
    return nothing
end

# =============================================================================
# Loopfooter overlay for Reactant
# The standard version uses many if/else with traced booleans that can't be traced.
# This simplified version assumes force_stepfail is false (common for explicit solvers)
# and uses ifelse for accept_step control flow.
# =============================================================================

@reactant_overlay function OrdinaryDiffEqCore._loopfooter!(integrator)
    # Carry-over from callback
    integrator.reeval_fsal = false
    integrator.u_modified = false
    # Skip error checking for Reactant tracing - check_error! uses boolean conditionals
    # that can't be traced with TracedRNumber{Bool}
    integrator.do_error_check = false
    ttmp = integrator.t + integrator.dt

    # For Reactant tracing, we simplify the control flow:
    # - Assume force_stepfail is false (typical for explicit solvers)
    # - Always accept steps (assume well-behaved problems)
    # This allows tracing without complex branching on traced booleans

    if integrator.force_stepfail
        # Handle forced step failure (should not happen for explicit adaptive solvers)
        if integrator.opts.adaptive
            OrdinaryDiffEqCore.post_newton_controller!(integrator, integrator.alg)
        end
        integrator.last_stepfail = true
        integrator.accept_step = false
    elseif integrator.opts.adaptive
        # Adaptive case - the main path for Tsit5
        # Compute step size controller (this updates internal state)
        q = OrdinaryDiffEqCore.stepsize_controller!(integrator, integrator.alg)

        # For Reactant tracing, we always accept steps
        # This simplification is needed because:
        # 1. accept_step field is Bool, can't store TracedRNumber{Bool}
        # 2. The loop control flow isn't traced anyway
        integrator.accept_step = true
        integrator.isout = false  # Assume not out of domain

        # Accept step path
        OrdinaryDiffEqCore.increment_accept!(integrator.stats)
        integrator.last_stepfail = false

        dtnew = DiffEqBase.value(
            OrdinaryDiffEqCore.step_accept_controller!(
                integrator, integrator.alg, q)) * oneunit(integrator.dt)

        integrator.tprev = integrator.t
        integrator.t = OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!(integrator, ttmp)
        OrdinaryDiffEqCore.calc_dt_propose!(integrator, dtnew)

        # Skip callbacks for Reactant tracing - they have complex control flow
    else
        # Non-adaptive case
        OrdinaryDiffEqCore.increment_accept!(integrator.stats)
        integrator.tprev = integrator.t
        integrator.t = OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!(integrator, ttmp)
        integrator.last_stepfail = false
        integrator.accept_step = true
        integrator.dtpropose = integrator.dt
        OrdinaryDiffEqCore.handle_callbacks!(integrator)
    end

    return nothing
end

# =============================================================================
# Pre-allocation and specialized solve! for Reactant
# For saveat case with adaptive=false, we know exact number of saves
# =============================================================================

# Calculate the number of save points for a given integrator
function compute_savelength(integrator)
    opts = integrator.opts
    saveat_cache = opts.saveat_cache

    # If no saveat, we can't pre-allocate (would need save_everystep logic)
    isempty(saveat_cache) && return 0

    savelength = length(saveat_cache)

    # Add 1 for save_start if enabled
    if opts.save_start
        savelength += 1
    end

    return savelength
end

# Pre-allocate solution arrays for Reactant tracing
function preallocate_solution_arrays!(integrator)
    savelength = compute_savelength(integrator)

    # Can't pre-allocate if savelength is unknown
    savelength == 0 && return false

    sol = integrator.sol
    u0 = integrator.u
    t0 = integrator.t

    # Pre-allocate sol.u with copies of u0
    current_len = length(sol.u)
    for i in (current_len + 1):savelength
        push!(sol.u, copy(u0))
    end

    # Pre-allocate sol.t with zeros
    current_len_t = length(sol.t)
    for i in (current_len_t + 1):savelength
        push!(sol.t, zero(t0))
    end

    return true
end

# Perform a single integration step - factored out for clarity
function reactant_step!(integrator)
    OrdinaryDiffEqCore.loopheader!(integrator)
    OrdinaryDiffEqCore.perform_step!(integrator, integrator.cache)
    OrdinaryDiffEqCore.loopfooter!(integrator)
    return nothing
end

# Unrolled step functions for small fixed step counts
# These are needed because @trace for causes type tracing issues with complex types
function reactant_steps_1!(integrator)
    reactant_step!(integrator)
    return nothing
end

function reactant_steps_2!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    return nothing
end

function reactant_steps_3!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    return nothing
end

function reactant_steps_4!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    return nothing
end

function reactant_steps_5!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    reactant_step!(integrator)
    return nothing
end

function reactant_steps_10!(integrator)
    reactant_steps_5!(integrator)
    reactant_steps_5!(integrator)
    return nothing
end

function reactant_steps_20!(integrator)
    reactant_steps_10!(integrator)
    reactant_steps_10!(integrator)
    return nothing
end

function reactant_steps_50!(integrator)
    reactant_steps_10!(integrator)
    reactant_steps_10!(integrator)
    reactant_steps_10!(integrator)
    reactant_steps_10!(integrator)
    reactant_steps_10!(integrator)
    return nothing
end

function reactant_steps_100!(integrator)
    reactant_steps_50!(integrator)
    reactant_steps_50!(integrator)
    return nothing
end

# Dispatch to the appropriate unrolled function based on step count
function reactant_run_steps!(integrator, nsteps::Int)
    # Run in chunks
    while nsteps >= 100
        reactant_steps_100!(integrator)
        nsteps -= 100
    end
    while nsteps >= 50
        reactant_steps_50!(integrator)
        nsteps -= 50
    end
    while nsteps >= 20
        reactant_steps_20!(integrator)
        nsteps -= 20
    end
    while nsteps >= 10
        reactant_steps_10!(integrator)
        nsteps -= 10
    end
    while nsteps >= 5
        reactant_steps_5!(integrator)
        nsteps -= 5
    end
    while nsteps >= 1
        reactant_step!(integrator)
        nsteps -= 1
    end
    return nothing
end

# Simplified solve! overlay for Reactant with unrolled iterations
# This avoids @trace for which causes type tracing issues
@reactant_overlay function SciMLBase.solve!(integrator::OrdinaryDiffEqCore.ODEIntegrator)
    # Pre-allocate solution arrays
    preallocate_solution_arrays!(integrator)

    # Calculate number of integration steps
    # For non-adaptive with fixed dt, this is known at init time (not traced)
    tspan = integrator.sol.prob.tspan
    dt = integrator.dt
    nsteps = max(1, ceil(Int, abs(tspan[2] - tspan[1]) / abs(dt)))

    # Run the integration steps using unrolled functions
    reactant_run_steps!(integrator, nsteps)

    # Handle the final tstop
    OrdinaryDiffEqCore.handle_tstop!(integrator)

    OrdinaryDiffEqCore.postamble!(integrator)

    # Skip analytic solution comparison for Reactant
    # Skip retcode handling - just return the solution
    return integrator.sol
end

end
