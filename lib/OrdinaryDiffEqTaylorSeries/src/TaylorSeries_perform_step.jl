using TaylorDiff: TaylorDiff, extract_derivative, extract_derivative!
using LinearAlgebra
using NLsolve

@inline make_taylor(all::Vararg{X, P}) where {P, X <: AbstractArray} = TaylorArray(Base.first(all), Base.tail(all))
@inline make_taylor(all::Vararg{X, P}) where {P, X} = TaylorScalar(all)

function initialize!(integrator, cache::ExplicitTaylor2ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    k1 = f(uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    k2 = f(u1, p, t1).partials[1]
    u = @.. uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylor2Cache)
    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k1, k2, k3, utilde, tmp = cache

    # The following code is written to be fully non-allocating
    f(k1, uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    out1 = make_taylor(k1, k2)
    f(out1, u1, p, t1)
    @.. u = uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end

function initialize!(integrator, cache::ExplicitTaylorConstantCache{P}) where P
    integrator.kshortsize = P
    integrator.k = typeof(integrator.k)(undef, P)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylorConstantCache{P}, repeat_step = false) where P
    @unpack t, dt, uprev, u, f, p = integrator
    us = typeof(u)[]
    integrator.k[1] = f(uprev, p, t)
    push!(us, integrator.k[1])
    u = @.. uprev + dt * us[1]
    dti = dt
    for i in 1:P-1
        ui = make_taylor(uprev, us...)
        ti = TaylorScalar{i}(t, one(t))
        integrator.k[i + 1] = f(ui, p, ti).partials[i]
        push!(us, integrator.k[i + 1] / (i + 1))
        dti *= dt
        u += dti * us[i + 1]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P)
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylorCache{P}) where P
    integrator.kshortsize = P
    resize!(integrator.k, P)
    # Setup k pointers
    for (i, k) in enumerate(cache.ks)
        integrator.k[i] = k
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylorCache{P}, repeat_step = false) where P
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack ks, us = cache

    # The following code is written to be fully non-allocating
    f(ks[1], uprev, p, t)
    @.. us[1] .= ks[1]
    @.. u = uprev + dt * us[1]
    dti = dt
    for i in 1:P-1
        ui = make_taylor(uprev, us[1:i]...)
        ti = TaylorScalar{i}(t, one(t))
        outi = make_taylor(ks[1:i+1]...)
        f(outi, ui, p, ti)
        us[i + 1] .= ks[i + 1] / (i + 1)
        dti *= dt
        @.. u += dti * us[i + 1]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end

# ------------------------------------------------------------------------------
# Differential Algebriac Equation Taylor Series
# ------------------------------------------------------------------------------

function initialize!(integrator, cache::DAETSCache)
    resize!(cache.xTS, length(integrator.u))
    fill!(cache.xTS, 0.0)
    cache.xtrial = zero(integrator.u)  # Projected solution
    cache.htrial = integrator.dt  #Step size
    cache.e = 0.0  # Error estimate
    cache.tmp = zero(integrator.u)  # Temp storage for intermediate calculations
    cache.atmp = zero(integrator.u)  # Temp storage for error compute
    f = integrator.f 
    vars = integrator.p 
    t = integrator.t 

    #Preprocessing
    cache.Σ = signature_matrix(f, vars, t)
    T, _ = highest_value_transversal(cache.Σ)
    cache.c, cache.d = find_offsets(cache.Σ, T)
    cache.J = system_jacobian(f, vars, t, cache.c, cache.d, cache.Σ)
    return nothing
end

# Below are more helper functions for DAETS. Couldn't get these to work with the tests in a different file so I put them here. Probably should be moved to DAETS_utils.jl somehow.
function compute_taylor_coefficients!(integrator, cache)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack Σ, c, d, J, xTS, xtrial, htrial, e, tmp = cache

    order = 3  # Set some default order
    if hasproperty(integrator.alg, :order)
        order = integrator.alg.order
    end

    # Debug information
    println("Computing Taylor coefficients for order p = ", order)
    println("Current time t = ", t)
    println("Current step size dt = ", dt)

    # Initialize array for Taylor coefficients
    n = length(u)
    xcur = zeros(n)  # First-order coefficients

    if f isa ODEFunction
        f(tmp, uprev, p, t)  # Evaluate f at current state
        xcur .= tmp          # First coefficient is the derivative
        #########################################################
        # TODO: Add higher-order coefficients ###################
        #########################################################
        # Store the coefficients in the cache
        cache.xcur = xcur
        
        println("Computed first-order coefficients: ", xcur)
        return xcur
    else
        try
            for i in 1:n
                xcur[i] = TaylorDiff.extract_derivative(f, u, t, dt, i)
                println("Coefficient for variable ", i, ": ", xcur[i])
            end
            cache.xcur = xcur
            return xcur
        catch e
            # Use finite differences?
            println("Warning: Using finite differences for derivatives")
            f(tmp, uprev, p, t)  # Evaluate f at current state
            xcur .= tmp          # First coefficient is the derivative
            cache.xcur = xcur
            return xcur
        end
    end
end

function sum_taylor_series!(integrator, cache)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack Σ, c, d, J, xTS, xtrial, htrial, e, tmp, xcur = cache

    order = 3
    if hasproperty(integrator.alg, :order)
        order = integrator.alg.order
    end

    println("Summing Taylor series for order p = ", order)
    println("Trial step size htrial = ", htrial)

    # Evaluate the Taylor series at trial time
    n = length(u)
    for i in 1:n
        xTS[i] = 0.0
    end
    ets = 0.0  # Error estimate

    if xcur isa Vector{<:Number}
        # First-order Taylor approximation: x(t+h) ≈ x(t) + h*x'(t)
        for i in 1:n
            xTS[i] = uprev[i] + htrial * xcur[i]
        end
        # Error estimate based on truncation error
        ets = norm(xcur) * htrial^2 / 2
        println("Using first-order Taylor approximation")
    else
        # For higher-order Taylor series (when xcur contains higher derivatives)
        # This would be used if xcur was a vector of Taylor series or similar
        println("Warning: Higher-order Taylor series not fully implemented")
        # Fallback to first-order approximation
        for i in 1:n
            xTS[i] = uprev[i] + htrial * xcur[i]
        end
        ets = norm(xcur) * htrial^2 / 2
    end

    println("Taylor series approximation: ", xTS)
    println("Error estimate: ", ets)
    
    # Update the cache
    cache.xTS .= xTS
    
    return xTS, ets
end

using LinearAlgebra

function project!(xtrial, xTS, J)
    println("Projecting solution onto constraints")
    println("Unprojected solution xTS = ", xTS)
    println("System Jacobian J = ", J)

    #########################################################
    # TODO: This is a placeholder ###########################
    #########################################################
    xtrial .= xTS
    
    # In a real implementation we need something like:
    # try
    #     xtrial .= J \ xTS
    # catch
    #     # Fallback if the system is singular
    #     xtrial .= xTS
    # end
    
    println("Projected solution xtrial = ", xtrial)
    return xtrial
end

function compute_error!(e_out, xtrial, xTS, ets)
    e = ets + norm(xtrial - xTS)  # Combine Taylor series error and projection error
    e_out = e
    println("Combined error estimate: ", e)
    return e
end

@muladd function perform_step!(integrator, cache::DAETSCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k1, k2, k3, utilde, tmp, diff_indices, jacobian, taylor_coeffs, error_estimate, 
            u_unprojected, residual, correction = cache

    # ===== Stepping Algorithm =====
    
    # Step 1: Setup current step parameters
    tcur = t + dt
    htrial = dt
    
    # Limit the order
    order = min(integrator.opts.maxiters, 20)
    
    # Step 2: ComputeTCs - Compute Taylor coefficients to required order
    duprev = tmp
    if integrator.iter == 0 && !isempty(integrator.sol.duals)
        duprev .= integrator.sol.duals[1]
    else
        # Use NLsolve for consistent derivatives
        function F!(res, x)
            f(res, x, uprev, p, t)
        end

        result = nlsolve(F!, zeros(length(uprev)), 
                         method=:newton, autodiff=:forward,
                         ftol=1e-10, iterations=10)

        if converged(result)
            duprev .= result.zero
        else
            error("Failed to find consistent initial derivatives.")
        end
    end
    
    # Might need to update the Jacobian here...
    compute_taylor_coefficients!(taylor_coeffs, f, uprev, duprev, p, t, htrial, order, jacobian)
    
    # Step 3: SumTS - Compute unprojected Taylor series solution
    evaluate_taylor_series!(u_unprojected, taylor_coeffs, htrial, order)
    
    # Step 4: Project - Project the Taylor series solution onto the constraints
    projection_success = project_solution!(u, u_unprojected, f, p, tcur, jacobian)
    
    # Check for projection
    if !projection_success
        @warn "Projection failure in DAE solution"
        integrator.dt = 0.5 * dt  # Reduce step size and try again
        return nothing
    end
    
    # Step 5: Compute error estimate
    compute_error_estimate!(error_estimate, u, u_unprojected, taylor_coeffs, order, htrial)
    
    # Add check to see if we should accept step...? TODO: Check if this is already done by the integrator.

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, order+1)
    return nothing
end

## Helper Functions
function compute_taylor_coefficients!(taylor_coeffs, f, u0, du0, p, t, h, order, jacobian)
    if length(taylor_coeffs) < order + 1
        resize!(taylor_coeffs, order + 1)
        for i in eachindex(taylor_coeffs)
            if i > length(taylor_coeffs) || taylor_coeffs[i] === nothing
                taylor_coeffs[i] = similar(u0)
            end
        end
    end
    taylor_coeffs[1] .= u0
    taylor_coeffs[2] .= du0
    
    n = length(u0)
    
    # Extract differential vs. algebraic information from the problem structure
    # This should ideally come from the cache... : TODO: Add preprocessing to cache. Check how DAEProblem does it.
    diff_indices = zeros(Int, n)
    
    # Try to extract differential/algebraic indices from the problem
    if hasfield(typeof(f), :sys) && applicable(hasfield, typeof(f.sys), :differential_vars)
        differential_vars = f.sys.differential_vars
        for i in 1:n
            diff_indices[i] = differential_vars[i] ? 0 : 1
        end
    elseif hasfield(typeof(p), :differential_vars)
        differential_vars = p.differential_vars
        for i in 1:n
            diff_indices[i] = differential_vars[i] ? 0 : 1
        end
    elseif hasfield(typeof(f), :differential_vars)
        differential_vars = f.differential_vars
        for i in 1:n
            diff_indices[i] = differential_vars[i] ? 0 : 1
        end
    else
        for i in 1:n
            diff_indices[i] = i <= div(n, 2) ? 0 : 1
        end
    end
    
    # Use TaylorDiff for higher-order coefficients
    for k in 3:order+1
        u_taylor_values = [taylor_coeffs[i] for i in 1:k-1]
        u_taylor = make_taylor(u_taylor_values...)
        t_taylor = TaylorScalar{k-2}(t, one(t))
        solve_for_coefficient!(taylor_coeffs[k], jacobian, f, u_taylor, p, t_taylor, k-1, diff_indices)
        taylor_coeffs[k] ./= factorial(k-1)
    end
end

function solve_for_coefficient!(coeff, jacobian, f, u_taylor, p, t_taylor, order, diff_indices)
    n = length(coeff)
    rhs = zeros(n)

    du_taylor_values = Vector{typeof(u_taylor.value)}(undef, length(u_taylor.partials) + 1)
    du_taylor_values[1] = zeros(n)
    
    # Higher-order terms are adjusted based on differential/algebraic structure
    for i in 1:length(u_taylor.partials)
        du_taylor_values[i+1] = similar(u_taylor.value)
        for j in 1:n
            if diff_indices[j] == 0  # Diff eqn.
                if i < length(u_taylor.partials)
                    du_taylor_values[i+1][j] = u_taylor.partials[i][j] * factorial(i)
                else
                    du_taylor_values[i+1][j] = 0.0
                end
            else # Alg eqn.
                # For algebraic constraints, higher derivatives are 0
                du_taylor_values[i+1][j] = 0.0
            end
        end
    end
    
    du_taylor = make_taylor(du_taylor_values...)
    res = similar(u_taylor.value)
    # Use TaylorDiff to automatically compute all needed derivatives. TODO: Check if we need to do this in place or out of place.
    if applicable(f, res, du_taylor, u_taylor, p, t_taylor)
        # In-place
        f(res, du_taylor, u_taylor, p, t_taylor)
    else
        # Out-of-place
        res = f(du_taylor, u_taylor, p, t_taylor)
    end
    
    # Extract the appropriate derivatives for the RHS of our linear system
    if hasfield(typeof(res), :partials) && !isempty(res.partials)
        for i in 1:n
            if diff_indices[i] == 0  # Diff eqn.     
                # For differential equations, use the order-th derivative
                if order < length(res.partials)
                    rhs[i] = -res.partials[order][i]
                end
            else  # Alg eqn.
                # For algebraic constraints, use -1 order
                if order-1 < length(res.partials)
                    rhs[i] = -res.partials[order-1][i]
                end
            end
        end
    else
        @warn "No partials found in result. Cannot compute higher-order coefficients."
        fill!(coeff, 0.0)
        return
    end

    # Solve the linear system using the provided Jacobian
    if isnothing(jacobian) || size(jacobian, 1) != n
        @warn "No valid Jacobian available for coefficient calculation."
        fill!(coeff, 0.0)
        return
    end
    
    try
        # J * coeff = rhs
        ldiv!(coeff, lu(jacobian), rhs)
    catch e
        error("Failed to solve for coefficient: $(e)")
    end
end

function evaluate_taylor_series!(u_next, taylor_coeffs, h, order)
    u_next .= taylor_coeffs[1]
    for i in 1:min(order, length(taylor_coeffs)-1)
        u_next .+= h^i .* taylor_coeffs[i+1]
    end
end

function project_solution!(u_projected, u_unprojected, f, p, t, jacobian)
    u_projected .= u_unprojected  # Start with unprojected solution as initial guess
    
    # Define residual function for NLsolve
    function F!(res, x)
        f(res, zeros(length(x)), x, p, t)  # evaluate with zero derivatives
    end
    
    # Use NLsolve with Newton's method
    result = nlsolve(F!, u_unprojected,  # Start from unprojected solution
                    method=:newton,
                    autodiff=:forward,
                    ftol=1e-10,
                    iterations=10)
    
    if converged(result)
        u_projected .= result.zero
        return true
    else
        # If NLsolve fails, projection failed
        @warn "NLsolve failed to project solution onto constraint manifold"
        return false
    end
end

# Compute error estimate
function compute_error_estimate!(error_estimate, u_projected, u_unprojected, taylor_coeffs, order, dt)
    max_order = length(taylor_coeffs) - 1
    order = min(order, max_order)
    
    truncation_error = 0.0
    if order+1 <= length(taylor_coeffs)
        truncation_error = norm(taylor_coeffs[order+1] * dt^order / factorial(order))
    end
    projection_error = norm(u_projected - u_unprojected)
    error_estimate = truncation_error + projection_error
    return error_estimate
end