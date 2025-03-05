using TaylorDiff: TaylorDiff, extract_derivative, extract_derivative!

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

# DAETS
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
    @unpack Σ, c, d, J, xTS, xtrial, htrial, e, tmp, xcur = cache

    ttrial = t + htrial
    compute_taylor_coefficients!(integrator, cache)
    xTS, ets = sum_taylor_series!(integrator, cache)
    project!(cache.xtrial, cache.xTS, cache.J)
    err = compute_error!(cache.e, cache.xtrial, cache.xTS, ets)
    tol = 1e-6  # Set some tolerance
    if hasproperty(integrator, :opts) && hasproperty(integrator.opts, :reltol)
        tol = integrator.opts.reltol
    end
    
    println("Using tolerance: ", tol)
    println("Current error: ", err)
    if err <= tol
        println("Step accepted")
        integrator.u .= cache.xtrial
        integrator.t = ttrial
        if hasproperty(integrator, :tprev)
            integrator.tprev = t
        end
        if hasproperty(integrator, :tcur)
            integrator.tcur = ttrial
        end
        if hasproperty(integrator, :uprev)
            integrator.uprev .= u
        end
        integrator.dt = cache.htrial
    else
        println("Step rejected, adjusting step size")
        # Temporary step size adjust
        new_htrial = htrial * 0.5 
        cache.htrial = max(new_htrial, 1e-10) 
    end
    if hasproperty(integrator, :stats)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    end

    return nothing
end