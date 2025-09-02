#!/usr/bin/env julia

# Simple test to verify our threading fix for issue #2612

println("Running on $(Threads.nthreads()) thread(s), maxthreadid = $(Threads.maxthreadid())")

# Directly test that our array sizes are correct
using OrdinaryDiffEqExtrapolation
using RecursiveFactorization

# Test: Create a simple ODE problem and solve with threading enabled
using OrdinaryDiffEq
import OrdinaryDiffEqCore

# Simple exponential growth problem
function simple_ode(u, p, t)
    return u
end

# Create problem
prob = ODEProblem(simple_ode, 1.0, (0.0, 1.0))

println("\nTesting ExtrapolationMidpointDeuflhard with threading...")
try
    sol = solve(prob, ExtrapolationMidpointDeuflhard(threading = true), reltol = 1e-3, abstol = 1e-6)
    if sol.retcode == :Success
        println("✓ ExtrapolationMidpointDeuflhard with threading: SUCCESS")
        println("  Final value: $(sol.u[end]) (expected ≈ $(exp(1.0)))")
    else
        println("✗ ExtrapolationMidpointDeuflhard failed with retcode: $(sol.retcode)")
    end
catch e
    println("✗ ExtrapolationMidpointDeuflhard with threading failed with error: $e")
end

println("\nTesting ImplicitEulerExtrapolation with threading...")
try
    sol = solve(prob, ImplicitEulerExtrapolation(threading = true), reltol = 1e-3, abstol = 1e-6)
    if sol.retcode == :Success
        println("✓ ImplicitEulerExtrapolation with threading: SUCCESS")
        println("  Final value: $(sol.u[end]) (expected ≈ $(exp(1.0)))")
    else
        println("✗ ImplicitEulerExtrapolation failed with retcode: $(sol.retcode)")
    end
catch e
    println("✗ ImplicitEulerExtrapolation with threading failed with error: $e")
end

println("\nTesting ImplicitDeuflhardExtrapolation with threading...")
try
    sol = solve(prob, ImplicitDeuflhardExtrapolation(threading = true), reltol = 1e-3, abstol = 1e-6)
    if sol.retcode == :Success
        println("✓ ImplicitDeuflhardExtrapolation with threading: SUCCESS")
        println("  Final value: $(sol.u[end]) (expected ≈ $(exp(1.0)))")
    else
        println("✗ ImplicitDeuflhardExtrapolation failed with retcode: $(sol.retcode)")
    end
catch e
    println("✗ ImplicitDeuflhardExtrapolation with threading failed with error: $e")
end

println("\n✅ All threading tests completed!")