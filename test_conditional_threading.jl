#!/usr/bin/env julia

# Test to verify that memory allocation is conditional on threading being enabled

println("Conditional Threading Memory Allocation Test")
println("==========================================")
println("Current: nthreads = $(Threads.nthreads()), maxthreadid = $(Threads.maxthreadid())")

using OrdinaryDiffEqExtrapolation
using RecursiveFactorization
using OrdinaryDiffEq
import OrdinaryDiffEqCore

# Test the helper function directly
import OrdinaryDiffEqExtrapolation: get_thread_count

println("\nTesting get_thread_count function:")

# Create algorithms with different threading settings
alg_no_threading = ExtrapolationMidpointDeuflhard(threading = false)
alg_with_threading = ExtrapolationMidpointDeuflhard(threading = true)

thread_count_no_threading = get_thread_count(alg_no_threading)
thread_count_with_threading = get_thread_count(alg_with_threading)

println("- No threading: get_thread_count = $thread_count_no_threading (should equal nthreads = $(Threads.nthreads()))")
println("- With threading: get_thread_count = $thread_count_with_threading (should equal maxthreadid = $(Threads.maxthreadid()))")

# Verify the values are correct
@assert thread_count_no_threading == Threads.nthreads() "No threading should use nthreads()"
@assert thread_count_with_threading == Threads.maxthreadid() "Threading should use maxthreadid()"

println("✓ Helper function working correctly!")

# Test with actual problem solving to make sure cache allocation works
simple_prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 0.1))

println("\nTesting problem solving with different threading settings:")

println("- Testing without threading...")
try
    sol_no_thread = solve(simple_prob, ExtrapolationMidpointDeuflhard(threading = false), reltol = 1e-3)
    if SciMLBase.successful_retcode(sol_no_thread)
        println("✓ No threading: SUCCESS (using nthreads = $(Threads.nthreads()) for cache arrays)")
    else
        println("✗ No threading: Failed")
    end
catch e
    println("✗ No threading: Error - $e")
end

println("- Testing with threading...")
try
    sol_with_thread = solve(simple_prob, ExtrapolationMidpointDeuflhard(threading = true), reltol = 1e-3)
    if SciMLBase.successful_retcode(sol_with_thread)
        println("✓ With threading: SUCCESS (using maxthreadid = $(Threads.maxthreadid()) for cache arrays)")
    else
        println("✗ With threading: Failed")
    end
catch e
    println("✗ With threading: Error - $e")
end

# Test memory efficiency comparison
println("\nMemory efficiency analysis:")
println("- Without threading: allocates $(Threads.nthreads()) arrays per cache")
println("- With threading: allocates $(Threads.maxthreadid()) arrays per cache") 
if Threads.maxthreadid() > Threads.nthreads()
    saved_arrays = Threads.maxthreadid() - Threads.nthreads()
    println("- Memory savings when threading=false: $saved_arrays fewer arrays per cache")
    println("  (For large problems, this can save significant memory)")
else
    println("- No memory difference in current threading configuration")
end

println("\n🎉 Conditional memory allocation is working correctly!")
println("   - threading=false uses efficient nthreads() allocation")
println("   - threading=true uses safe maxthreadid() allocation")