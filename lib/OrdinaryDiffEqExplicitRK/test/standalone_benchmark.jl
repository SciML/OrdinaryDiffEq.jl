#!/usr/bin/env julia
#=
Standalone Benchmark: eval_poly_derivative - Old vs New Implementation

Run with:
    julia standalone_benchmark.jl

Requirements:
    Only BenchmarkTools.jl is needed. Install with:
    julia -e 'using Pkg; Pkg.add("BenchmarkTools")'
=#

using BenchmarkTools

println("="^70)
println("Standalone Benchmark: eval_poly_derivative Performance Comparison")
println("="^70)
println()

# ============================================================================
# OLD IMPLEMENTATION (original slow version using prod)
# ============================================================================
"""
Original implementation using `prod` with lambda - causes allocations.
Polynomial format: p(Θ) = coeffs[1]*Θ + coeffs[2]*Θ² + ... + coeffs[n]*Θⁿ
"""
function eval_poly_derivative_OLD(Θ, coeffs, order::Int)
    n = length(coeffs)
    if n <= order
        return zero(eltype(coeffs))
    end
    result = zero(eltype(coeffs))
    for i in (order+1):n
        # SLOW: prod with lambda creates allocations
        coeff = (order == 0) ? coeffs[i] : coeffs[i] * prod(j -> i-j, 0:(order-1))
        Θpow = Θ^(i-order)  # SLOW: computes power separately each iteration
        result += coeff * Θpow
    end
    return result
end

# ============================================================================
# NEW IMPLEMENTATION (optimized with Horner's method + inline falling factorial)
# ============================================================================
"""
Fast inline falling factorial: n × (n-1) × ... × (n-k+1)
No allocations, uses simple loop.
"""
@inline function falling_factorial(n::Integer, k::Integer)
    result = one(n)
    for j in 0:(k-1)
        result *= (n - j)
    end
    return result
end

"""
Optimized implementation using Horner's method.
Polynomial format: p(Θ) = coeffs[1]*Θ + coeffs[2]*Θ² + ... + coeffs[n]*Θⁿ
                        = Θ * (coeffs[1] + Θ*(coeffs[2] + ... ))

- For order=0: Uses Horner's method with muladd (fused multiply-add)
- For order>0: Computes derivative coefficients with falling_factorial, then Horner's
"""
function eval_poly_derivative_NEW(Θ, coeffs, order::Int)
    n = length(coeffs)
    if n <= order
        return zero(eltype(coeffs)) * zero(Θ)
    end

    if order == 0
        # Polynomial: p(Θ) = c[1]*Θ + c[2]*Θ² + ... + c[n]*Θⁿ
        #           = Θ * (c[1] + Θ*(c[2] + Θ*(c[3] + ... + Θ*c[n])))
        # Use Horner's method on the inner part, then multiply by Θ
        result = coeffs[n]
        for i in (n-1):-1:1
            result = muladd(Θ, result, coeffs[i])
        end
        return Θ * result
    else
        # For the k-th derivative of c[i] * Θ^i:
        #   d^k/dΘ^k [c[i] * Θ^i] = c[i] * i*(i-1)*...*(i-k+1) * Θ^(i-k)
        #                        = c[i] * falling_factorial(i, k) * Θ^(i-k)
        #
        # For i from (order+1) to n, the powers Θ^(i-order) range from Θ^1 to Θ^(n-order)
        # So the result is: Θ * (sum of terms with Horner's method)

        result = coeffs[n] * falling_factorial(n, order)
        for i in (n-1):-1:(order+1)
            multiplier = falling_factorial(i, order)
            result = muladd(Θ, result, coeffs[i] * multiplier)
        end
        return Θ * result
    end
end

# ============================================================================
# TEST DATA
# ============================================================================
# Typical Tsit5 interpolation coefficients (5 terms)
const COEFFS_TSIT5 = [0.0, 1.0, -2.763706197274826, 2.9132554618219126, -1.0530884977290216]

# Medium size polynomial (10 terms)
const COEFFS_MEDIUM = Float64.(1:10)

# Large polynomial (20 terms)
const COEFFS_LARGE = Float64.(1:20)

# ============================================================================
# VERIFICATION
# ============================================================================
function verify()
    println("🔍 Verifying correctness...")

    test_cases = [
        (0.5, COEFFS_TSIT5, 0),
        (0.5, COEFFS_TSIT5, 1),
        (0.3, COEFFS_MEDIUM, 0),
        (0.3, COEFFS_MEDIUM, 2),
        (0.7, COEFFS_LARGE, 0),
        (0.7, COEFFS_LARGE, 3),
    ]

    all_passed = true
    for (Θ, coeffs, order) in test_cases
        old = eval_poly_derivative_OLD(Θ, coeffs, order)
        new = eval_poly_derivative_NEW(Θ, coeffs, order)

        if !isapprox(old, new, rtol=1e-10)
            println("   ❌ FAIL: Θ=$Θ, n=$(length(coeffs)), order=$order")
            println("      Old: $old, New: $new")
            all_passed = false
        end
    end

    if all_passed
        println("   ✅ All tests passed - implementations are equivalent\n")
    end
    return all_passed
end

# ============================================================================
# BENCHMARK
# ============================================================================
function run_benchmark()
    println("⏱️  Running benchmarks...\n")
    println("-"^70)

    benchmarks = [
        # (name, Θ, coeffs, order)
        ("Tsit5 (5 terms), order=0  [MOST COMMON]", 0.5, COEFFS_TSIT5, 0),
        ("Tsit5 (5 terms), order=1  [derivative]",  0.5, COEFFS_TSIT5, 1),
        ("Medium (10 terms), order=0",              0.5, COEFFS_MEDIUM, 0),
        ("Medium (10 terms), order=2",              0.5, COEFFS_MEDIUM, 2),
        ("Large (20 terms), order=0",               0.5, COEFFS_LARGE, 0),
        ("Large (20 terms), order=3",               0.5, COEFFS_LARGE, 3),
    ]

    results = []

    for (name, Θ, coeffs, order) in benchmarks
        # Warm-up
        eval_poly_derivative_OLD(Θ, coeffs, order)
        eval_poly_derivative_NEW(Θ, coeffs, order)

        # Benchmark with enough samples for accuracy
        old_b = @benchmark eval_poly_derivative_OLD($Θ, $coeffs, $order)
        new_b = @benchmark eval_poly_derivative_NEW($Θ, $coeffs, $order)

        old_time = median(old_b.times)  # nanoseconds
        new_time = median(new_b.times)
        speedup = old_time / new_time

        old_alloc = old_b.allocs
        new_alloc = new_b.allocs

        push!(results, (; name, old_time, new_time, speedup, old_alloc, new_alloc))

        println("📊 $name")
        println("   OLD: $(round(old_time, digits=1)) ns | $(old_alloc) allocs")
        println("   NEW: $(round(new_time, digits=1)) ns | $(new_alloc) allocs")

        if speedup >= 1.0
            color = speedup >= 2.0 ? "🚀" : "✅"
            println("   $color  $(round(speedup, digits=1))x FASTER")
        else
            println("   ⚠️  $(round(1/speedup, digits=1))x slower")
        end
        println()
    end

    return results
end

# ============================================================================
# SUMMARY
# ============================================================================
function print_summary(results)
    println("="^70)
    println("SUMMARY")
    println("="^70)

    speedups = [r.speedup for r in results]

    println()
    println("📈 Speedup Statistics:")
    println("   Minimum: $(round(minimum(speedups), digits=1))x")
    println("   Average: $(round(sum(speedups)/length(speedups), digits=1))x")
    println("   Maximum: $(round(maximum(speedups), digits=1))x")

    old_allocs = sum(r.old_alloc for r in results)
    new_allocs = sum(r.new_alloc for r in results)

    println()
    println("📦 Allocations:")
    println("   OLD total: $old_allocs allocations")
    println("   NEW total: $new_allocs allocations")

    if new_allocs == 0 && old_allocs > 0
        println("   🎉 ELIMINATED ALL ALLOCATIONS!")
    end

    println()
    println("="^70)
    println()
end

# ============================================================================
# DETAILED ANALYSIS
# ============================================================================
function detailed_analysis()
    println("\n📋 DETAILED ANALYSIS: Why is the new version faster?\n")

    println("1️⃣  ALLOCATION ELIMINATION")
    println("   Old: prod(j -> i-j, 0:(order-1))")
    println("        └─ Creates anonymous function (heap allocation)")
    println("        └─ Creates range object 0:(order-1) (potential allocation)")
    println("        └─ Generic prod() has overhead")
    println()
    println("   New: falling_factorial(n, k) with @inline")
    println("        └─ Simple loop, no allocations")
    println("        └─ Compiler can inline and optimize")
    println()

    println("2️⃣  HORNER'S METHOD (for order=0)")
    println("   Old: Computes Θ^i separately for each term")
    println("        └─ n terms → O(n²) multiplications for powers")
    println()
    println("   New: Uses Horner's scheme")
    println("        └─ p(Θ) = Θ × (c₁ + Θ×(c₂ + Θ×(...)))")
    println("        └─ n terms → O(n) multiplications")
    println()

    println("3️⃣  FUSED MULTIPLY-ADD (muladd)")
    println("   Old: result += coeff * Θpow  (2 operations)")
    println("   New: result = muladd(Θ, result, coeffs[i])  (1 FMA operation)")
    println("        └─ Uses hardware FMA instruction when available")
    println("        └─ More accurate (single rounding)")
    println()
end

# ============================================================================
# MAIN
# ============================================================================
function main()
    if !verify()
        println("❌ Verification failed! Aborting.")
        return
    end

    results = run_benchmark()
    print_summary(results)
    detailed_analysis()

    println("💡 To run this benchmark yourself:")
    println("   julia standalone_benchmark.jl")
    println()
end

main()
