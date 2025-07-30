"""
GPU-safe automatic differentiation wrapper for OrdinaryDiffEq.

This implements the solution suggested in https://github.com/JuliaDiff/DifferentiationInterface.jl/issues/820
where we automatically wrap AutoForwardDiff with AutoForwardFromPrimitive when dealing with 
arrays that don't support fast scalar indexing (like GPU arrays).
"""

using ArrayInterface
using ADTypes
using DifferentiationInterface

"""
    make_gpu_safe_autodiff(backend, u)

Automatically wrap the AD backend with AutoForwardFromPrimitive if the array type `u` 
doesn't support fast scalar indexing (e.g., GPU arrays).

# Arguments
- `backend`: The original AD backend (e.g., AutoForwardDiff())
- `u`: The array to differentiate with respect to

# Returns
The original backend if fast scalar indexing is supported, otherwise wrapped with AutoForwardFromPrimitive.
"""
function make_gpu_safe_autodiff(backend::AutoForwardDiff, u)
    if ArrayInterface.fast_scalar_indexing(u)
        # CPU arrays - use original backend
        return backend
    else
        # GPU arrays - wrap with AutoForwardFromPrimitive for safety
        println("GPU array detected (fast_scalar_indexing=false), using AutoForwardFromPrimitive wrapper")
        return DifferentiationInterface.AutoForwardFromPrimitive(backend)
    end
end

# Fallback for other backends - return as-is
function make_gpu_safe_autodiff(backend, u)
    return backend
end

"""
Test function to demonstrate the GPU-safe wrapper behavior.
"""
function test_gpu_safe_autodiff()
    println("=== Testing GPU-Safe AutoDiff Wrapper ===")
    
    # Test with CPU array (regular Array)
    cpu_array = [1.0, 2.0, 3.0]
    println("CPU array fast_scalar_indexing: ", ArrayInterface.fast_scalar_indexing(cpu_array))
    cpu_backend = make_gpu_safe_autodiff(AutoForwardDiff(), cpu_array)
    println("CPU backend: ", typeof(cpu_backend))
    
    # Test with GPU array (JLArray)  
    gpu_array = jl([1.0, 2.0, 3.0])
    println("GPU array fast_scalar_indexing: ", ArrayInterface.fast_scalar_indexing(gpu_array))
    gpu_backend = make_gpu_safe_autodiff(AutoForwardDiff(), gpu_array)
    println("GPU backend: ", typeof(gpu_backend))
    
    println()
    return cpu_backend, gpu_backend
end

# Export the main function
export make_gpu_safe_autodiff