using Test, RecursiveArrayTools, StaticArrays, ForwardDiff

using DiffEqBase: UNITLESS_ABS2, recursive_length, ODE_DEFAULT_NORM

@test recursive_length(1.0) == 1

n = UNITLESS_ABS2(3.0 + 4.0im)
@test n == 25.0
@test typeof(n) <: Real

@test ODE_DEFAULT_NORM(3.0 + 4.0im, 0.0) == 5.0

u1 = ones(3)
@test UNITLESS_ABS2(u1) == 3.0
@test recursive_length(u1) == 3
@test ODE_DEFAULT_NORM(u1, 0.0) == 1.0

u2 = [SA[1.0 1.0; 1.0 1.0] for i in 1:3]
@test UNITLESS_ABS2(u2) == 12.0
@test recursive_length(u2) == 12
@test ODE_DEFAULT_NORM(u2, 0.0) == 1.0

u3 = VectorOfArray([ones(5), ones(5)])
@test UNITLESS_ABS2(u3) == 10.0
@test recursive_length(u3) == 10
@test ODE_DEFAULT_NORM(u3, 0.0) == 1.0

u4 = ArrayPartition(u1, u2, u3)
@test UNITLESS_ABS2(u4) == 25.0
@test recursive_length(u4) == 25
@test ODE_DEFAULT_NORM(u4, 0.0) == 1.0

u5 = ArrayPartition(u4, u4)
@test UNITLESS_ABS2(u5) == 50.0
@test recursive_length(u5) == 50
@test ODE_DEFAULT_NORM(u5, 0.0) == 1.0

u6 = ArrayPartition(1.0, 1.0)
@test UNITLESS_ABS2(u6) == 2.0
@test recursive_length(u6) == 2
@test ODE_DEFAULT_NORM(u6, 0.0) == 1.0

u7 = ArrayPartition(u1, ones(0))
@test UNITLESS_ABS2(u7) == 3.0
@test recursive_length(u7) == 3
@test ODE_DEFAULT_NORM(u7, 0.0) == 1.0

@test ODE_DEFAULT_NORM(Float64[], 0.0) == 0.0

# https://github.com/SciML/DiffEqBase.jl/issues/1023
u8 = ForwardDiff.Dual{:b}.(ForwardDiff.Dual{:a}.([1.0, 2.0, 3.0], true), true)
u8_ref = 1.2909944487358056
@test ODE_DEFAULT_NORM(u8, 4.0) isa Float64
@test ODE_DEFAULT_NORM(u8, 4.0) ≈ u8_ref
@test ODE_DEFAULT_NORM(u8, ForwardDiff.Dual{:b}(4.0, true)) isa Float64
@test ODE_DEFAULT_NORM(u8, ForwardDiff.Dual{:b}(4.0, true)) ≈ u8_ref
