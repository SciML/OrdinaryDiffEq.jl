using DiffEqDevTools, Test
using OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK: RK4
using OrdinaryDiffEqSDIRK: ImplicitEuler
using OrdinaryDiffEqSSPRK: SSPRK33, SSPRK104

# Test stability_region with algorithm-based interface
@test @inferred(stability_region(Tsit5())) ≈ stability_region(Tsit5()) rtol = 1.0e-3
@test @inferred(stability_region(SSPRK104())) ≈ stability_region(SSPRK104()) rtol = 1.0e-3
@test @inferred(stability_region(ImplicitEuler())) ≈ stability_region(ImplicitEuler()) rtol = 1.0e-3
# Use -1e16 rather than nextfloat(typemin(Float64)): the extreme value produces subnormal
# Newton corrections (1/|z| < floatmin) that FTZ-enabled BLAS flushes to zero, causing the
# Newton solve to fail and the integrator to return u_old = 1 instead of the correct ~0.
# With z = -1e16 the correction 1e-16 is a normal float and the solve converges robustly.
@test @inferred(abs(stability_region(-1.0e16, ImplicitEuler()))) < eps(Float64)

@test @inferred(imaginary_stability_interval(SSPRK33())) ≈ sqrt(3)
@test @inferred(imaginary_stability_interval(SSPRK33(), initial_guess = 5.0f0)) ≈ sqrt(3.0f0)
@test @inferred(imaginary_stability_interval(RK4())) ≈ 2.8284271247
