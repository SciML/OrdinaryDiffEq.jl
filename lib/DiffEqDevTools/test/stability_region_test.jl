using DiffEqDevTools, Test
using OrdinaryDiffEq
# v7 OrdinaryDiffEq does not blanket-reexport every algorithm; pull in the
# specific solvers this file uses from their sublibraries.
using OrdinaryDiffEqLowOrderRK: RK4
using OrdinaryDiffEqSDIRK: ImplicitEuler
using OrdinaryDiffEqSSPRK: SSPRK33, SSPRK104

# Test stability_region with algorithm-based interface
@test @inferred(stability_region(Tsit5())) ≈ stability_region(Tsit5()) rtol = 1.0e-3
@test @inferred(stability_region(SSPRK104())) ≈ stability_region(SSPRK104()) rtol = 1.0e-3
@test @inferred(stability_region(ImplicitEuler())) ≈ stability_region(ImplicitEuler()) rtol = 1.0e-3
@test @inferred(abs(stability_region(nextfloat(typemin(Float64)), ImplicitEuler()))) < eps(Float64)

@test @inferred(imaginary_stability_interval(SSPRK33())) ≈ sqrt(3)
@test @inferred(imaginary_stability_interval(SSPRK33(), initial_guess = 5.0f0)) ≈ sqrt(3.0f0)
@test @inferred(imaginary_stability_interval(RK4())) ≈ 2.8284271247
