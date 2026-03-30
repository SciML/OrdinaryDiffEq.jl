using DiffEqDevTools, Test
using OrdinaryDiffEq

@test stability_region(constructDormandPrince6(), initial_guess = -3.5)≈-3.95413 rtol=1e-3
@test stability_region(constructTsitourasPapakostas6(), initial_guess = -3.5)≈-3.95413 rtol=1e-3
@test stability_region(constructRadauIIA5(), initial_guess = 12.0)≈11.84 rtol=1e-2

@test @inferred(stability_region(constructTsitouras5()))≈@inferred(stability_region(Tsit5()))
@test @inferred(stability_region(constructSSPRK104()))≈@inferred(stability_region(SSPRK104()))
@test @inferred(stability_region(constructImplicitEuler()))≈@inferred(stability_region(ImplicitEuler()))
@test @inferred(abs(stability_region(nextfloat(typemin(Float64)), ImplicitEuler()))) < eps(Float64)

@test @inferred(imaginary_stability_interval(constructSSPRK33()))≈sqrt(3)
@test @inferred(imaginary_stability_interval(constructSSPRK33(Float32)))≈sqrt(3.0f0)
@test @inferred(imaginary_stability_interval(constructKutta3()))≈sqrt(3)
@test @inferred(imaginary_stability_interval(constructKutta3(Float32)))≈sqrt(3.0f0)
@test @inferred(imaginary_stability_interval(constructRK4()))≈2.8284271247

@test @inferred(imaginary_stability_interval(SSPRK33()))≈sqrt(3)
@test @inferred(imaginary_stability_interval(SSPRK33(), initial_guess = 5.0f0))≈sqrt(3.0f0)
@test @inferred(imaginary_stability_interval(RK4()))≈2.8284271247
