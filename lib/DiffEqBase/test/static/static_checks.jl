using DiffEqBase, ComponentArrays, AllocCheck, Test

u = ComponentArray(x = 1.0, y = 0.0, z = 0.0)
t = 0.0
@test length(check_allocs(DiffEqBase.ODE_DEFAULT_NORM, (typeof(u), typeof(t)))) == 0
