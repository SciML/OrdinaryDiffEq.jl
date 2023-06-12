using OrdinaryDiffEq, LinearAlgebra, Test

function dynamics!(dx, x, θ, t)
    zmin, zmax, M, g = θ
    z, dz, p = x

    dx[1] = dz
    dx[2] = -M * g - z + p
    dx[3] = p - w(t) # the last state (algebraic) p has to be equal to w

    if z <= zmin
        dx[2] = max(dx[2], 0.0)
    elseif z >= zmax
        dx[2] = min(dx[2], 0.0)
    end
end

w(t) = 0.1 * (7.0 + 3.0sin(t))
zmin_cond(x, t, integrator) = x[1] - zmin
zmax_cond(x, t, integrator) = zmax - x[1]

function zmin_affect_neg!(integrator)
    integrator.u[1] = zmin
    integrator.u[2] = 0.0
end

function zmax_affect_neg!(integrator)
    integrator.u[1] = zmax
    integrator.u[2] = 0.0
end

cbs = CallbackSet(ContinuousCallback(zmin_cond, zmin_affect_neg!),
    ContinuousCallback(zmax_cond, zmax_affect_neg!))

tf = 20.0
tspan = (0.0, tf)

zmin = 1.0e-3
zmax = 2.0e-2
M = 0.05
g = 9.81
θ = zmin, zmax, M, g

x0 = [(zmin + zmax) / 2, 0.0, w(0.0)]

E = diagm([1.0, M, 0.0])

f = ODEFunction(dynamics!, mass_matrix = E)
prob = ODEProblem(f, x0, tspan, θ)

sol1 = solve(prob, Rodas4(), callback = cbs, reltol = 1e-6)
@test sol1(0.06692341688237893)[3]≈0.72 atol=1e-2
sol1 = solve(prob, Rodas5(), callback = cbs, reltol = 1e-6)
@test sol1(0.06692341688237893)[3]≈0.72 atol=1e-2
sol1 = solve(prob, Rodas5P(), callback = cbs, reltol = 1e-6)
@test sol1(0.06692341688237893)[3]≈0.72 atol=1e-2

#=
sol1 = solve(prob,Rosenbrock23(),callback=cbs, reltol=1e-6)
@test sol1(1.0)[3] ≈ 0.95 atol=1e-2
=#
