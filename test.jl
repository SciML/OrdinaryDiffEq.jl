using DiffEqBase, OrdinaryDiffEq
function lorenz!(du,u,p,t)
    du[1] = 10.0*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

u0 = [1.0;0.0;0.0]
du=similar(u0)
tspan = (0.0,100.0)

prob = ODEProblem(lorenz!,u0,tspan)
integR4 = init(prob, Rodas4P())
integR5 = init(prob, Rodas5())

step!(integR4)
step!(integR5)

get_du(integR4).*integR4.dt
integR4.dt
lorenz!(du,integR5.u,0,integR5.t)
all(get_du(integR5) .≈ du) 
lorenz!(du,integR4.u,0,integR4.t)
all(get_du(integR4) .≈ du) 
(du .- get_du(integR4)).^2 |>sum
du
get_du(integR4)
du
get_du(integR5)
integR5.t
integR4.t
integR4(integR4.t,Val{1})
using ForwardDiff
ForwardDiff.derivative(t->integR4(t,Val{0}),integR4.t).-integR4(integR4.t,Val{1})
integ2.cache |>typeof
OrdinaryDiffEq.perform_step!(integ2,integ2.cache)
integ2.fsallast
integ2.cache.dense1
using BenchmarkTools
@benchmark step!(integ2) samples=1     
OrdinaryDiffEq.perform_step!(integ2,integ2.cache) 