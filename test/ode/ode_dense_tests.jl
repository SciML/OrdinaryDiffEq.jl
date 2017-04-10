using OrdinaryDiffEq, DiffEqProblemLibrary,Base.Test, DiffEqBase

bools = Vector{Bool}(0)
prob = prob_ode_linear

sol =solve(prob,Euler(),dt=1//2^(2),dense=true)

interpd_1d = sol(0:1//2^(4):1)

sol2 =solve(prob,Euler(),dt=1//2^(4),dense=true)

sol3 =solve(prob,Euler(),dt=1//2^(5),dense=true)

prob = prob_ode_2Dlinear
sol =solve(prob,Euler(),dt=1//2^(2),dense=true)

interpd = sol(0:1//2^(4):1)

interpd_idxs = sol(0:1//2^(4):1,idxs=1:2:5)

@test minimum([interpd_idxs[i] == interpd[i][1:2:5] for i in 1:length(interpd)])

interpd_single = sol(0:1//2^(4):1,idxs=1)

@test typeof(interpd_single) <: Vector{Float64}

@test typeof(sol(0.5,idxs=1)) <: Float64

A = rand(4,2)
sol(A,0.777)
A == sol(0.777)
sol2 =solve(prob,Euler(),dt=1//2^(4),dense=true)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < .2

sol(interpd,0:1//2^(4):1)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < .2


sol =solve(prob,Euler(),dt=1//2^(2),dense=false)

@test_throws ErrorException sol(0.5)

prob = prob_ode_linear

sol =solve(prob,Midpoint(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Midpoint(),dt=1//2^(4),dense=true)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 2e-2

prob = prob_ode_2Dlinear

sol =solve(prob,Midpoint(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Midpoint(),dt=1//2^(4),dense=true)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 2.3e-2

# SSPRK22
prob = prob_ode_linear
sol =solve(prob,SSPRK22(),dt=1//2^(2),dense=true)
sol(interpd_1d,0:1//2^(4):1)
sol2 =solve(prob,SSPRK22(),dt=1//2^(4),dense=true)
@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1.5e-2

prob = prob_ode_2Dlinear
sol =solve(prob,SSPRK22(),dt=1//2^(2),dense=true)
sol(interpd,0:1//2^(4):1)
sol2 =solve(prob,SSPRK22(),dt=1//2^(4),dense=true)
@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 2.5e-2

# SSPRK33
prob = prob_ode_linear
sol =solve(prob,SSPRK33(),dt=1//2^(2),dense=true)
sol(interpd_1d,0:1//2^(4):1)
sol2 =solve(prob,SSPRK33(),dt=1//2^(4),dense=true)
@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 7.5e-4

prob = prob_ode_2Dlinear
sol =solve(prob,SSPRK33(),dt=1//2^(2),dense=true)
sol(interpd,0:1//2^(4):1)
sol2 =solve(prob,SSPRK33(),dt=1//2^(4),dense=true)
@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1.5e-3

# SSPRK104
prob = prob_ode_linear
sol =solve(prob,SSPRK104(),dt=1//2^(2),dense=true)
sol(interpd_1d,0:1//2^(4):1)
sol2 =solve(prob,SSPRK104(),dt=1//2^(4),dense=true)
@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1.5e-5

prob = prob_ode_2Dlinear
sol =solve(prob,SSPRK104(),dt=1//2^(2),dense=true)
sol(interpd,0:1//2^(4):1)
sol2 =solve(prob,SSPRK104(),dt=1//2^(4),dense=true)
@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 3.e-5

# RK4
prob = prob_ode_linear

sol =solve(prob,RK4(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,RK4(),dt=1//2^(4),dense=true)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 5e-5

prob = prob_ode_2Dlinear

sol =solve(prob,RK4(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,RK4(),dt=1//2^(4),dense=true)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-4

prob = prob_ode_linear

sol =solve(prob,DP5(),dense=true)

sol2 =solve(prob,DP5(),dt=1//2^(4),dense=true,adaptive=false)

sol(interpd_1d,0:1//2^(4):1)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 5e-6

prob = prob_ode_2Dlinear

sol =solve(prob,DP5(),dense=true)

sol2 =solve(prob,DP5(),dt=1//2^(4),dense=true,adaptive=false)

sol(interpd,0:1//2^(4):1)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-5

prob = prob_ode_linear

sol =solve(prob,BS3(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,BS3(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1e-3

prob = prob_ode_2Dlinear

sol =solve(prob,BS3(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,BS3(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-3


prob = prob_ode_linear

sol =solve(prob,Tsit5(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Tsit5(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1e-5

prob = prob_ode_2Dlinear

sol =solve(prob,Tsit5(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Tsit5(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-5

prob = prob_ode_linear

sol =solve(prob,TanYam7(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,TanYam7(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1e-3

prob = prob_ode_2Dlinear

sol =solve(prob,TanYam7(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,TanYam7(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-3


prob = prob_ode_linear

sol =solve(prob,TsitPap8(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,TsitPap8(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 2e-2

prob = prob_ode_2Dlinear

sol =solve(prob,TsitPap8(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,TsitPap8(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-2


prob = prob_ode_linear

sol =solve(prob,Feagin10(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Feagin10(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1e-3

prob = prob_ode_2Dlinear

sol =solve(prob,Feagin10(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Feagin10(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-3


### Vern6()
const linear_bigα = parse(BigFloat,"1.01")
f = (t,u) -> (linear_bigα*u)
prob_ode_bigfloatlinear = ODEProblem(f,parse(BigFloat,"0.5"),(0.0,1.0))
prob = prob_ode_bigfloatlinear

sol =solve(prob,Vern6(),dt=1//2^(2),dense=true)

interpd_1d_big = sol(0:1//2^(7):1)

sol2 =solve(prob,Vern6(),dt=1//2^(7),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_big)) < 1e-7

#plot(sol2.t,interpd)
#plot!(sol.t,sol[:])
#scatter!(sol.t,sol[:])

prob_ode_bigfloatveclinear = ODEProblem(f,[parse(BigFloat,"0.5")],(0.0,1.0))
prob = prob_ode_bigfloatveclinear

sol =solve(prob,Vern6(),dt=1//2^(2),dense=true)

interpd_big = sol(0:1//2^(4):1)

sol2 =solve(prob,Vern6(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_big)) < 1e-7

### BS5()

prob = prob_ode_linear

sol =solve(prob,BS5(),dt=1//2^(1),dense=true,adaptive=false)

interpd_1d_long = sol(0:1//2^(7):1)

sol2 =solve(prob,BS5(),dt=1//2^(7),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_long)) < 2e-7

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,BS5(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,BS5(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 2e-7

### Vern7()

prob = prob_ode_linear

sol =solve(prob,Vern7(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Vern7(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 3e-9

prob = prob_ode_2Dlinear

sol =solve(prob,Vern7(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Vern7(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 5e-9

### Vern8()

prob = prob_ode_linear

sol =solve(prob,Vern8(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Vern8(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1e-7

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,Vern8(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Vern8(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-7

### Vern9()

prob = prob_ode_linear

sol =solve(prob,Vern9(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Vern9(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 1e-9

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,Vern9(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Vern9(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 2e-9

### Rosenbrock23()

prob = prob_ode_linear

sol =solve(prob,Rosenbrock23(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Rosenbrock23(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 3e-3

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,Rosenbrock23(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Rosenbrock23(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 6e-3

### Rosenbrock32()

prob = prob_ode_linear

sol =solve(prob,Rosenbrock32(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Rosenbrock32(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 3e-4

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,Rosenbrock32(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Rosenbrock32(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 6e-4

### ImplicitEuler()

prob = prob_ode_linear

sol =solve(prob,ImplicitEuler(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,ImplicitEuler(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 0.2

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,ImplicitEuler(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,ImplicitEuler(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < .354

### Trapezoid()

prob = prob_ode_linear

sol =solve(prob,Trapezoid(),dt=1//2^(2),dense=true)

sol(interpd_1d,0:1//2^(4):1)

sol2 =solve(prob,Trapezoid(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d)) < 7e-3

# plot(sol2.t,interpd)

prob = prob_ode_2Dlinear

sol =solve(prob,Trapezoid(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,Trapezoid(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1.4e-2

### DP8()

prob = prob_ode_linear

sol =solve(prob,DP8(),dt=1//2^(2),dense=true)

sol(interpd_1d_long,0:1//2^(7):1)

sol2 =solve(prob,DP8(),dt=1//2^(7),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_long)) < 2e-7

prob = prob_ode_2Dlinear

sol =solve(prob,DP8(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,DP8(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 2.01e-7

### ExplicitRK()

prob = prob_ode_linear

sol =solve(prob,ExplicitRK(),dt=1//2^(2),dense=true)

sol(interpd_1d_long,0:1//2^(7):1)

sol2 =solve(prob,ExplicitRK(),dt=1//2^(7),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_long)) < 6e-5

prob = prob_ode_2Dlinear

sol =solve(prob,ExplicitRK(),dt=1//2^(2),dense=true)

sol(interpd,0:1//2^(4):1)

sol2 =solve(prob,ExplicitRK(),dt=1//2^(4),dense=true,adaptive=false)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 1e-3
