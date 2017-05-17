using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test, DiffEqBase
using Calculus, ForwardDiff


# use `PRINT_TESTS = true` to print the tests, including results
const PRINT_TESTS = false
print_results(x) = if PRINT_TESTS; @printf("%s \n", x) end


# points and storage arrays used in the interpolation tests
const interpolation_points = 0:1//2^(4):1
const interpolation_results_1d = zeros(typeof(prob_ode_linear.u0), length(interpolation_points))
const interpolation_results_2d = Vector{typeof(prob_ode_2Dlinear.u0)}(length(interpolation_points))
for idx in eachindex(interpolation_results_2d)
  interpolation_results_2d[idx] = zeros(prob_ode_2Dlinear.u0)
end

# perform the regression tests
# NOTE: If you want to add new tests (for new algorithms), you have to run the
#       commands below to get numerical values for `tol_ode_linear` and
#       `tol_ode_2Dlinear`.
function regression_test(alg, tol_ode_linear, tol_ode_2Dlinear)
  PRINT_TESTS && println("\n", alg)

  sol  = solve(prob_ode_linear, alg, dt=1//2^(2), dense=true)
  sol(interpolation_results_1d, interpolation_points)
  sol(interpolation_points[1])
  sol2 = solve(prob_ode_linear, alg, dt=1//2^(4), dense=true, adaptive=false)
  print_results( @test maximum(abs.(sol2[:] - interpolation_results_1d)) < tol_ode_linear )

  sol  = solve(prob_ode_2Dlinear, alg, dt=1//2^(2), dense=true)
  sol(interpolation_results_2d,  interpolation_points)
  sol(interpolation_points[1])
  sol2 = solve(prob_ode_2Dlinear, alg, dt=1//2^(4), dense=true, adaptive=false)
  print_results( @test maximum(maximum.(abs.(sol2[:] - interpolation_results_2d))) < tol_ode_2Dlinear)
end


# Some extra tests using Euler()
prob = prob_ode_linear
sol  = solve(prob, Euler(), dt=1//2^(2), dense=true)
interpd_1d = sol(0:1//2^(4):1)
sol2 = solve(prob, Euler(), dt=1//2^(4), dense=true)
sol3 = solve(prob, Euler(), dt=1//2^(5), dense=true)

prob = prob_ode_2Dlinear
sol  = solve(prob, Euler(), dt=1//2^(2), dense=true)
interpd = sol(0:1//2^(4):1)

interpd_idxs = sol(0:1//2^(4):1,idxs=1:2:5)

@test minimum([interpd_idxs[i] == interpd[i][1:2:5] for i in 1:length(interpd)])

interpd_single = sol(0:1//2^(4):1,idxs=1)

@test typeof(interpd_single) <: Vector{Float64}

@test typeof(sol(0.5,idxs=1)) <: Float64

A = rand(4,2)
sol(A,0.777)
A == sol(0.777)
sol2 = solve(prob, Euler(), dt=1//2^(4), dense=true)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < .2

sol(interpd, 0:1//2^(4):1)

@test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < .2

sol  = solve(prob, Euler(), dt=1//2^(2), dense=false)

@test !(sol(0.6)[4,2] ≈ 0)


# Euler
regression_test(Euler(), 0.2, 0.2)

# Midpoint
regression_test(Midpoint(), 1.5e-2, 2.3e-2)

# SSPRK22
regression_test(SSPRK22(), 1.5e-2, 2.5e-2)

# SSPRK33
regression_test(SSPRK33(), 7.5e-4, 1.5e-3)

# SSPRK104
regression_test(SSPRK104(), 1.5e-5, 3e-5)

# RK4
regression_test(RK4(), 4.5e-5, 1e-4)

# DP5
regression_test(DP5(), 5e-6, 1e-5)

# BS3
regression_test(BS3(), 5e-4, 8e-4)

# Tsit5
regression_test(Tsit5(), 2e-6, 4e-6)

# TanYam7
regression_test(TanYam7(), 4e-4, 6e-4)

# TsitPap8
regression_test(TsitPap8(), 1e-3, 3e-3)

# Feagin10
regression_test(Feagin10(), 6e-4, 9e-4)

# Vern6
regression_test(Vern6(), 7e-8, 7e-8)

const linear_bigα = parse(BigFloat, "1.01")
f = (t,u) -> (linear_bigα*u)
prob_ode_bigfloatlinear = ODEProblem(f,parse(BigFloat,"0.5"),(0.0,1.0))
prob = prob_ode_bigfloatlinear

sol  = solve(prob, Vern6(), dt=1//2^(2), dense=true)
interpd_1d_big = sol(0:1//2^(7):1)
sol2 = solve(prob, Vern6(), dt=1//2^(7), dense=true, adaptive=false)
print_results( @test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_big)) < 5e-8 )

prob_ode_bigfloatveclinear = ODEProblem(f,[parse(BigFloat,"0.5")],(0.0,1.0))
prob = prob_ode_bigfloatveclinear
sol  = solve(prob, Vern6(), dt=1//2^(2), dense=true)
interpd_big = sol(0:1//2^(4):1)
sol2 = solve(prob, Vern6(), dt=1//2^(4), dense=true, adaptive=false)
print_results( @test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_big)) < 5e-8 )

# BS5
regression_test(BS5(), 4e-8, 6e-8)

prob = prob_ode_linear
sol  = solve(prob, BS5(), dt=1//2^(1), dense=true, adaptive=false)
interpd_1d_long = sol(0:1//2^(7):1)
sol2 = solve(prob, BS5(), dt=1//2^(7), dense=true, adaptive=false)
print_results( @test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_long)) < 2e-7 )

# Vern7
regression_test(Vern7(), 3e-9, 5e-9)

# Vern8
regression_test(Vern8(), 3e-8, 5e-8)

# Vern9
regression_test(Vern9(), 1e-9, 2e-9)

# Rosenbrock23
regression_test(Rosenbrock23(), 3e-3, 6e-3)

prob = prob_ode_linear
sol  = solve(prob, Rosenbrock23(), dt=1//2^(12), dense=true)
sol(0:1//2^(4):1)
sol(0:1//2^(4):1,Val{1})

const deriv_test_points = linspace(0,1,10)
for t in deriv_test_points
  deriv = sol(t,Val{1})
  if t == 0
    #@test deriv ≈ derivative(sol,0.00,:forward)
  elseif t != 1
    #@test deriv ≈ derivative(sol,t)
  end
  @test deriv ≈ ForwardDiff.derivative(sol,t)
end

# Rosenbrock32
regression_test(Rosenbrock32(), 4e-4, 6e-4)

# ImplicitEuler
regression_test(ImplicitEuler(), 0.2, 0.354)

# Trapezoid
regression_test(Trapezoid(), 7e-3, 1.4e-2)

# DP8
regression_test(DP8(), 2e-7, 3e-7)

prob = prob_ode_linear
sol  = solve(prob, DP8(), dt=1//2^(2), dense=true)
sol(interpd_1d_long,0:1//2^(7):1)
sol2 = solve(prob, DP8(), dt=1//2^(7), dense=true, adaptive=false)
print_results( @test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_long)) < 2e-7 )

# ExplicitRK
regression_test(ExplicitRK(), 7e-5, 2e-4)

prob = prob_ode_linear
sol  = solve(prob, ExplicitRK(), dt=1//2^(2), dense=true)
sol(interpd_1d_long,0:1//2^(7):1)
sol2 = solve(prob, ExplicitRK(), dt=1//2^(7), dense=true, adaptive=false)
print_results( @test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd_1d_long)) < 6e-5 )

prob = prob_ode_2Dlinear
sol  = solve(prob, ExplicitRK(), dt=1//2^(2), dense=true)
sol(interpd, 0:1//2^(4):1)
sol2 = solve(prob, ExplicitRK(), dt=1//2^(4), dense=true, adaptive=false)
print_results( @test maximum(map((x)->maximum(abs.(x)),sol2[:] - interpd)) < 2e-4 )
