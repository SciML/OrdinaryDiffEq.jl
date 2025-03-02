using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, van

testTol = 0.5

for prob in [prob_ode_linear, prob_ode_2Dlinear]
    sim21 = test_convergence(1 .// 2 .^ (6:-1:3), prob, RadauIIA5(), dense_errors = true)
    @test sim21.𝒪est[:final]≈5 atol=testTol
    @test sim21.𝒪est[:L2]≈4 atol=testTol
end

sim21 = test_convergence(1 .// 2 .^ (6:-1:3), prob_ode_linear, RadauIIA5(), dense_errors = true)
@test sim21.𝒪est[:final]≈5 atol=testTol
@test sim21.𝒪est[:L2]≈4 atol=testTol


sim21 = test_convergence(1 ./ 2 .^ (2.5:-1:0.5), prob_ode_linear, RadauIIA9(), dense_errors = true)
@test sim21.𝒪est[:final]≈8 atol=testTol
@test sim21.𝒪est[:L2]≈6 atol=testTol

sim21 = test_convergence(1 ./ 2 .^ (2.5:-1:0.5), prob_ode_2Dlinear, RadauIIA9(), dense_errors = true)
@test sim21.𝒪est[:final]≈8 atol=testTol
@test sim21.𝒪est[:L2]≈6 atol=testTol

using GenericSchur

prob_ode_linear_big = remake(
    prob_ode_linear, u0 = big.(prob_ode_linear.u0), tspan = big.(prob_ode_linear.tspan))
prob_ode_2Dlinear_big = remake(prob_ode_2Dlinear, u0 = big.(prob_ode_2Dlinear.u0),
    tspan = big.(prob_ode_2Dlinear.tspan))

#non-threaded tests
for i in [5, 9, 13, 17, 21, 25], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    local sim21 = test_convergence(dts, prob, AdaptiveRadau(min_order = i, max_order = i), dense_errors = true)
    @test sim21.𝒪est[:final] ≈ i atol=testTol
    @test sim21.𝒪est[:L2] ≈ ((i + 3) ÷ 2) atol=testTol
end

dts = 1 ./ 2 .^ (4.25:-1:0.25)
local sim21 = test_convergence(dts, prob_ode_2Dlinear_big, AdaptiveRadau(min_order = 5, max_order = 5), dense_errors = true)

#threaded tests
using OrdinaryDiffEqCore
for i in [5, 9, 13, 17, 21, 25], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    local sim21 = test_convergence(dts, prob, AdaptiveRadau(min_order = i, max_order = i, threading = OrdinaryDiffEqCore.PolyesterThreads()))
    @test sim21.𝒪est[:final] ≈ i atol=testTol
end

# test adaptivity
for iip in (true, false)
    if iip
        vanstiff = ODEProblem{iip}(van, [0; sqrt(3)], (0.0, 1.0), 1e6)
    else
        vanstiff = ODEProblem{false}((u, p, t) -> van(u, p, t), [0; sqrt(3)], (0.0, 1.0),
            1e6)
    end
    sol = solve(vanstiff, RadauIIA5())
    if iip
        @test sol.stats.naccept + sol.stats.nreject > sol.stats.njacs # J reuse
        @test sol.stats.njacs < sol.stats.nw # W reuse
    end
    @test length(sol) < 150
    @test length(solve(remake(vanstiff, p = 1e7), RadauIIA5())) < 150
    @test length(solve(remake(vanstiff, p = 1e7), reltol = [1e-4, 1e-6], RadauIIA5())) < 170
    @test length(solve(remake(vanstiff, p = 1e7), RadauIIA5(), reltol = 1e-9,
        abstol = 1e-9)) < 870
    @test length(solve(remake(vanstiff, p = 1e9), RadauIIA5())) < 170
    @test length(solve(remake(vanstiff, p = 1e10), RadauIIA5())) < 190
end

##Tests for RadauIIA3
for prob in [prob_ode_linear, prob_ode_2Dlinear]
    dts = 1 ./ 2 .^ (8:-1:1)
    sim = test_convergence(dts, prob, RadauIIA3(), dense_errors = true)
    @test sim.𝒪est[:final]≈3 atol=0.25
    @test sim.𝒪est[:L2]≈2 atol=0.25
end

# test adaptivity
for iip in (true, false)
    if iip
        vanstiff = ODEProblem{iip}(van, [0; sqrt(3)], (0.0, 1.0), 1e6)
    else
        vanstiff = ODEProblem{false}((u, p, t) -> van(u, p, t), [0; sqrt(3)], (0.0, 1.0),
            1e6)
    end
    sol = solve(vanstiff, RadauIIA3())
    if iip
        @test sol.stats.naccept + sol.stats.nreject > sol.stats.njacs # J reuse
        @test sol.stats.njacs < sol.stats.nw # W reuse
    end
    @test length(sol) < 5000 # the error estimate is not very good
end

using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, 
      ODEInterfaceDiffEq, LSODA, SparseArrays, LinearSolve,
      LinearAlgebra, IncompleteLU, AlgebraicMultigrid, Symbolics, ModelingToolkit, OrdinaryDiffEqFIRK, Plots
using OrdinaryDiffEqCore
using DifferentialEquations
gr()
const N = 8
xyd_brusselator = range(0,stop=1,length=N)
brusselator_f(x, y, t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.
limit(a, N) = a == N+1 ? 1 : a == 0 ? N : a
function brusselator_2d_loop(du, u, p, t)
  A, B, alpha, dx = p
  alpha = alpha/dx^2
  @inbounds for I in CartesianIndices((N, N))
    i, j = Tuple(I)
    x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
    ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
    du[i,j,1] = alpha*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) +
                B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1] + brusselator_f(x, y, t)
    du[i,j,2] = alpha*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
end
p = (3.4, 1., 10., step(xyd_brusselator))
input = rand(N,N,2)
output = similar(input)
sparsity_pattern = Symbolics.jacobian_sparsity(brusselator_2d_loop,output,input,p,0.0)
jac_sparsity = Float64.(sparse(sparsity_pattern))
f = ODEFunction{true, SciMLBase.FullSpecialize}(brusselator_2d_loop;jac_prototype=jac_sparsity)
function init_brusselator_2d(xyd)
  N = length(xyd)
  u = zeros(N, N, 2)
  for I in CartesianIndices((N, N))
    x = xyd[I[1]]
    y = xyd[I[2]]
    u[I,1] = 22*(y*(1-y))^(3/2)
    u[I,2] = 27*(x*(1-x))^(3/2)
  end
  u
end
u0 = init_brusselator_2d(xyd_brusselator)
prob = ODEProblem(f,u0,(0.,5.),p);
test_sol = solve(prob,RadauIIA5(),abstol=1/10^25,reltol=1/10^20)
abstols = 1.0 ./ 10.0 .^ (7:12)
reltols = 1.0 ./ 10.0 .^ (4:9)
label = ["FBDF", "CVODE_BDF", "ddebdf", "Rodas5P", "radau", "AdaptiveRadau (No threading)", "AdaptiveRadau (Threaded)"]
setups = [
          Dict(:alg=>FBDF()),
          Dict(:alg=>CVODE_BDF()),
          Dict(:alg=>ddebdf()),
          Dict(:alg=>Rodas5P()),
          Dict(:alg=>radau()),
          Dict(:alg=>AdaptiveRadau(max_order = 25)),
          Dict(:alg=>AdaptiveRadau(max_order=25, threading = OrdinaryDiffEqCore.PolyesterThreads()))
]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;verbose=false, dense=false,
                    save_everystep=false,appxsol=test_sol, names = label, maxiters=Int(1e5))
plot(wp)

setups = [
          Dict(:alg=>radau()),
          Dict(:alg=>AdaptiveRadau(max_order = 25)),
]

setups = [
    Dict(:alg=>radau()),
    Dict(:alg=>AdaptiveRadau(max_order = 25)),
    Dict(:alg=>AdaptiveRadau(max_order = 25, threading = OrdinaryDiffEqCore.PolyesterThreads()))
]

setups = [
          Dict(:alg=>Rodas5P()),
          Dict(:alg=>AdaptiveRadau(max_order = 25)),
]

u = [0.,10.]
g(du, u,p,t) = du.=[-9.8, u[1]]
prob = ODEProblem(g, u, (0,5))
condition(u, t, integrator) = u[2]
affect!(integrator) = integrator.u[1] *= -1
cb = ContinuousCallback(condition, affect!)
sol = solve(prob, RadauIIA3(), callback = cb)
plot(sol)