using OrdinaryDiffEq, Test

f(u,p,t) = 0.98u
u0 = 1.0
tspan = (0.0, 1.0)
prob=ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5())

# Make sure various differentiation forms work on scalars

sol1 = solve(prob,Rosenbrock23(),abstol=1e-12,reltol=1e-12)
sol2 = solve(prob,Rosenbrock23(autodiff=false),abstol=1e-12,reltol=1e-12)
sol3 = solve(prob,Rosenbrock23(autodiff=false,diff_type=Val{:central}),abstol=1e-12,reltol=1e-12)
sol4 = solve(prob,Rosenbrock23(autodiff=false,diff_type=Val{:complex}),abstol=1e-12,reltol=1e-12)
sol5 = solve(prob,KenCarp4(),abstol=1e-12,reltol=1e-12)
sol6 = solve(prob,KenCarp4(autodiff=false),abstol=1e-12,reltol=1e-12)
sol7 = solve(prob,KenCarp4(autodiff=false,diff_type=Val{:central}),abstol=1e-12,reltol=1e-12)
sol8 = solve(prob,KenCarp4(autodiff=false,diff_type=Val{:complex}),abstol=1e-12,reltol=1e-12)

ts = 0.0:0.1:1.0
@test sol1(ts) ≈ sol2(ts)
@test sol1(ts) ≈ sol3(ts)
@test sol1(ts) ≈ sol4(ts)
@test sol1(ts) ≈ sol5(ts)
@test sol5(ts) ≈ sol6(ts)
@test sol5(ts) ≈ sol7(ts)
@test sol5(ts) ≈ sol8(ts)

# Test that in-place and out-of-place gets the same results

function f_ip(du, u, p, t)
    du[1] = - u[1]
    nothing
end
f_scalar(u, p, t) = -u

prob_ip = ODEProblem(f_ip, [1.0], (0.0, 10.0))
prob_scalar = ODEProblem(f_scalar, 1.0, (0.0, 10.0))
ts = 0:0.1:10.0

rk_algs = [Euler(),Midpoint(),Heun(),Ralston(),RK4(),SSPRK104(),SSPRK22(),SSPRK33(),
        SSPRK432(),BS3(),BS5(),DP5(),DP8(),Feagin10(),Feagin12(),
        Feagin14(),TanYam7(),Tsit5(),TsitPap8(),Vern6(),Vern7(),Vern8(),Vern9()]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in rk_algs
  println(nameof(typeof(alg)))
  sol_ip = solve(prob_ip, alg)
  sol_scalar = solve(prob_scalar, alg)

  @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
  @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
end

sdirk_algs = [GenericImplicitEuler(), GenericTrapezoid(),
              ImplicitEuler(), ImplicitMidpoint(), Trapezoid(),
              TRBDF2(), SDIRK2(), SSPSDIRK2(),
              Kvaerno3(), KenCarp3(),
              Cash4(), Hairer4(), Hairer42(), Kvaerno4(), KenCarp4(),
              Kvaerno5(), KenCarp5()]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in sdirk_algs
  println(nameof(typeof(alg)))
  sol_ip = solve(prob_ip, alg)
  sol_scalar = solve(prob_scalar, alg)

  @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
  @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
end

rosenbrock_algs = [Rosenbrock23(), Rosenbrock32(), ROS3P(), Rodas3(),
              RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
              Ros4LStab(), Rodas4(), Rodas42(), Rodas4P(), Rodas5()]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in rosenbrock_algs
  println(nameof(typeof(alg)))
  sol_ip = solve(prob_ip, alg)
  sol_scalar = solve(prob_scalar, alg)

  @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
  @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
end

rkc_algs = [RKC(),ROCK2(),ROCK4(),SERK2()]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in rkc_algs
  println(nameof(typeof(alg)))
  sol_ip = solve(prob_ip, alg)
  sol_scalar = solve(prob_scalar, alg)

  @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
  @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
end
