using OrdinaryDiffEq, DiffEqProblemLibrary, Test

@testset "saveat Tests" begin
  prob = prob_ode_linear

  sol =solve(prob,DP5(),dt=1//2^(2),save_everystep=false,dense=false)
  sol2=solve(prob,DP5(),dt=1//2^(2),save_everystep=false,dense=false,saveat=[1/2])

  @test symdiff(sol.t,sol2.t) == [1/2]

  sol2=solve(prob,DP5(),dt=1//2^(2),save_everystep=false,dense=false,saveat=1/2)

  @test symdiff(sol.t,sol2.t) == [1/2]

  sol3=solve(prob,DP5(),dt=1//2^(2),save_everystep=false,dense=false,saveat=[1/2],tstops=[1/2])

  @test sol3.t == [0.0,0.5,1.00]

  sol3=solve(prob,DP5(),dt=1//2^(2),saveat=[1/2],tstops=[1/2])

  @test sol3.t == [0.0,0.5,1.00]

  sol3=solve(prob,DP5(),dt=1//2^(2),saveat=1/10,tstops=[1/2])

  @test sol3.t == collect(0.0:0.1:1.00)

  #plot(0:1//2^(4):1,interpd)

  sol =solve(prob,RK4(),dt=1/2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,RK4(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[.125,.6,.61,.8])

  @test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

  sol =solve(prob,Rosenbrock32(),dt=1/2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,Rosenbrock32(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[.125,.6,.61,.8])

  @test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

  sol =solve(prob,GenericTrapezoid(),dt=1/2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,GenericTrapezoid(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[.125,.6,.61,.8])

  @test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

  prob = prob_ode_2Dlinear

  sol =solve(prob,DP5(),dt=1//2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,DP5(),dt=1//2^(2),save_everystep=true,dense=false,saveat=[1/2])

  @test symdiff(sol.t,sol2.t) == [1/2]

  sol =solve(prob,RK4(),dt=1/2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,RK4(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[.125,.6,.61,.8])

  @test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

  sol =solve(prob,Rosenbrock32(),dt=1/2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,Rosenbrock32(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[.125,.6,.61,.8])

  @test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

  sol =solve(prob,GenericTrapezoid(),dt=1/2^(2),save_everystep=true,dense=false)
  sol2=solve(prob,GenericTrapezoid(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[.125,.6,.61,.8])

  @test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

  sol=solve(prob,GenericTrapezoid(),dt=1/2^(2),save_everystep=true,dense=false,saveat=[0,.125,.6,.61,.8])

  @test !(sol.t[2] ≈ 0)

  # Test Iterators

  sol2=solve(prob,DP5(),dt=1//2^(2),save_everystep=false,dense=false,saveat=0:1//100:1)

  @test sol2.t ≈ collect(0:1//100:1)

  sol2=solve(prob,DP5(),dt=1//2^(2),save_everystep=false,dense=false,saveat=linspace(0,1,100))

  @test sol2.t ≈ linspace(0,1,100)

  f = (du,u,p,t) -> prob.f(du,u,p,t)
  prob2 = ODEProblem(f,vec(prob.u0),prob.tspan)

  sol2=solve(prob2,DP5(),dt=1//2^(2),saveat=.1,save_idxs=1:2:5)

  for u in sol2.u
    @test length(u) == 3
  end

  sol2=solve(prob2,DP5(),dt=1//2^(2),saveat=.1,save_idxs=1:2:5,save_everystep=true)

  sol=solve(prob2,DP5(),dt=1//2^(2),save_start=false)

  @test sol.t[1] == 1//2^(2)
end
