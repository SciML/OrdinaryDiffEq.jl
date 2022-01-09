using OrdinaryDiffEq, Test, Random, LinearAlgebra
Random.seed!(123)

A = 0.01*rand(3, 3)
rn = (du, u, p, t) -> begin
    mul!(du, A, u)
end
u0 = rand(3)
prob = ODEProblem(rn, u0, (0, 10.))

function precsl(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = lu(convert(AbstractMatrix,W))
  else
    Pl = Plprev
  end
  Pl,nothing
end

function precsr(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pr = lu(convert(AbstractMatrix,W))
  else
    Pr = Prprev
  end
  nothing,Pr
end

function precslr(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pr = lu(convert(AbstractMatrix,W))
  else
    Pr = Prprev
  end
  Pr,Pr
end

sol = @test_nowarn solve(prob, TRBDF2(autodiff=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=IterativeSolversJL_GMRES()));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=IterativeSolversJL_GMRES(), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=IterativeSolversJL_GMRES(precs=precsl), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=IterativeSolversJL_GMRES(precs=precsr), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=IterativeSolversJL_GMRES(precs=precslr), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, QNDF(autodiff=false, linsolve=IterativeSolversJL_GMRES()));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, Rosenbrock23(autodiff=false, linsolve=IterativeSolversJL_GMRES(precs=precslr), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, Rodas4(autodiff=false, linsolve=IterativeSolversJL_GMRES(precs=precslr), smooth_est=false));
@test length(sol.t) < 20


sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=KrylovJL_GMRES()));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=KrylovJL_GMRES(), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=KrylovJL_GMRES(precs=precsl), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=KrylovJL_GMRES(precs=precsr), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=KrylovJL_GMRES(precs=precslr), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, QNDF(autodiff=false, linsolve=KrylovJL_GMRES()));
@test length(sol.t) < 20
ol = @test_nowarn solve(prob, Rosenbrock23(autodiff=false, linsolve=KrylovJL_GMRES(precs=precslr), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, Rodas4(autodiff=false, linsolve=KrylovJL_GMRES(precs=precslr), smooth_est=false));
@test length(sol.t) < 20
