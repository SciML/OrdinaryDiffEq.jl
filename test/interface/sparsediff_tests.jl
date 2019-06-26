using Test
using OrdinaryDiffEq

#https://github.com/JuliaDiffEq/SparseDiffTools.jl/blob/master/test/test_integration.jl
function f(dx,x,p,t)
    for i in 2:length(x)-1
      dx[i] = x[i-1] - 2x[i] + x[i+1]
    end
    dx[1] = -2x[1] + x[2]
    dx[end] = x[end-1] - 2x[end]
    nothing
  end

function second_derivative_stencil(N)
  A = zeros(N,N)
  for i in 1:N, j in 1:N
    (j-i==-1 || j-i==1) && (A[i,j]=1)
    j-i==0 && (A[i,j]=-2)
  end
  A
end
  

colors = repeat(1:3,10)
odefun_sp = ODEFunction(f,colorvec=colors)
prob_sp = ODEProblem(odefun_sp,[1.,2.,3,4,5,4,3,2,1],(0.,10.))
prob_std = ODEProblem(f,[1.,2.,3,4,5,4,3,2,1],(0.,10.))
sol_sp=solve(prob_sp,Rodas5(autodiff=false),abstol=1e-10,reltol=1e-10)
@test sol_sp.retcode==:Success#test sparse finitediff
sol=solve(prob_std,Rodas5(autodiff=false),abstol=1e-10,reltol=1e-10)
@test sol_sp.u[end]≈sol.u[end] atol=1e-10
@test length(sol_sp.t)==length(sol.t)
sol_sp=solve(prob_sp,Rodas5(autodiff=false))
sol=solve(prob_std,Rodas5(autodiff=false))
@test sol_sp.u[end]≈sol.u[end] atol=1e-3
@test length(sol_sp.t)==length(sol.t)
sol_sp=solve(prob_sp,Rodas5())
sol=solve(prob_std,Rodas5())
@test sol_sp.u[end]≈sol.u[end] atol=1e-3
@test length(sol_sp.t)==length(sol.t)