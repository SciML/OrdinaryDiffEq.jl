using Test
using OrdinaryDiffEq
using SparseArrays
using LinearAlgebra

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

function generate_sparsity_pattern(N::Integer)
  dl = repeat([1.0],N-1)
  du = repeat([1.0],N-1)
  d = repeat([-2.0],N)
  return Tridiagonal(dl,d,du)
end

jac_sp = sparse(generate_sparsity_pattern(10))
#jac = second_derivative_stencil(10)
colors = repeat(1:3,10)[1:10]
u0=[1.,2.,3,4,5,5,4,3,2,1]
tspan=(0.,10.)
odefun_sp = ODEFunction(f,colorvec=colors)
odefun_sp2= ODEFunction(f,colorvec=colors,jac_prototype=jac_sp)
prob_sp = ODEProblem(odefun_sp,u0,tspan)
prob_sp2 = ODEProblem(odefun_sp2,u0,tspan)
prob_std = ODEProblem(f,u0,tspan)

sol_sp=solve(prob_sp,Rodas5(autodiff=false),abstol=1e-10,reltol=1e-10)
sol_sp2=solve(prob_sp2,Rodas5(autodiff=false),abstol=1e-10,reltol=1e-10)
@test sol_sp.retcode==:Success#test sparse finitediff
@test sol_sp2.retcode==:Success
sol=solve(prob_std,Rodas5(autodiff=false),abstol=1e-10,reltol=1e-10)
@test sol_sp.u[end]≈sol.u[end] atol=1e-10
@test sol_sp2.u[end]≈sol.u[end] atol=1e-10
@test length(sol_sp.t)==length(sol.t)
@test length(sol_sp2.t)==length(sol.t)

sol_sp=solve(prob_sp,Rodas5(autodiff=false))
sol_sp2=solve(prob_sp2,Rodas5(autodiff=false))
sol=solve(prob_std,Rodas5(autodiff=false))
@test sol_sp.u[end]≈sol_sp2.u[end] atol=1e-10
@test length(sol_sp.t)==length(sol_sp2.t)

@test sol_sp.u[end]≈sol.u[end] atol=1e-3
@test sol_sp2.u[end]≈sol.u[end] atol=1e-3
@test length(sol_sp.t)==length(sol.t)

sol_sp=solve(prob_sp,Rodas5())
sol_sp2=solve(prob_sp2,Rodas5())
sol=solve(prob_std,Rodas5())
@test sol_sp.u[end]≈sol_sp2.u[end] atol=1e-10
@test length(sol_sp.t)==length(sol_sp2.t)

@test sol_sp.u[end]≈sol.u[end] atol=1e-3
@test length(sol_sp.t)==length(sol.t)