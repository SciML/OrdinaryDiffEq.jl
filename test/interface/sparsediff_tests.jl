using Test
using OrdinaryDiffEq
using SparseArrays
using LinearAlgebra

## in-place
#https://github.com/JuliaDiffEq/SparseDiffTools.jl/blob/master/test/test_integration.jl
function f_ip(dx,x,p,t)
  for i in 2:length(x)-1
    dx[i] = x[i-1] - 2x[i] + x[i+1]
  end
  dx[1] = -2x[1] + x[2]
  dx[end] = x[end-1] - 2x[end]
  nothing
end

## out-of-place
function f_oop(x,p,t)
  dx = similar(x)
  for i in 2:length(x)-1
    dx[i] = x[i-1] - 2x[i] + x[i+1]
  end
  dx[1] = -2x[1] + x[2]
  dx[end] = x[end-1] - 2x[end]
  return dx
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

for f in [f_oop, f_ip]
  odefun_std = ODEFunction(f)
  prob_std = ODEProblem(odefun_std,u0,tspan)

  for ad in [true, false]
    for Solver in [Rodas5, Trapezoid, KenCarp4]
      for tol in [nothing, 1e-10]
        sol_std=solve(prob_std,Solver(autodiff=ad),reltol=tol,abstol=tol)
        @test sol_std.retcode==:Success
        for (i,prob) in enumerate(map(f->ODEProblem(f,u0,tspan),
                        [ODEFunction(f,colorvec=colors,jac_prototype=jac_sp),
                         ODEFunction(f,jac_prototype=jac_sp),
                         ODEFunction(f,colorvec=colors)
                         ]))
          # TODO: these broken test-cases need to be investigated.
          #       Note they only happen for prob=ODEFunction(f,colorvec=colors).
          isbroken = i==3 && (
            (f, ad, tol) == (f_oop, true, nothing) ||
            (f, ad, Solver, tol) == (f_oop, false, Trapezoid, nothing) ||
            (f, Solver, tol) == (f_ip, Trapezoid, nothing) ||
            (f, Solver) == (f_oop, Rodas5) ||
            (f, Solver) == (f_ip, Rodas5) ||
            (f, Solver) == (f_oop, KenCarp4) ||
            (f, Solver) == (f_ip, KenCarp4)
          )
          sol=solve(prob,Solver(autodiff=ad),reltol=tol,abstol=tol)
          @test sol.retcode==:Success
          if tol !=nothing
            if isbroken
              @test_broken sol_std.u[end]≈sol.u[end] atol=tol
            else
              @test sol_std.u[end]≈sol.u[end] atol=tol
            end
          else
            if isbroken
              @test_broken sol_std.u[end]≈sol.u[end]
            else
              @test sol_std.u[end]≈sol.u[end]
            end
          end
          if isbroken
            @test_broken length(sol_std.t)==length(sol.t)
          else
            @test length(sol_std.t)==length(sol.t)
          end
        end
      end
    end
  end
end
