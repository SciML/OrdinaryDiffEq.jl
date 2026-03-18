using OrdinaryDiffEq, DiffEqBase, CUDA, LinearAlgebra, Test, StaticArrays, ADTypes
function f(u, p, t)
    return A * u
end
function f(du, u, p, t)
    return mul!(du, A, u)
end
function jac(J, u, p, t)
    return J .= A
end
function jac(u, p, t)
    return A
end
function tgrad(du, u, p, t)
    return du .= 0
end
function tgrad(u, p, t)
    return zero(u)
end
ff = ODEFunction(f, jac = jac, tgrad = tgrad)
CUDA.allowscalar(false)
A = cu(-rand(3, 3))
u0 = cu([1.0; 0.0; 0.0])
tspan = (0.0f0, 100.0f0)

prob = ODEProblem(ff, u0, tspan)
sol = solve(prob, Tsit5())
@test solve(prob, Rosenbrock23()).retcode == ReturnCode.Success
solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()));

prob_oop = ODEProblem{false}(ff, u0, tspan)
CUDA.allowscalar(false)
sol = solve(prob_oop, Tsit5())
@test solve(prob_oop, Rosenbrock23()).retcode == ReturnCode.Success
@test solve(prob_oop, Rosenbrock23(autodiff = AutoFiniteDiff())).retcode ==
    ReturnCode.Success

prob_nojac = ODEProblem(f, u0, tspan)
@test solve(prob_nojac, Rosenbrock23()).retcode == ReturnCode.Success
@test solve(prob_nojac, Rosenbrock23(autodiff = AutoFiniteDiff())).retcode ==
    ReturnCode.Success
@test solve(
    prob_nojac,
    Rosenbrock23(autodiff = AutoFiniteDiff(; fdtype = Val(:central)))
).retcode ==
    ReturnCode.Success
@test solve(
    prob_nojac,
    Rosenbrock23(autodiff = AutoFiniteDiff(; fdtype = Val(:complex)))
).retcode ==
    ReturnCode.Success

#=
prob_nojac_oop = ODEProblem{false}(f,u0,tspan)
DiffEqBase.prob2dtmin(prob_nojac_oop)
@test_broken solve(prob_nojac_oop,Rosenbrock23()).retcode == ReturnCode.Success
@test_broken solve(prob_nojac_oop,Rosenbrock23(autodiff=AutoFiniteDiff())).retcode == ReturnCode.Success
@test_broken solve(prob_nojac_oop,Rosenbrock23(autodiff=AutoFiniteDiff(; fdtype = Val(:central))).retcode == ReturnCode.Success
@test_broken solve(prob_nojac_oop,Rosenbrock23(autodiff=AutoFiniteDiff(; fdtype = Val(:complex))).retcode == ReturnCode.Success
=#

# Complex Numbers Adaptivity DifferentialEquations.jl#460
f_complex(u, nothing, t) = 5.0f-1 .* u
u0 = cu(rand(32, 32) .+ 1im * rand(32, 32));
prob = ODEProblem(f_complex, u0, (0.0f0, 1.0f0))
@test_nowarn sol = solve(prob, Tsit5())

# Calculating norm of Static Arrays in GPU kernel DiffEqBase.jl#864

function test_SA_norm(u::T) where {T <: AbstractArray}
    @cushow DiffEqBase.ODE_DEFAULT_NORM(u, 1.0)
    return nothing
end

u = @SVector rand(100)

@testset "Static arrays norm on GPU" begin
    @cuda test_SA_norm(u)
end
