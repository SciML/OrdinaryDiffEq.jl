using OrdinaryDiffEq, Base.Test, DiffEqDevTools

# Ad-hoc functor type for matrices
struct MatrixOperator{T <: AbstractMatrix}
    A::T
end
Base.eltype(op::MatrixOperator) = eltype(op.A)
Base.size(op::MatrixOperator, args...) = size(op.A, args...)
Base.norm(op::MatrixOperator, args...) = norm(op.A, args...)
Base. *(op::MatrixOperator, B) = op.A * B
Base. *(B, op::MatrixOperator) = B * op.A
Base.A_mul_B!(Y, op::MatrixOperator, B) = A_mul_B!(Y, op.A, B)
(op::MatrixOperator)(u,p,t) = op.A * u
(op::MatrixOperator)(du,u,p,t) = A_mul_B!(du, op.A, u)

N = 20
b = 0.2
tspan = (0.0, 1.0)
dts = 1./2.^(7:-1:4) #14->7 good plot
#########################################
# Test problem 1 (heat equation with radiative loss, u' = Au - bu)
dd = -2 * ones(N); du = ones(N - 1)
A = spdiagm((du,dd,du), (-1,0,1))
f1 = MatrixOperator(A)
f2 = (u,p,t) -> -b * u

srand(0); u0 = rand(N)
prob = SplitODEProblem(f1,f2,u0,tspan)
function (::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t)
    tmp = (A - b*I)*t
    return expm(full(tmp)) * u0
end

sim = test_convergence(dts,prob,LawsonEulerKrylov())
@test abs(sim.𝒪est[:l2]-1) < 0.1
sim = test_convergence(dts,prob,ExpEulerKrylov())
@test abs(sim.𝒪est[:l2]-1) < 0.1
#########################################
# Test problem 2 (Schrodinger Equation, iu' = -Au - bu)
# Uses in-place functions
f1 = MatrixOperator(1im * A)
f2 = (du,u,p,t) -> du .= (1im*b) * u
u0 = complex(u0)
prob = SplitODEProblem(f1,f2,u0,tspan)
function (::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t)
    tmp = (A + b*I) * (1im*t)
    return expm(full(tmp)) * u0
end

sim = test_convergence(dts,prob,LawsonEulerKrylov())
@test abs(sim.𝒪est[:l2]-1) < 0.1
sim = test_convergence(dts,prob,ExpEulerKrylov())
@test abs(sim.𝒪est[:l2]-1) < 0.1
