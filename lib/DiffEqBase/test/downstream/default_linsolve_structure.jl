using LinearAlgebra, OrdinaryDiffEq, Test, ADTypes
f = (du, u, p, t) -> du .= u ./ t
jac = (J, u, p, t) -> (J[1, 1] = 1 / t; J[2, 2] = 1 / t; J[1, 2] = 0; J[2, 1] = 0)

jp_diag = Diagonal(zeros(2))
fun = ODEFunction(f; jac = jac, jac_prototype = jp_diag)
prob = ODEProblem(fun, ones(2), (1.0, 10.0))
sol = solve(prob, Rosenbrock23())
@test sol.u[end] ≈ [10.0, 10.0]
@test length(sol) < 60

sol = solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
@test sol.u[end] ≈ [10.0, 10.0]
@test length(sol) < 60

jp = Tridiagonal(jp_diag)
fun = ODEFunction(f; jac = jac, jac_prototype = jp)
prob = ODEProblem(fun, ones(2), (1.0, 10.0))

sol = solve(prob, Rosenbrock23())
@test sol.u[end] ≈ [10.0, 10.0]
@test length(sol) < 60

sol = solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
@test sol.u[end] ≈ [10.0, 10.0]
@test length(sol) < 60

#=
jp = SymTridiagonal(jp_diag)
fun = ODEFunction(f; jac=jac, jac_prototype=jp)
prob = ODEProblem(fun,ones(2),(1.0,10.0))
sol = solve(prob,Rosenbrock23())
@test sol[end] ≈ [10.0,10.0]
@test length(sol) < 60
=#

# Don't test the autodiff=false version here because it's not as numerically stable,
# so lack of optimizations would lead to unsymmetric which causes an error:
# LoadError: ArgumentError: broadcasted assignment breaks symmetry between locations (1, 2) and (2, 1)

@test_broken begin
    local jp = Hermitian(jp_diag)
    local fun = ODEFunction(f; jac = jac, jac_prototype = jp)
    local prob = ODEProblem(fun, ones(2), (1.0, 10.0))
    local sol = solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
    @test sol.u[end] ≈ [10.0, 10.0]
    @test length(sol) < 60
end

@test_broken begin
    local jp = Symmetric(jp_diag)
    local fun = ODEFunction(f; jac = jac, jac_prototype = jp)
    local prob = ODEProblem(fun, ones(2), (1.0, 10.0))
    local sol = solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
    @test sol.u[end] ≈ [10.0, 10.0]
    @test length(sol) < 60
end
