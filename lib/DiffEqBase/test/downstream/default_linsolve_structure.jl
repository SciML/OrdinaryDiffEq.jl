using LinearAlgebra, OrdinaryDiffEq, Test, ADTypes
using LinearSolve: LUFactorization
using OrdinaryDiffEqRosenbrock: Rosenbrock23
f = (du, u, p, t) -> du .= u ./ t
function jac(J, u, p, t)
    fill!(J, zero(eltype(J)))
    for i in axes(J, 1)
        J[i, i] = 1 / t
    end
    return nothing
end

jp_diag = Diagonal(zeros(3))
structured_prototypes = (
    (jp_diag, Diagonal(ones(3))),
    (Bidiagonal(zeros(3), zeros(2), :U), Bidiagonal(ones(3), ones(2), :U)),
    (Tridiagonal(jp_diag), Tridiagonal(ones(2), ones(3), ones(2))),
    (SymTridiagonal(zeros(3), zeros(2)), SymTridiagonal(ones(3), ones(2))),
    (UpperTriangular(zeros(3, 3)), UpperTriangular(ones(3, 3))),
)

for (jp, expected_cache) in structured_prototypes
    fun = ODEFunction(f; jac = jac, jac_prototype = jp)
    prob = ODEProblem(fun, ones(3), (1.0, 10.0))

    alg = jp isa UpperTriangular ? Rosenbrock23(linsolve = LUFactorization()) : Rosenbrock23()
    integ = init(prob, alg)
    @test integ.cache.J isa typeof(jp)
    @test integ.cache.W isa typeof(jp)
    @test integ.cache.J == expected_cache
    @test integ.cache.W == expected_cache

    if jp isa Union{Diagonal, Tridiagonal}
        sol = solve!(integ)
        @test sol.u[end] ≈ [10.0, 10.0, 10.0]
        @test length(sol.t) < 60

        integ = init(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
        sol = solve!(integ)
        @test sol.u[end] ≈ [10.0, 10.0, 10.0]
        @test length(sol.t) < 60
    end
end

#=
jp = SymTridiagonal(Diagonal(zeros(2)))
fun = ODEFunction(f; jac=jac, jac_prototype=jp)
prob = ODEProblem(fun,ones(2),(1.0,10.0))
sol = solve(prob,Rosenbrock23())
@test sol[end] ≈ [10.0,10.0]
@test length(sol.t) < 60
=#

# Don't test the autodiff=false version here because it's not as numerically stable,
# so lack of optimizations would lead to unsymmetric which causes an error:
# LoadError: ArgumentError: broadcasted assignment breaks symmetry between locations (1, 2) and (2, 1)

@test_broken begin
    local jp = Hermitian(Diagonal(zeros(2)))
    local fun = ODEFunction(f; jac = jac, jac_prototype = jp)
    local prob = ODEProblem(fun, ones(2), (1.0, 10.0))
    local sol = solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
    @test sol.u[end] ≈ [10.0, 10.0]
    @test length(sol.t) < 60
end

@test_broken begin
    local jp = Symmetric(Diagonal(zeros(2)))
    local fun = ODEFunction(f; jac = jac, jac_prototype = jp)
    local prob = ODEProblem(fun, ones(2), (1.0, 10.0))
    local sol = solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
    @test sol.u[end] ≈ [10.0, 10.0]
    @test length(sol.t) < 60
end
