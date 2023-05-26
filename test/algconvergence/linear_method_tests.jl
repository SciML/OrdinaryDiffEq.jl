using OrdinaryDiffEq, Test, DiffEqDevTools
using LinearAlgebra, Random

# Linear exponential solvers
A = MatrixOperator([2.0 -1.0; -1.0 2.0])
u0 = ones(2)
prob = ODEProblem(A, u0, (0.0, 1.0))
solve(prob, LinearExponential(krylov = :off))

sol1 = solve(prob, LinearExponential(krylov = :off))(1.0)
sol2 = solve(prob, LinearExponential(krylov = :simple))(1.0)
sol3 = solve(prob, LinearExponential(krylov = :adaptive))(1.0)
sol4 = solve(prob, Rosenbrock23(), reltol = 1e-12, abstol = 1e-12)(1.0)
sol_analytic = exp(1.0 * Matrix(A)) * u0

@test isapprox(sol1, sol_analytic, rtol = 1e-10)
@test isapprox(sol2, sol_analytic, rtol = 1e-10)
@test isapprox(sol3, sol_analytic, rtol = 1e-10)
@test isapprox(sol4, sol_analytic, rtol = 1e-8)

# u' = A(t)u solvers
function update_func!(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = sin(u[1])
    A[1, 2] = -1
    A[2, 2] = 0
end
A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (10, 50.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.RKMK2(), dt = 1 / 4)
dts = 1 ./ 2 .^ (10:-1:5)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, RKMK2(), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.RKMK4(), dt = 1 / 4)
dts = (0.38) .^ (6:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, RKMK4(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.22

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.LieRK4(), dt = 1 / 4)
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, LieRK4(), test_setup)
@test sim.𝒪est[:l2]≈5 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern9(), dt = 1 / 4)
sol2 = solve(prob, OrdinaryDiffEq.CG2(), dt = 1 / 4)
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG2(), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 20.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern6(), dt = 1 / 8)
sol2 = solve(prob, OrdinaryDiffEq.CG3(), dt = 1 / 8)
dts = 1 ./ 2 .^ (10:-1:3)
test_setup = Dict(:alg => Vern6(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG3(), test_setup)
@test sim.𝒪est[:l2]≈3 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 30.0))
sol1 = solve(prob, Vern9(), dt = 1 / 4)
sol2 = solve(prob, CG4a(), dt = 1 / 4)
dts = (0.38) .^ (6:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG4a(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2

function update_func!(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = 1
    A[1, 2] = -2 * (1 - cos(u[2]) - u[2] * sin(u[2]))
    A[2, 2] = 0
end
A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (30, 150.0))
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Tsit5(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusAdapt4(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2

function update_func!(A, u, p, t)
    A[1, 1] = cos(t)
    A[2, 1] = sin(t)
    A[1, 2] = -sin(t)
    A[2, 2] = cos(t)
end
A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0.0, 5.0))
dts = 1 ./ 2 .^ (10:-1:1)
sol = solve(prob, OrdinaryDiffEq.MagnusMidpoint(), dt = 1 / 4)

dts = 1 ./ 2 .^ (10:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusMidpoint(), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusMidpoint(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0.0, 5.0))
dts = 1 ./ 2 .^ (10:-1:1)
sol = solve(prob, OrdinaryDiffEq.MagnusLeapfrog(), dt = 1 / 4)

dts = 1 ./ 2 .^ (10:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusLeapfrog(), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusLeapfrog(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈2 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0.5, 5.0))
dts = 1 ./ 2 .^ (10:-1:1)
sol = solve(prob, OrdinaryDiffEq.LieEuler(), dt = 1 / 4)

dts = 1 ./ 2 .^ (10:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, LieEuler(), test_setup)
@test sim.𝒪est[:l2]≈1 atol=0.2
sim = analyticless_test_convergence(dts, prob, LieEuler(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈1 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (1.0, 6.0))
dts = 1 ./ 2 .^ (10:-1:1)
sol = solve(prob, MagnusGauss4(), dt = 1 / 4)

dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Vern6(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusGauss4(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusGauss4(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (1.0, 6.0))
dts = 1 ./ 2 .^ (4:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusNC6(), test_setup)
@test sim.𝒪est[:l2]≈6 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusNC6(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈6 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (1.0, 6.0))
sol = solve(prob, MagnusGL6(), dt = 1 / 10)
dts = 1 ./ 2 .^ (4:-1:1)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusGL6(), test_setup)
@test sim.𝒪est[:l2]≈6 atol=0.3
sim = analyticless_test_convergence(dts, prob, MagnusGL6(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈6 atol=0.3

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0.0, 100.0))
dts = 1.775 .^ (5:-1:0)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusGL8(), test_setup)
@test sim.𝒪est[:l2]≈8 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusGL8(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈8 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0.0, 100.0))
dts = 1.773 .^ (5:-1:0)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusNC8(), test_setup)
@test sim.𝒪est[:l2]≈8 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusNC8(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈8 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (1.0, 6.0))
dts = 1 ./ 2 .^ (7:-1:1)
test_setup = Dict(:alg => Vern6(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, MagnusGL4(), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2
sim = analyticless_test_convergence(dts, prob, MagnusGL4(krylov = true), test_setup)
@test sim.𝒪est[:l2]≈4 atol=0.2

A = MatrixOperator(ones(2, 2), update_func! = update_func!)
prob = ODEProblem(A, ones(2), (0, 20.0))
sol1 = solve(prob, OrdinaryDiffEq.Vern6(), dt = 1 / 8)
sol2 = solve(prob, OrdinaryDiffEq.CG3(), dt = 1 / 8)
dts = 1 ./ 2 .^ (10:-1:3)
test_setup = Dict(:alg => Vern6(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, CG3(), test_setup)
@test sim.𝒪est[:l2]≈3 atol=0.2

function B(y::AbstractMatrix)
    b = similar(y)
    N = size(b, 1)
    for l in 1:N
        for k in 1:l
            if k < l
                @inbounds b[k, l] = y[k, l - 1] - y[k + 1, l]
            else
                @inbounds b[k, l] = 0
            end
        end
    end
    for l in 1:N
        for k in (l + 1):N
            @inbounds b[k, l] = -b[l, k]
        end
    end

    return b
end

function update_func(A, u, p, t)
    B(u)
end

function update_func!(A, u, p, t)
    A .= B(u)
    return nothing
end

η = diagm([1.0, 2, 3, 4, 5])
A = MatrixOperator(Matrix{eltype(η)}(I(size(η, 1))), update_func = update_func, update_func! = update_func!)
dts = 1 ./ 2 .^ (10:-1:2)
tspan = (0.0, 20.0)

# IIP
f = SplitFunction(A, (du, u, p, t) -> du .= -u * B(u), _func_cache = similar(η))
prob = SplitODEProblem(f, η, tspan)
sol = solve(prob, CayleyEuler(), dt = 1 / 10)

@test sol.retcode == ReturnCode.Success
eig_err = [norm(eigvals(sol[i]) - eigvals(η)) for i in eachindex(sol)]
@test all(≈(e, 0, atol = 1e-13) for e in eig_err)

test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)

sim = analyticless_test_convergence(dts, prob, CayleyEuler(), test_setup)
@test sim.𝒪est[:l2]≈1 atol=0.2

# OOP
f = SplitFunction(A, (u, p, t) -> -u * B(u), _func_cache = similar(η))
prob = SplitODEProblem(f, η, tspan)
sol = solve(prob, CayleyEuler(), dt = 1 / 10)

@test sol.retcode == ReturnCode.Success
eig_err = [norm(eigvals(sol[i]) - eigvals(η)) for i in eachindex(sol)]
@test all(≈(e, 0, atol = 1e-13) for e in eig_err)

test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)

sim = analyticless_test_convergence(dts, prob, CayleyEuler(), test_setup)
@test sim.𝒪est[:l2]≈1 atol=0.2
