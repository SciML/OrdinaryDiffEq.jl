using OrdinaryDiffEq, Test, LinearAlgebra
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqFIRK, OrdinaryDiffEqRKN,
    OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqLowStorageRK, OrdinaryDiffEqLinear
using RecursiveArrayTools, SciMLOperators

@testset "remake with du0 on SecondOrderODEProblem (#2138)" begin
    f(du, u, p, t) = 0.0
    prob = SecondOrderODEProblem(f, 0.0, 0.0, (0.0, 1.0), [0.0, 0.0])
    new_prob = remake(prob; u0 = 1.0, du0 = 1.0)
    sol = solve(new_prob, DPRKN6())
    @test sol.retcode == ReturnCode.Success
end

@testset "NamedArrayPartition with stiff implicit solvers (#2707)" begin
    function int!(du, u, p, t)
        du.y .= u.y .+ t / 10
        return du.x .= u.x .+ 1.0
    end
    u0 = NamedArrayPartition((x = zeros(10), y = zeros(10, 10)))
    prob = ODEProblem(int!, u0, (0.0, 2.0))
    sol = solve(prob, Trapezoid())
    @test sol.retcode == ReturnCode.Success
end

@testset "RDPK3Sp solvers with Float32 (#2894)" begin
    foo(du, u, p, t) = (du .= 0.9f0 .* u)
    prob = ODEProblem(foo, fill(0.5f0, 10), (0.0f0, 1.0f0))
    @testset "$alg" for alg in (RDPK3Sp35(), RDPK3Sp49(), RDPK3Sp510())
        sol = solve(prob, alg)
        @test sol.retcode == ReturnCode.Success
    end
end

@testset "stats.nf on SecondOrderODEProblem with implicit methods (#2537)" begin
    function test_nf(alg)
        f_counter = Ref(0)
        function f(ddu, du, u, p, t)
            f_counter[] += 1
            ddu .= p .* u
            return nothing
        end
        u0 = [1.0]
        du0 = [0.0]
        p = [-1.0]
        prob = SecondOrderODEProblem(f, du0, u0, (0.0, 1.0), p)
        sol = solve(prob, alg, save_everystep = false, dense = false)
        return sol.stats.nf == f_counter[]
    end
    @testset "$alg" for alg in (RadauIIA5(), Trapezoid(), ImplicitEuler())
        @test test_nf(alg)
    end
end

@testset "explicit jac + jac_prototype skips AD Jacobian prep (#2671)" begin
    f!(du, u, p, t) = (du .= u)
    N = 10
    jac_prototype = Matrix{Float64}(I, N, N)
    function jac!(J, u, p, t)
        J .= 0
        for i in 1:N
            J[i, i] = 1
        end
        return nothing
    end
    odef = ODEFunction(f!; jac = jac!, jac_prototype = jac_prototype)
    prob = ODEProblem(odef, ones(N), (0.0, 1.0))
    sol = solve(prob, Trapezoid())
    @test sol.retcode == ReturnCode.Success
    @test sol(1.0) ≈ exp.(ones(N)) rtol = 1.0e-2
end

@testset "MagnusGL6 with MatrixOperator (#3232)" begin
    function update_func(A, u, p, t)
        A[1, 1] = cos(t)
        A[2, 1] = sin(t)
        A[1, 2] = -sin(t)
        return A[2, 2] = cos(t)
    end
    A = MatrixOperator(ones(2, 2), update_func! = update_func)
    prob = ODEProblem(A, ones(2), (1.0, 6.0))
    sol = solve(prob, MagnusGL6(), dt = 1 / 10)
    @test sol.retcode == ReturnCode.Success
end
