using OrdinaryDiffEqExponentialRK, Test, LinearAlgebra
using SciMLOperators: MatrixOperator, ScalarOperator
using SciMLBase: successful_retcode

@testset "ExpRK6s16 numeric types" begin
    for iip in (false, true), krylov in (false, true)
        A = MatrixOperator(reshape(Float32[-1], 1, 1))
        g = if iip
            (du::Vector{Float32}, u::Vector{Float32}, p, t::Float32) -> (du .= t)
        else
            (u::Vector{Float32}, p, t::Float32) -> Float32[t]
        end
        prob = SplitODEProblem{iip}(A, g, Float32[1], (0.0f0, 1.0f0))
        sol = solve(prob, ExpRK6s16(; krylov, m = 1); dt = 0.25f0)
        @test successful_retcode(sol)
        @test eltype(sol.u[end]) === Float32
        @test sol.u[end][1] ≈ 2exp(-1.0f0) rtol = sqrt(eps(Float32))
    end
    setprecision(BigFloat, 192) do
        prob = SplitODEProblem(
            ScalarOperator(BigFloat(0)),
            (u, p, t) -> 6t^5, BigFloat(1), (BigFloat(0), BigFloat(1))
        )
        sol = solve(prob, ExpRK6s16(); dt = BigFloat(1) / 4)
        @test successful_retcode(sol)
        @test abs(sol.u[end] - 2) < 10^6 * eps(BigFloat)
    end
end

@testset "ExpRK6s16 RHS statistics" begin
    for iip in (false, true), split in (false, true)
        calls = Ref(0)
        f = if iip
            (du, u, p, t) -> (calls[] += 1; du .= -u)
        else
            (u, p, t) -> (calls[] += 1; -u)
        end
        prob = if split
            SplitODEProblem{iip}(MatrixOperator(zeros(1, 1)), f, [1.0], (0.0, 1.0))
        else
            jac = iip ? ((J, u, p, t) -> (J .= -1)) : ((u, p, t) -> -ones(1, 1))
            ODEProblem(ODEFunction{iip}(f; jac), [1.0], (0.0, 1.0))
        end
        integrator = init(prob, ExpRK6s16(; krylov = !split, m = 1); dt = 0.25)
        @test calls[] == (split ? integrator.stats.nf2 : integrator.stats.nf)
        step!(integrator)
        @test calls[] == (split ? integrator.stats.nf2 : integrator.stats.nf)
    end
end

@testset "ExpRK6s16 shared operators" begin
    prob = SplitODEProblem(
        MatrixOperator(-Matrix{Float64}(I, 2, 2)),
        (u, p, t) -> zero(u), ones(2), (0.0, 1.0)
    )
    integrator = init(prob, ExpRK6s16(); dt = 0.25)
    ops = integrator.cache.ops
    @test ops[1] === ops[2] === ops[4] === ops[7] === ops[11]
    @test ops[3] === ops[10] === ops[14]
    @test ops[5] === ops[8] === ops[12]
    @test ops[6] === ops[9] === ops[13]
    @test ops[15] === ops[16]
    sol = solve(prob, ExpRK6s16(); dt = 0.25)
    @test sol.u[end] ≈ fill(exp(-1), 2)
end
