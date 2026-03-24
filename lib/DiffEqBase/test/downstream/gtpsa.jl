using OrdinaryDiffEq, GTPSA, Test
using DifferentiationInterface
using ADTypes: AutoForwardDiff
using ForwardDiff

# GTPSA is itself an AD engine - these tests compare GTPSA jacobian/hessian
# results against ForwardDiff as a reference implementation

# ODEProblem 1 =======================

f!(du, u, p, t) = du .= p .* u

# Initial variables and parameters
x = [1.0, 2.0, 3.0]
p = [4.0, 5.0, 6.0]

prob = ODEProblem(f!, x, (0.0, 1.0), p)
sol = solve(prob, Tsit5(), reltol = 1.0e-16, abstol = 1.0e-16)

# Parametric GTPSA map
desc = Descriptor(3, 2, 3, 2) # 3 variables 3 parameters, both to 2nd order
dx = @vars(desc)
dp = @params(desc)
prob_GTPSA = ODEProblem(f!, x .+ dx, (0.0, 1.0), p .+ dp)
sol_GTPSA = solve(prob_GTPSA, Tsit5(), reltol = 1.0e-16, abstol = 1.0e-16)

@test sol.u[end] ≈ scalar.(sol_GTPSA.u[end]) # scalar gets 0th order part

# Compare Jacobian against AD backends using DifferentiationInterface
function sol_end_problem1(t)
    prob = ODEProblem(f!, t[1:3], (0.0, 1.0), t[4:6])
    sol = solve(prob, Tsit5(), reltol = 1.0e-16, abstol = 1.0e-16)
    return sol.u[end]
end

@testset "GTPSA Problem 1 Jacobian tests" begin
    J_AD = DifferentiationInterface.jacobian(sol_end_problem1, AutoForwardDiff(), [x..., p...])
    @test J_AD ≈ GTPSA.jacobian(sol_GTPSA.u[end], include_params = true)
end

@testset "GTPSA Problem 1 Hessian tests" begin
    for i in 1:3
        function sol_end_i_problem1(t)
            prob = ODEProblem(f!, t[1:3], (0.0, 1.0), t[4:6])
            sol = solve(prob, Tsit5(), reltol = 1.0e-16, abstol = 1.0e-16)
            return sol.u[end][i]
        end
        Hi_AD = DifferentiationInterface.hessian(sol_end_i_problem1, AutoForwardDiff(), [x..., p...])
        @test Hi_AD ≈ GTPSA.hessian(sol_GTPSA.u[end][i], include_params = true)
    end
end

# ODEProblem 2 =======================
pdot!(dq, p, q, params, t) = dq .= [0.0, 0.0, 0.0]
function qdot!(dp, p, q, params, t)
    return dp .= [
        p[1] / sqrt((1 + p[3])^2 - p[1]^2 - p[2]^2),
        p[2] / sqrt((1 + p[3])^2 - p[1]^2 - p[2]^2),
        p[3] / sqrt(1 + p[3]^2) - (p[3] + 1) / sqrt((1 + p[3])^2 - p[1]^2 - p[2]^2),
    ]
end

prob2 = DynamicalODEProblem(pdot!, qdot!, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], (0.0, 25.0))
sol2 = solve(prob2, Yoshida6(), dt = 1.0, reltol = 1.0e-16, abstol = 1.0e-16)

desc2 = Descriptor(6, 2) # 6 variables to 2nd order
dx2 = @vars(desc2) # identity map
prob_GTPSA2 = DynamicalODEProblem(pdot!, qdot!, dx2[1:3], dx2[4:6], (0.0, 25.0))
sol_GTPSA2 = solve(prob_GTPSA2, Yoshida6(), dt = 1.0, reltol = 1.0e-16, abstol = 1.0e-16)

@test sol2.u[end] ≈ scalar.(sol_GTPSA2.u[end]) # scalar gets 0th order part

# Compare Jacobian against AD backends using DifferentiationInterface
function sol_end_problem2(t)
    prob = DynamicalODEProblem(pdot!, qdot!, t[1:3], t[4:6], (0.0, 25.0))
    sol = solve(prob, Yoshida6(), dt = 1.0, reltol = 1.0e-16, abstol = 1.0e-16)
    return sol.u[end]
end

@testset "GTPSA Problem 2 Jacobian tests" begin
    J_AD = DifferentiationInterface.jacobian(sol_end_problem2, AutoForwardDiff(), zeros(6))
    @test J_AD ≈ GTPSA.jacobian(sol_GTPSA2.u[end], include_params = true)
end

@testset "GTPSA Problem 2 Hessian tests" begin
    for i in 1:6
        function sol_end_i_problem2(t)
            prob = DynamicalODEProblem(pdot!, qdot!, t[1:3], t[4:6], (0.0, 25.0))
            sol = solve(prob, Yoshida6(), dt = 1.0, reltol = 1.0e-16, abstol = 1.0e-16)
            return sol.u[end][i]
        end
        Hi_AD = DifferentiationInterface.hessian(sol_end_i_problem2, AutoForwardDiff(), zeros(6))
        @test Hi_AD ≈ GTPSA.hessian(sol_GTPSA2.u[end][i], include_params = true)
    end
end
