using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqRosenbrock
using Test

const prob_ip = prob_dde_constant_1delay_ip
const prob_scalar = prob_dde_constant_1delay_scalar
const ts = 0:0.1:10

# ODE algorithms
const algs = [
    Rosenbrock23(), Rosenbrock32(), ROS3P(), Rodas3(),
    RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
    Ros4LStab(), Rodas4(), Rodas42(), Rodas4P(), Rodas5(),
]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in algs
    println(nameof(typeof(alg)))

    stepsalg = MethodOfSteps(alg)
    sol_ip = solve(prob_ip, stepsalg)
    sol_scalar = solve(prob_scalar, stepsalg)

    @test isapprox(sol_ip(ts, idxs = 1), sol_scalar(ts), rtol = 1.0e-3)
    @test isapprox(sol_ip.t, sol_scalar.t, rtol = 1.0e-3)
    @test isapprox(sol_ip[1, :], sol_scalar.u, rtol = 1.0e-3)
end
