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

    # Rosenbrock32 IIP vs scalar DDEs can diverge more due to
    # different initdt paths, so use a wider tolerance
    _rtol = alg isa Rosenbrock32 ? 5.0e-2 : 1.0e-2
    @test isapprox(sol_ip(ts, idxs = 1), sol_scalar(ts), rtol = _rtol)
    # Compare endpoints: in-place and scalar may take different step counts
    @test isapprox(sol_ip.t[end], sol_scalar.t[end], rtol = _rtol)
    @test isapprox(sol_ip[1, end], sol_scalar.u[end], rtol = _rtol)
end
