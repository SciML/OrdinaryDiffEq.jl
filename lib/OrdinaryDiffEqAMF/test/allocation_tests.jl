using Pkg
Pkg.add("AllocCheck")
Pkg.add("OrdinaryDiffEqRosenbrock")

using OrdinaryDiffEqAMF
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using SciMLOperators: MatrixOperator
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqAMF (Approximate Matrix Factorization) using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.

AMF requires build_amf_function with a MatrixOperator for the Jacobian and
a Rosenbrock algorithm (e.g. ROS34PW1a) as the underlying integrator.

AMF is marked broken=true as it currently allocates in perform_step!.
"""

@testset "AMF Allocation Tests" begin
    function f!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    J_op = MatrixOperator([-0.5 0.0; 0.0 -1.5])
    amf_func = build_amf_function(f!; jac = J_op)
    prob = ODEProblem(amf_func, [1.0, 1.0], (0.0, 1.0))

    @testset "AMF(ROS34PW1a) perform_step! Static Analysis" begin
        integrator = init(
            prob, AMF(ROS34PW1a), dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )
        step!(integrator)

        cache = integrator.cache
        allocs = check_allocs(
            OrdinaryDiffEqCore.perform_step!,
            (typeof(integrator), typeof(cache))
        )

        @test length(allocs) == 0 broken = true

        if length(allocs) > 0
            println(
                "AllocCheck found $(length(allocs)) allocation sites in AMF(ROS34PW1a) perform_step!"
            )
        else
            println("AMF(ROS34PW1a) perform_step! appears allocation-free with AllocCheck")
        end
    end
end
