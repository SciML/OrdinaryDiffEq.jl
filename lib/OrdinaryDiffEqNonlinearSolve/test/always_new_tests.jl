using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqBDF
using OrdinaryDiffEqSDIRK
using SciMLBase
using Test

@testset "NLNewton(always_new = true) on Robertson" begin
    function rober!(du, u, p, t)
        y1, y2, y3 = u
        k1, k2, k3 = p
        du[1] = -k1 * y1 + k3 * y2 * y3
        du[2] = k1 * y1 - k2 * y2^2 - k3 * y2 * y3
        du[3] = k2 * y2^2
        return nothing
    end
    rober(u, p, t) = (du = similar(u); rober!(du, u, p, t); du)
    p = [0.04, 3.0e7, 1.0e4]
    u0 = [1.0, 0.0, 0.0]
    tspan = (0.0, 1.0e5)

    @testset "$(nameof(alg)) $(iip ? "iip" : "oop")" for alg in (FBDF, QNDF, TRBDF2, ImplicitEuler),
            iip in (true, false)

        prob = ODEProblem{iip}(iip ? rober! : rober, u0, tspan, p)
        ref = solve(prob, alg(), abstol = 1.0e-8, reltol = 1.0e-8)
        sol = solve(
            prob, alg(nlsolve = NLNewton(always_new = true)),
            abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-4
    end

    # The Jacobian must be taken at the Newton iterate's stage value, so its last
    # evaluation point lies next to the converged step.
    @testset "Jacobian evaluation point, $(nameof(alg))" for alg in (FBDF, TRBDF2)
        jac_points = Vector{Vector{Float64}}()
        function rober_jac!(J, u, p, t)
            push!(jac_points, copy(u))
            y1, y2, y3 = u
            k1, k2, k3 = p
            J[1, 1] = -k1
            J[1, 2] = k3 * y3
            J[1, 3] = k3 * y2
            J[2, 1] = k1
            J[2, 2] = -2k2 * y2 - k3 * y3
            J[2, 3] = -k3 * y2
            J[3, 1] = 0
            J[3, 2] = 2k2 * y2
            J[3, 3] = 0
            return nothing
        end
        prob = ODEProblem(ODEFunction(rober!; jac = rober_jac!), u0, (0.0, 1.0), p)
        integrator = init(
            prob, alg(nlsolve = NLNewton(always_new = true)),
            abstol = 1.0e-8, reltol = 1.0e-8
        )
        for _ in 1:20
            step!(integrator)
        end
        @test SciMLBase.check_error(integrator) == SciMLBase.ReturnCode.Success
        @test jac_points[end] ≈ integrator.u rtol = 1.0e-3
    end
end
