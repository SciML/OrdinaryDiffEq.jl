using StochasticDiffEq, LinearAlgebra, SparseArrays, Random, LinearSolve, Test
using OrdinaryDiffEqDifferentiation: OrdinaryDiffEqDifferentiation, calc_W!, calc_W
using SciMLOperators: MatrixOperator, WOperator
using OrdinaryDiffEq

#horid nasty hack to deal with temporary calc_W refactor
# if there is a method that takes a W_transform argument, define the version that doesn't to set W_transform to true
if hasmethod(calc_W, (Any, Any, Any, Any, Any))
    function OrdinaryDiffEqDifferentiation.calc_W(integ, nlsolver, dgamma, repeat_step::Bool)
        return OrdinaryDiffEqDifferentiation.calc_W(integ, nlsolver, dgamma, repeat_step, true)
    end
    function OrdinaryDiffEqDifferentiation.calc_W!(
            integ, nlsolver, cache, dgamma, repeat_step::Bool
        )
        return OrdinaryDiffEqDifferentiation.calc_W(
            integ, nlsolver, dgamma, cacherepeat_step, true
        )
    end
end
@testset "Derivative Utilities" begin
    @testset "calc_W!" begin
        A = [-1.0 0.0; 0.0 -0.5]
        σ = [0.9 0.0; 0.0 0.8]
        mm = [2.0 0.0; 0.0 1.0]
        u0 = [1.0, 1.0]
        tmp = zeros(2)
        tspan = (0.0, 1.0)
        dt = 0.01
        dtgamma = 0.5dt
        concrete_W = -mm / dtgamma + A

        # Out-of-place
        _f = (u, p, t) -> A * u
        _g = (u, p, t) -> σ * u
        fun = SDEFunction(
            _f, _g;
            mass_matrix = mm,
            jac = (u, p, t) -> A
        )
        prob = SDEProblem(fun, u0, tspan)
        integrator = init(prob, ImplicitEM(theta = 1); adaptive = false, dt = dt)
        W = calc_W(integrator, integrator.cache.nlsolver, dtgamma, #=repeat_step=# false)
        @test convert(AbstractMatrix, W) ≈ concrete_W
        @test W \ u0 ≈ concrete_W \ u0

        # In-place
        _f = (du, u, p, t) -> mul!(du, A, u)
        _g = (du, u, p, t) -> mul!(du, σ, u)
        fun = SDEFunction(
            _f, _g;
            mass_matrix = mm,
            jac_prototype = MatrixOperator(A)
        )
        prob = SDEProblem(fun, u0, tspan)
        integrator = init(prob, ImplicitEM(theta = 1); adaptive = false, dt = dt)
        W = integrator.cache.nlsolver.cache.W
        calc_W!(
            W, integrator, integrator.cache.nlsolver,
            integrator.cache, dtgamma,
            #=repeat_step=#
            false
        )

        # Did not update because it's an array operator
        # We don't want to build Jacobians when we have operators!
        @test convert(AbstractMatrix, integrator.cache.nlsolver.cache.W) != concrete_W
        ldiv!(tmp, lu!(integrator.cache.nlsolver.cache.W), u0)
        @test tmp != concrete_W \ u0

        # But jacobian2W! will update the cache
        OrdinaryDiffEqDifferentiation.jacobian2W!(
            integrator.cache.nlsolver.cache.W._concrete_form,
            mm, dtgamma, integrator.cache.nlsolver.cache.W.J.A
        )
        @test convert(AbstractMatrix, integrator.cache.nlsolver.cache.W) == concrete_W
        ldiv!(tmp, lu!(integrator.cache.nlsolver.cache.W), u0)
        @test tmp == concrete_W \ u0
    end

    @testset "Implicit solver with lazy W" begin
        A = sparse([-1.0 0.0; 0.0 -0.5])
        σ = sparse([0.9 0.0; 0.0 0.8])
        mm = sparse([2.0 0.0; 0.0 1.0])
        u0 = [1.0, 1.0]
        tspan = (0.0, 1.0)

        _f = (u, p, t) -> t * (A * u)
        _f_ip = (du, u, p, t) -> lmul!(t, mul!(du, A, u))
        _g = (u, p, t) -> σ * u
        _g_ip = (du, u, p, t) -> mul!(du, σ, u)
        prob1 = SDEProblem(SDEFunction(_f, _g; mass_matrix = mm), u0, tspan)
        prob2 = SDEProblem(SDEFunction(_f, _g; mass_matrix = mm, jac = (u, p, t) -> t * A), u0, tspan)
        prob1_ip = SDEProblem(SDEFunction(_f_ip, _g_ip; mass_matrix = mm), u0, tspan)
        jac_prototype = MatrixOperator(
            similar(A); update_func! = (
                J, u, p, t,
            ) -> (J .= t .* A; J)
        )
        prob2_ip = SDEProblem(
            SDEFunction(_f_ip, _g_ip; mass_matrix = mm, jac_prototype = jac_prototype), u0, tspan
        )

        for Alg in [ImplicitEM, ISSEM]
            println(Alg)
            Random.seed!(0)
            sol1 = solve(prob1, Alg(theta = 1); adaptive = false, dt = 0.01)
            Random.seed!(0)
            sol2 = solve(prob2, Alg(theta = 1); adaptive = false, dt = 0.01)
            @test sol1(1.0) ≈ sol2(1.0) rtol = 1.0e-2
            Random.seed!(0)
            sol1_ip = solve(prob1_ip, Alg(theta = 1); adaptive = false, dt = 0.01)
            Random.seed!(0)
            sol2_ip = solve(prob2_ip, Alg(theta = 1); adaptive = false, dt = 0.01)
            @test sol1_ip(1.0) ≈ sol2_ip(1.0) rtol = 1.0e-4
        end

        σ = 1.0
        _g = (u, p, t) -> σ
        _g_ip = (du, u, p, t) -> (du .= σ)
        prob1 = SDEProblem(SDEFunction(_f, _g), u0, tspan)
        prob2 = SDEProblem(SDEFunction(_f, _g; jac = (u, p, t) -> t * A), u0, tspan)
        prob1_ip = SDEProblem(SDEFunction(_f_ip, _g_ip), u0, tspan)
        prob2_ip = SDEProblem(SDEFunction(_f_ip, _g_ip; jac_prototype = jac_prototype), u0, tspan)

        println(SKenCarp)
        Random.seed!(0)
        sol1 = solve(prob1, SKenCarp(); adaptive = false, dt = 0.01)
        Random.seed!(0)
        sol2 = solve(prob2, SKenCarp(); adaptive = false, dt = 0.01)
        @test sol1(1.0) ≈ sol2(1.0) rtol = 1.0e-2
        Random.seed!(0)
        sol1_ip = solve(prob1_ip, SKenCarp(); adaptive = false, dt = 0.01)
        Random.seed!(0)
        sol2_ip = solve(prob2_ip, SKenCarp(); adaptive = false, dt = 0.01)
        @test sol1_ip(1.0) ≈ sol2_ip(1.0) rtol = 1.0e-3
    end
end
