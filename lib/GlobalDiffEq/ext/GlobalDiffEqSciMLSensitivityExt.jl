module GlobalDiffEqSciMLSensitivityExt

import GlobalDiffEq, LinearAlgebra, QuadGK, SciMLBase, SciMLSensitivity

function GlobalDiffEq._default_quadrature_sensealg()
    return SciMLSensitivity.QuadratureAdjoint(autojacvec = true)
end

GlobalDiffEq._is_quadrature_adjoint(::SciMLSensitivity.QuadratureAdjoint) = true

function GlobalDiffEq._adjoint_solution(
        sol, sensealg::SciMLSensitivity.QuadratureAdjoint, adjoint_alg, direction;
        abstol, reltol
    )
    terminal_gradient! = let direction = direction
        function (out, u, p, t, i)
            copyto!(out, direction)
            return nothing
        end
    end
    terminal_time = sol.prob.tspan[2]
    adjoint_prob = SciMLSensitivity.ODEAdjointProblem(
        sol, sensealg, adjoint_alg, [terminal_time], terminal_gradient!
    )
    adjoint_sol = SciMLBase.solve(
        adjoint_prob, adjoint_alg;
        abstol, reltol, dense = true, save_everystep = true
    )
    SciMLBase.successful_retcode(adjoint_sol) ||
        throw(ErrorException("the adjoint solve failed with retcode $(adjoint_sol.retcode)"))
    return adjoint_sol
end

function GlobalDiffEq._defect_projection(sol, adjoint_sol; abstol, reltol)
    prob = sol.prob
    isinplace = SciMLBase.isinplace(prob)

    function integrand(t)
        u = sol(t, continuity = :right)
        du = sol(t, Val{1}, continuity = :right)
        rhs = if isinplace
            value = similar(u)
            prob.f(value, u, prob.p, t)
            value
        else
            prob.f(u, prob.p, t)
        end
        lambda = adjoint_sol(t, continuity = :right)
        return LinearAlgebra.dot(lambda, du - rhs)
    end

    projection = zero(eltype(prob.u0))
    for i in 1:(length(sol.t) - 1)
        left = sol.t[i]
        right = sol.t[i + 1]
        left == right && continue
        interval_projection, _ = QuadGK.quadgk(
            integrand, left, right; atol = abstol, rtol = reltol
        )
        projection += interval_projection
    end
    return projection
end

end
