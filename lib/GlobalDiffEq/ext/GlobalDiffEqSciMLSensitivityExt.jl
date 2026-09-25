module GlobalDiffEqSciMLSensitivityExt

import GlobalDiffEq, SciMLBase, SciMLSensitivity
import Accessors: @set

# Endpoint global-error projection along `direction`, computed entirely by the
# adjoint system. The dual-weighted residual ∫ λ(t)ᵀ (f(P(t)) − P'(t)) dt is
# obtained as the sensitivity of the discrete cost g = ⟨direction, u(T)⟩ to a
# scalar forcing ε that adds the numerical solution's defect r(t) = f(P) − P' to
# the dynamics: with `u' = f(u,p,t) + ε·r(t)`, dG/dε|₀ = ∫ λᵀ r dt exactly. The
# adjoint (λ' = −(∂f/∂u)ᵀλ, λ(T) = direction) is unchanged because r depends only
# on t, so `adjoint_sensitivities` returns that integral as the parameter
# gradient — no separate adjoint solve or hand-rolled quadrature.
function GlobalDiffEq._adjoint_defect_projection(
        sol, sensealg, adjoint_alg, direction;
        abstol, reltol, terminal_time = sol.prob.tspan[2]
    )
    prob = sol.prob
    f_orig = SciMLBase.unwrapped_f(prob.f)
    realp = prob.p
    inplace = SciMLBase.isinplace(prob)

    forced_f = if inplace
        let f = f_orig, realp = realp, sol = sol
            function (du, u, forcing, t)
                f(du, u, realp, t)
                residual = similar(du)
                f(residual, sol(t), realp, t)
                deriv = sol(t, Val{1})
                @. du += forcing[1] * (residual - deriv)
                return nothing
            end
        end
    else
        let f = f_orig, realp = realp, sol = sol
            function (u, forcing, t)
                residual = f(sol(t), realp, t) .- sol(t, Val{1})
                return f(u, realp, t) .+ forcing[1] .* residual
            end
        end
    end

    forced_prob = SciMLBase.remake(
        prob;
        f = SciMLBase.ODEFunction{inplace}(forced_f),
        p = [zero(eltype(prob.u0))]
    )
    forced_sol = @set sol.prob = forced_prob

    terminal_gradient! = let direction = direction
        function (out, u, p, t, i)
            copyto!(out, direction)
            return nothing
        end
    end
    # A `nothing` sensealg defers to `adjoint_sensitivities`' own default; any
    # user-passed adjoint method is forwarded unchanged.
    sensealg_kwargs = sensealg === nothing ? (;) : (; sensealg)
    _, defect_gradient = SciMLSensitivity.adjoint_sensitivities(
        forced_sol, adjoint_alg;
        sensealg_kwargs...,
        t = [terminal_time], dgdu_discrete = terminal_gradient!,
        abstol, reltol
    )
    return only(defect_gradient)
end

end
