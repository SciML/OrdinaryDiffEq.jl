"""
    GlobalRichardson(alg)

Wrap an ODE algorithm with global Richardson extrapolation.
"""
struct GlobalRichardson{A <: SciMLBase.AbstractODEAlgorithm} <: GlobalDiffEqAlgorithm
    alg::A
end

# Forward algorithm traits to the wrapped algorithm
# This allows GlobalRichardson to inherit capabilities from the inner algorithm
SciMLBase.allows_arbitrary_number_types(alg::GlobalRichardson) =
    SciMLBase.allows_arbitrary_number_types(alg.alg)
SciMLBase.allowscomplex(alg::GlobalRichardson) =
    SciMLBase.allowscomplex(alg.alg)
SciMLBase.isautodifferentiable(alg::GlobalRichardson) =
    SciMLBase.isautodifferentiable(alg.alg)

function SciMLBase.__solve(
        prob::Union{SciMLBase.AbstractODEProblem, SciMLBase.AbstractDAEProblem},
        alg::GlobalRichardson, args...;
        dt, kwargs...
    )
    opt = Dict(kwargs)
    otheropts = delete!(copy(opt), :dt)
    tstops = get(opt, :tstops, range(prob.tspan[1], stop = prob.tspan[2], step = dt))
    local sol
    val,
        err = Richardson.extrapolate(
        dt, rtol = get(opt, :reltol, 1.0e-3),
        atol = get(opt, :abstol, 1.0e-6), contract = 0.5
    ) do _dt
        sol = solve(prob, alg.alg, args...; dt = _dt, adaptive = false, otheropts...)
        # Convert Vector{Vector{T}} to Matrix{T} for Richardson.jl compatibility
        reduce(hcat, sol.(tstops))
    end
    return sol
end
