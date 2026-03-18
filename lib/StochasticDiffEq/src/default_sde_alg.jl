# Default algorithm selection for SDEs
# Moved from DifferentialEquations.jl as part of modularization effort

using LinearAlgebra: I

# Helper function to extract alg_hints from keyword arguments
function get_alg_hints(o)
    return :alg_hints ∈ keys(o) ? alg_hints = o[:alg_hints] : alg_hints = Symbol[:auto]
end

function default_algorithm(
        prob::DiffEqBase.AbstractSDEProblem{uType, tType, isinplace, ND};
        kwargs...
    ) where {uType, tType, isinplace, ND}
    o = Dict{Symbol, Any}(kwargs)
    alg = SOSRI() # Standard default

    alg_hints = get_alg_hints(o)

    if :commutative ∈ alg_hints
        alg = RKMilCommute()
    end

    is_stiff = :stiff ∈ alg_hints
    is_stratonovich = :stratonovich ∈ alg_hints
    if is_stiff || prob.f.mass_matrix !== I
        alg = ImplicitRKMil(autodiff = false)
    end

    if is_stratonovich
        if is_stiff || prob.f.mass_matrix !== I
            alg = ImplicitRKMil(
                autodiff = false,
                interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich
            )
        else
            alg = RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich)
        end
    end

    if prob.noise_rate_prototype != nothing || prob.noise != nothing
        if is_stratonovich
            if is_stiff || prob.f.mass_matrix !== I
                alg = ImplicitEulerHeun(autodiff = false)
            else
                alg = LambaEulerHeun()
            end
        else
            if is_stiff || prob.f.mass_matrix !== I
                alg = ISSEM(autodiff = false)
            else
                alg = LambaEM()
            end
        end
    end

    if :additive ∈ alg_hints
        if is_stiff || prob.f.mass_matrix !== I
            alg = SKenCarp(autodiff = false)
        else
            alg = SOSRA()
        end
    end

    return alg
end

SciMLBase.supports_solve_rng(::SciMLBase.AbstractSDEProblem, ::Nothing) = true

# Dispatch for __init with Nothing algorithm - use default
function DiffEqBase.__init(
        prob::DiffEqBase.AbstractSDEProblem, ::Nothing, args...; kwargs...
    )
    alg = default_algorithm(prob; kwargs...)
    return DiffEqBase.__init(prob, alg, args...; kwargs...)
end

# Dispatch for __solve with Nothing algorithm - use default
function DiffEqBase.__solve(
        prob::DiffEqBase.AbstractSDEProblem, ::Nothing, args...; kwargs...
    )
    alg = default_algorithm(prob; kwargs...)
    return DiffEqBase.__solve(prob, alg, args...; kwargs...)
end
