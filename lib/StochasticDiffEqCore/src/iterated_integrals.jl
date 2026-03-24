"""
All stochastic iterated integrals are written in Stratonovich sense as indicated
by the J.
"""

abstract type AbstractJ end

# Iterated stochastic integral (for diagonal and commutative noise processes where the Levy area approximation is not required.)
abstract type AbstractJDiagonal <: AbstractJ end
abstract type AbstractJCommute <: AbstractJ end

struct JDiagonal_oop <: AbstractJDiagonal end

mutable struct JDiagonal_iip{JType} <: AbstractJDiagonal
    J::JType
    JDiagonal_iip(ΔW) = new{typeof(ΔW)}(false .* ΔW .* ΔW)
end

struct JCommute_oop <: AbstractJCommute end

mutable struct JCommute_iip{JType} <: AbstractJCommute
    J::JType
    function JCommute_iip(ΔW)
        J = false .* vec(ΔW) .* vec(ΔW)'
        return new{typeof(J)}(J)
    end
end

function get_iterated_I(dt, dW, dZ, alg::JDiagonal_oop, p = nothing, c = 1, γ = 1 // 1)
    return 1 // 2 .* dW .* dW
end

function get_iterated_I!(dt, dW, dZ, alg::JDiagonal_iip, p = nothing, c = 1, γ = 1 // 1)
    (; J) = alg
    if dW isa Number
        alg.J = 1 // 2 .* dW .^ 2
    else
        @.. J = 1 // 2 * dW^2
    end
    return nothing
end

function get_iterated_I(dt, dW, dZ, alg::JCommute_oop, p = nothing, c = 1, γ = 1 // 1)
    J = 1 // 2 .* vec(dW) .* vec(dW)'
    return J
end

function get_iterated_I!(dt, dW, dZ, alg::JCommute_iip, p = nothing, c = 1, γ = 1 // 1)
    (; J) = alg
    mul!(J, vec(dW), vec(dW)')
    @.. J *= 1 // 2
    return nothing
end

# algs from LevyArea.jl # LevyArea.levyarea allocates random variables and then mutates these, see e.g.
# https://github.com/stochastics-uni-luebeck/LevyArea.jl/blob/68c5cb08ab103b4dcd3178651f7a5dd9ce8c666d/src/milstein.jl#L25
function get_iterated_I(
        dt, dW, dZ, alg::LevyArea.AbstractIteratedIntegralAlgorithm,
        p = nothing, c = 1, γ = 1 // 1
    )
    if isnothing(p)
        ε = c * dt^(γ + 1 // 2)
        p = terms_needed(length(dW), dt, ε, alg, MaxL2())
    end
    I = LevyArea.levyarea(dW / √dt, p, alg)
    return I .= 1 // 2 * dW .* dW' .+ dt .* I
end

mutable struct IteratedIntegralAlgorithm_iip{JType, levyalgType} <:
    LevyArea.AbstractIteratedIntegralAlgorithm
    J::JType
    levyalg::levyalgType
    function IteratedIntegralAlgorithm_iip(ΔW, levyalg)
        J = false .* vec(ΔW) .* vec(ΔW)'
        return new{typeof(J), typeof(levyalg)}(J, levyalg)
    end
end

function get_iterated_I!(
        dt, dW, dZ, alg::IteratedIntegralAlgorithm_iip, p = nothing, c = 1, γ = 1 // 1
    )
    (; J, levyalg) = alg
    if isnothing(p)
        ε = c * dt^(γ + 1 // 2)
        p = terms_needed(length(dW), dt, ε, levyalg, MaxL2())
    end
    J .= LevyArea.levyarea(dW / √dt, p, levyalg)
    J .= 1 // 2 * dW .* dW' .+ dt .* J
    return nothing
end

# Default algorithms, keep KPWJ_oop() to have a non-mutating version
function get_Jalg(ΔW, dt, prob, alg)
    return if alg.ii_approx isa IILevyArea
        if isinplace(prob)
            if ΔW isa Number || is_diagonal_noise(prob)
                return JDiagonal_iip(ΔW)
            else
                # optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm=MaxL2())
                return IteratedIntegralAlgorithm_iip(ΔW, LevyArea.optimal_algorithm(length(ΔW), dt))
            end
        else
            if ΔW isa Number || is_diagonal_noise(prob)
                return JDiagonal_oop()
            else
                return LevyArea.optimal_algorithm(length(ΔW), dt)
            end
        end
    elseif alg.ii_approx isa IICommutative
        if isinplace(prob)
            if ΔW isa Number || is_diagonal_noise(prob)
                return JDiagonal_iip(ΔW)
            else
                return JCommute_iip(ΔW)
            end
        else
            if ΔW isa Number || is_diagonal_noise(prob)
                return JDiagonal_oop()
            else
                return JCommute_oop()
            end
        end
    else
        if isinplace(prob)
            IteratedIntegralAlgorithm_iip(ΔW, alg.ii_approx)
        else
            return alg.ii_approx
        end
    end
end

# Specific Levy area alg for an SDE solver
# function StochasticDiffEq.get_Jalg(ΔW,prob,alg::SOLVER)
#  return MronRoe()
# end

"""
    compute_iterated_I_from_noise(W, t, dt)

Compute Stratonovich iterated integrals J_{jk} = ∫∫ ∘dW_j ∘dW_k over [t, t+dt]
from the sub-grid W values stored in the noise process. Returns the m×m matrix J
where J_{jk} = (1/2)*dW_j*dW_k + dt*A_{jk} (A = Lévy area).

This is computed via the Riemann sum over sub-grid increments:
  I_{jk}^{Ito} = Σ_{n} dW_k^{(n)} * Σ_{l<n} dW_j^{(l)}
then converted to Stratonovich: J_{jk} = I_{jk} + (1/2)*δ_{jk}*dt.

Falls back to `nothing` if the noise process doesn't have accessible sub-grid data.
"""
function compute_iterated_I_from_noise(W, t, dt)
    # Only works for NoiseGrid and NoiseWrapper with accessible grid
    source = _get_noise_source(W)
    source === nothing && return nothing

    t_grid = source.t
    W_grid = source.W

    m = length(W_grid[1])
    t_end = t + dt

    # Find grid indices covering [t, t+dt]
    i_start = searchsortedfirst(t_grid, t)
    i_end = searchsortedlast(t_grid, t_end)

    # Need at least 2 sub-steps to get useful iterated integrals
    if i_end - i_start < 2
        return nothing
    end

    # Compute Ito iterated integrals from sub-grid increments
    T = eltype(eltype(W_grid))
    I = zeros(T, m, m)
    W_cumsum = zeros(T, m)  # running sum of dW from t_start

    for n in (i_start + 1):i_end
        dWn = W_grid[n] .- W_grid[n - 1]
        for k in 1:m
            for j in 1:m
                I[j, k] += W_cumsum[j] * dWn[k]
            end
        end
        W_cumsum .+= dWn
    end

    # Convert Ito to Stratonovich: J_{jk} = I_{jk} + (1/2)*δ_{jk}*dt
    for j in 1:m
        I[j, j] += dt / 2
    end

    return I
end

function _get_noise_source(W)
    if hasproperty(W, :source)
        # NoiseWrapper — get the underlying source
        return _get_noise_source(W.source)
    elseif W isa DiffEqNoiseProcess.NoiseGrid || W isa DiffEqNoiseProcess.NoiseProcess
        return W
    else
        return nothing
    end
end
