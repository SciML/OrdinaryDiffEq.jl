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

Compute Stratonovich iterated integrals J_{jk} over [t, t+dt] from sub-step
dW increments available in the noise process. Returns m×m matrix J or `nothing`.

Two sources of sub-step data are checked:

1. **RSwM3 S₂ stack** (adaptive stepping): During `setup_next_step!`, RSwM3
   decomposes the current step's dW into sub-intervals stored in S₂. These
   are the finest available decomposition of the current noise increment.

2. **Saved trajectory grid** (NoiseWrapper/NoiseGrid): When re-solving on a
   previously generated noise path, the source process has fine-grid W values
   stored in W.t/W.W from the original solve.

The Stratonovich integral is computed as:
  J_{jk} = I_{jk}^{Ito} + (1/2)*δ_{jk}*dt
where I_{jk}^{Ito} = Σ_{n} dW_k^{(n)} * Σ_{l<n} dW_j^{(l)} (Riemann sum).
"""
function compute_iterated_I_from_noise(W, t, dt)
    # Strategy 1: Use RSwM3 S₂ stack sub-intervals (adaptive stepping)
    J = _compute_II_from_S2(W, dt)
    J !== nothing && return J

    # Strategy 2: Use saved trajectory grid (NoiseWrapper/NoiseGrid)
    return _compute_II_from_grid(W, t, dt)
end

"""
Compute iterated integrals from RSwM3's S₂ stack, which holds the sub-interval
decomposition (dt_i, dW_i, dZ_i) of the current step's noise increment.
"""
function _compute_II_from_S2(W, dt)
    # S₂ is only populated for RSwM3
    !hasproperty(W, :S₂) && return nothing
    S₂ = W.S₂
    n_sub = length(S₂)
    # Need ≥2 sub-intervals for a meaningful Lévy area estimate.
    # With only 1 entry (fresh step, no prior rejections), the iterated
    # integral would just be the commutative approximation (zero Lévy area),
    # which is worse than LevyArea.jl's random approximation.
    n_sub < 2 && return nothing

    # Extract sub-interval dW values from S₂
    # S₂.data[1:S₂.cur] contains (dt_i, dW_i, dZ_i) in time order
    first_dW = S₂.data[1][2]
    m = length(first_dW)
    T = eltype(first_dW)

    return _iterated_I_from_increments(S₂.data, n_sub, m, T, dt)
end

"""
Compute iterated integrals from saved trajectory W values on a fine grid.
Used for NoiseWrapper/NoiseGrid convergence testing.
"""
function _compute_II_from_grid(W, t, dt)
    source = _get_noise_source(W)
    source === nothing && return nothing

    t_grid = source.t
    W_grid = source.W

    m = length(W_grid[1])
    t_end = t + dt

    # Find grid indices covering [t, t+dt]
    i_start = searchsortedfirst(t_grid, t)
    i_end = searchsortedlast(t_grid, t_end)

    # Need at least 2 sub-steps
    n_sub = i_end - i_start
    n_sub < 2 && return nothing

    T = eltype(eltype(W_grid))
    I = zeros(T, m, m)
    W_cumsum = zeros(T, m)

    for n in (i_start + 1):i_end
        dWn = W_grid[n] .- W_grid[n - 1]
        for k in 1:m
            for j in 1:m
                I[j, k] += W_cumsum[j] * dWn[k]
            end
        end
        W_cumsum .+= dWn
    end

    # Ito → Stratonovich
    for j in 1:m
        I[j, j] += dt / 2
    end
    return I
end

"""
Compute Stratonovich iterated integrals from a sequence of (dt_i, dW_i, dZ_i) tuples.
"""
function _iterated_I_from_increments(data, n_sub, m, T, dt)
    I = zeros(T, m, m)
    W_cumsum = zeros(T, m)

    for n in 1:n_sub
        dWn = data[n][2]  # dW_i from the tuple
        for k in 1:m
            for j in 1:m
                I[j, k] += W_cumsum[j] * dWn[k]
            end
        end
        W_cumsum .+= dWn
    end

    # Ito → Stratonovich
    for j in 1:m
        I[j, j] += dt / 2
    end
    return I
end

function _get_noise_source(W)
    if hasproperty(W, :source)
        return _get_noise_source(W.source)
    elseif W isa DiffEqNoiseProcess.NoiseGrid || W isa DiffEqNoiseProcess.NoiseProcess
        return W
    else
        return nothing
    end
end
