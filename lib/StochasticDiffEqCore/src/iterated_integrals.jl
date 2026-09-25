"""
    AbstractJ

Supertype for the objects that evaluate the iterated stochastic integrals `J` needed
by higher order SDE solvers.

All stochastic iterated integrals are written in the Stratonovich sense, as indicated
by the `J`. An `AbstractJ` is produced from a problem and an algorithm by
[`get_Jalg`](@ref) and is then evaluated with [`get_iterated_I`](@ref) (out-of-place)
or [`get_iterated_I!`](@ref) (in-place, storing the result in the object's `J` field).

Subtypes: [`AbstractJDiagonal`](@ref), [`AbstractJCommute`](@ref), and the Lévy area
algorithms of StochasticDiffEqLevyArea together with
[`IteratedIntegralAlgorithm_iip`](@ref).
"""
abstract type AbstractJ end

"""
    AbstractJDiagonal <: AbstractJ

Iterated-integral evaluators for diagonal (or scalar) noise, where the integrals
reduce to `1/2 ΔWᵢ²` componentwise and no Lévy area approximation is required.

Subtypes: [`JDiagonal_oop`](@ref), [`JDiagonal_iip`](@ref).
"""
abstract type AbstractJDiagonal <: AbstractJ end

"""
    AbstractJCommute <: AbstractJ

Iterated-integral evaluators for commutative noise, where the Lévy area terms cancel
and the integrals reduce to the outer product `1/2 ΔW ΔWᵀ`.

Subtypes: [`JCommute_oop`](@ref), [`JCommute_iip`](@ref).
"""
abstract type AbstractJCommute <: AbstractJ end

"""
    JDiagonal_oop()

Out-of-place iterated-integral evaluator for diagonal noise.

[`get_iterated_I`](@ref) returns a freshly allocated `1/2 .* ΔW .* ΔW`.
"""
struct JDiagonal_oop <: AbstractJDiagonal end

"""
    JDiagonal_iip(ΔW)

In-place iterated-integral evaluator for diagonal noise, preallocated to the shape of
the Brownian increment `ΔW`.

[`get_iterated_I!`](@ref) writes `1/2 * ΔW^2` into the `J` field.
"""
mutable struct JDiagonal_iip{JType} <: AbstractJDiagonal
    J::JType
    JDiagonal_iip(ΔW) = new{typeof(ΔW)}(false .* ΔW .* ΔW)
end

"""
    JCommute_oop()

Out-of-place iterated-integral evaluator for commutative noise.

[`get_iterated_I`](@ref) returns a freshly allocated `1/2 .* vec(ΔW) .* vec(ΔW)'`.
"""
struct JCommute_oop <: AbstractJCommute end

"""
    JCommute_iip(ΔW)

In-place iterated-integral evaluator for commutative noise, preallocated to a
`length(ΔW) × length(ΔW)` matrix.

[`get_iterated_I!`](@ref) writes `1/2 * vec(ΔW) * vec(ΔW)'` into the `J` field.
"""
mutable struct JCommute_iip{JType} <: AbstractJCommute
    J::JType
    function JCommute_iip(ΔW)
        J = false .* vec(ΔW) .* vec(ΔW)'
        return new{typeof(J)}(J)
    end
end

"""
    get_iterated_I(dt, dW, dZ, alg, p = nothing, c = 1, γ = 1 // 1)

Evaluate and return the Stratonovich iterated stochastic integrals for the step, out
of place.

## Arguments

  - `dt`: The step size.
  - `dW`, `dZ`: The Brownian increment and, where the method needs one, the auxiliary
    increment for the step.
  - `alg`: The [`AbstractJ`](@ref) evaluator, as produced by [`get_Jalg`](@ref).
  - `p`: Number of terms in the Lévy area series. `nothing` (the default) selects it
    automatically from the requested accuracy `ε = c * dt^(γ + 1/2)`.
  - `c`, `γ`: Constant and exponent of that accuracy target, normally the constant and
    strong order of the calling solver.

See [`get_iterated_I!`](@ref) for the in-place form.
"""
function get_iterated_I(dt, dW, dZ, alg::JDiagonal_oop, p = nothing, c = 1, γ = 1 // 1)
    return 1 // 2 .* dW .* dW
end

"""
    get_iterated_I!(dt, dW, dZ, alg, p = nothing, c = 1, γ = 1 // 1)

In-place form of [`get_iterated_I`](@ref): compute the Stratonovich iterated
stochastic integrals for the step and store them in `alg.J`, returning `nothing`.
"""
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

# algs from StochasticDiffEqLevyArea (originally LevyArea.jl), see e.g.
# https://github.com/stochastics-uni-luebeck/LevyArea.jl/blob/68c5cb08ab103b4dcd3178651f7a5dd9ce8c666d/src/milstein.jl#L25
function get_iterated_I(
        dt, dW, dZ, alg::StochasticDiffEqLevyArea.AbstractIteratedIntegralAlgorithm,
        p = nothing, c = 1, γ = 1 // 1
    )
    if isnothing(p)
        ε = c * dt^(γ + 1 // 2)
        p = terms_needed(length(dW), dt, ε, alg, MaxL2())
    end
    I = StochasticDiffEqLevyArea.levyarea(dW / √dt, p, alg)
    return I .= 1 // 2 * dW .* dW' .+ dt .* I
end

"""
    IteratedIntegralAlgorithm_iip(ΔW, levyalg)

In-place wrapper around a StochasticDiffEqLevyArea iterated-integral algorithm.

StochasticDiffEqLevyArea's algorithms are allocating; this wrapper pairs one of them
(`levyalg`) with a preallocated `length(ΔW) × length(ΔW)` buffer `J`, so that
[`get_iterated_I!`](@ref) can serve in-place solvers with general non-commutative
noise.
"""
mutable struct IteratedIntegralAlgorithm_iip{JType, levyalgType} <:
    StochasticDiffEqLevyArea.AbstractIteratedIntegralAlgorithm
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
    J .= StochasticDiffEqLevyArea.levyarea(dW / √dt, p, levyalg)
    J .= 1 // 2 * dW .* dW' .+ dt .* J
    return nothing
end

"""
    get_Jalg(ΔW, dt, prob, alg) -> AbstractJ

Select the iterated-integral evaluator to use for `prob` with `alg`.

The choice combines the algorithm's `ii_approx` field — an
[`IteratedIntegralApprox`](@ref) — with the problem's noise structure and whether it
is in-place:

  - [`IILevyArea`](@ref) uses a diagonal evaluator when the noise is diagonal or
    scalar, and otherwise the StochasticDiffEqLevyArea algorithm chosen by
    `optimal_algorithm(length(ΔW), dt)` (wrapped in
    [`IteratedIntegralAlgorithm_iip`](@ref) for in-place problems).
  - [`IICommutative`](@ref) uses a diagonal evaluator for diagonal/scalar noise and a
    commutative evaluator otherwise.
  - Any other `ii_approx` is taken to be a concrete Lévy area algorithm and is used
    directly.

A solver subpackage can override the choice for a specific algorithm by adding a
method, e.g. `get_Jalg(ΔW, dt, prob, alg::MySolver) = MronRoe()`.
"""
function get_Jalg(ΔW, dt, prob, alg)
    return if alg.ii_approx isa IILevyArea
        if isinplace(prob)
            if ΔW isa Number || is_diagonal_noise(prob)
                return JDiagonal_iip(ΔW)
            else
                # optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm=MaxL2())
                return IteratedIntegralAlgorithm_iip(ΔW, StochasticDiffEqLevyArea.optimal_algorithm(length(ΔW), dt))
            end
        else
            if ΔW isa Number || is_diagonal_noise(prob)
                return JDiagonal_oop()
            else
                return StochasticDiffEqLevyArea.optimal_algorithm(length(ΔW), dt)
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
