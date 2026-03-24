"""
    calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm, t)

Calculate element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
"""
@inline @muladd function calculate_residuals(
        ũ::Number, u₀::Number, u₁::Number,
        α, ρ, internalnorm, t
    )
    @fastmath ũ / (α + max(internalnorm(u₀, t), internalnorm(u₁, t)) * ρ)
end

@inline function calculate_residuals(
        ũ::Array{T}, u₀::Array{T}, u₁::Array{T}, α::T2,
        ρ::Real, internalnorm,
        t
    ) where
    {T <: Number, T2 <: Number}
    out = similar(ũ)
    calculate_residuals!(out, ũ, u₀, u₁, α, ρ, internalnorm, t)
    return out
end

@inline function calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm, t)
    return @.. broadcast = false calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm, t)
end

"""
    calculate_residuals(u₀, u₁, α, ρ, internalnorm, t)

Calculate element-wise residuals
```math
\\frac{u₁ - u₀}{α+\\max{|u₀|,|u₁|}*ρ}
```
"""

@inline @muladd function calculate_residuals(
        u₀::Number, u₁::Number,
        α, ρ, internalnorm, t
    )
    @fastmath (u₁ - u₀) / (α + max(internalnorm(u₀, t), internalnorm(u₁, t)) * ρ)
end

@inline function calculate_residuals(
        u₀::Array{T}, u₁::Array{T}, α::T2,
        ρ::Real, internalnorm,
        t
    ) where {T <: Number, T2 <: Number}
    out = similar(u₀)
    calculate_residuals!(out, u₀, u₁, α, ρ, internalnorm, t)
    return out
end

@inline function calculate_residuals(u₀, u₁, α, ρ, internalnorm, t)
    return @.. broadcast = false calculate_residuals(u₀, u₁, α, ρ, internalnorm, t)
end

"""
    calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)

Return element-wise residuals
```math
\\frac{δ E₁ + E₂}{α+\\max{scalarnorm(u₀),scalarnorm(u₁)}*ρ}.
```
"""
@inline @muladd function calculate_residuals(
        E₁::Number, E₂::Number, u₀::Number, u₁::Number,
        α::Real, ρ::Real, δ::Number, scalarnorm, t
    )
    @fastmath (δ * E₁ + E₂) / (α + max(scalarnorm(u₀, t), scalarnorm(u₁, t)) * ρ)
end

@inline function calculate_residuals(
        E₁::Array{<:Number}, E₂::Array{<:Number},
        u₀::Array{<:Number}, u₁::Array{<:Number}, α::Real,
        ρ::Real, δ::Number, scalarnorm, t
    )
    out = similar(u₀)
    calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
    return out
end

@inline function calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
    return @.. broadcast = false calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
end

# Inplace Versions

"""
    DiffEqBase.calculate_residuals!(out, ũ, u₀, u₁, α, ρ, thread=False())

Save element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
in `out`.

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = False()`, default)
or use multiple threads (`thread = True()`) when Julia is started
with multiple threads.
"""
@inline function calculate_residuals!(
        out, ũ, u₀, u₁, α, ρ, internalnorm, t,
        thread::Union{False, True} = False()
    )
    @.. broadcast = false thread = thread out = calculate_residuals(
        ũ, u₀, u₁, α, ρ, internalnorm,
        t
    )
    return nothing
end

@inline function calculate_residuals!(
        out::Array, ũ::Array, u₀::Array, u₁::Array, α::Number,
        ρ::Number, internalnorm::F, t, ::False
    ) where {F}
    @inbounds @simd ivdep for i in eachindex(out, ũ, u₀, u₁)
        out[i] = calculate_residuals(ũ[i], u₀[i], u₁[i], α, ρ, internalnorm, t)
    end
    return nothing
end

"""
    calculate_residuals!(out, u₀, u₁, α, ρ, thread=False())

Save element-wise residuals
```math
\\frac{u₁ - u₀}{α+\\max{|u₀|,|u₁|}*ρ}
```
in `out`.

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = False()`, default)
or use multiple threads (`thread = True()`) when Julia is started
with multiple threads.
"""
@inline function calculate_residuals!(
        out, u₀, u₁, α, ρ, internalnorm, t,
        thread::Union{False, True} = False()
    )
    return @.. broadcast = false thread = thread out = calculate_residuals(u₀, u₁, α, ρ, internalnorm, t)
end

"""
    calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, thread=False())

Calculate element-wise residuals
```math
\\frac{δ E₁ + E₂}{α+\\max{scalarnorm(u₀),scalarnorm(u₁)}*ρ}.
```

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = False()`, default)
or use multiple threads (`thread = True()`) when Julia is started
with multiple threads.
"""
@inline function calculate_residuals!(
        out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t,
        thread::Union{False, True} = False()
    )
    @.. broadcast = false thread = thread out = calculate_residuals(
        E₁, E₂, u₀, u₁, α, ρ, δ,
        scalarnorm, t
    )
    return out
end
