struct DiffEqNLSolveTag end

struct DiffCache{T<:AbstractArray, S<:AbstractArray}
    du::T
    dual_du::S
end

Base.@pure function DiffCache(T, size, ::Type{Val{chunk_size}}) where chunk_size
  DiffCache(fill(zero(T), size...), fill(zero(Dual{typeof(ForwardDiff.Tag(DiffEqNLSolveTag(),T)),T,chunk_size}), size...))
end

Base.@pure DiffCache(u::AbstractArray) = DiffCache(eltype(u),size(u),Val{ForwardDiff.pickchunksize(length(u))})
Base.@pure DiffCache(u::AbstractArray,nlsolve) = DiffCache(eltype(u),size(u),Val{get_chunksize(nlsolve)})
Base.@pure DiffCache(u::AbstractArray,T::Type{Val{CS}}) where {CS} = DiffCache(eltype(u),size(u),T)

get_du(dc::DiffCache, ::Type{T}) where {T<:Dual} = dc.dual_du
get_du(dc::DiffCache, T) = dc.du

# Default nlsolve behavior, should move to DiffEqDiffTools.jl

Base.@pure determine_chunksize(u,alg::DiffEqBase.DEAlgorithm) = determine_chunksize(u,get_chunksize(alg))
Base.@pure function determine_chunksize(u,CS)
  if CS != 0
    return CS
  else
    return ForwardDiff.pickchunksize(length(u))
  end
end

struct NLSOLVEJL_SETUP{CS,AD} end
Base.@pure NLSOLVEJL_SETUP(;chunk_size=0,autodiff=true) = NLSOLVEJL_SETUP{chunk_size,autodiff}()
(::NLSOLVEJL_SETUP)(f,u0; kwargs...) = (res=NLsolve.nlsolve(f,u0; kwargs...); res.zero)
function (p::NLSOLVEJL_SETUP{CS,AD})(::Type{Val{:init}},f,u0_prototype) where {CS,AD}
  AD ? autodiff = :forward : autodiff = :central
  OnceDifferentiable(f, u0_prototype, u0_prototype, autodiff,
                     ForwardDiff.Chunk(determine_chunksize(u0_prototype,CS)))
end

get_chunksize(x) = 0
get_chunksize(x::NLSOLVEJL_SETUP{CS,AD}) where {CS,AD} = CS

@inline @muladd DiffEqBase.calculate_residuals(ũ::Number, u₀::Number, u₁::Number,
                                    α, ρ, internalnorm,t) = ũ / (α + max(internalnorm(u₀,t), internalnorm(u₁,t)) * ρ)

@inline @muladd DiffEqBase.calculate_residuals(u₀::Number, u₁::Number,
                                    α, ρ, internalnorm,t) = (u₁ - u₀) / (α + max(internalnorm(u₀,t), internalnorm(u₁,t)) * ρ)

"""
    DiffEqBase.calculate_residuals!(out, ũ, u₀, u₁, α, ρ)

Save element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
in `out`.
"""
@inline function DiffEqBase.calculate_residuals!(out, ũ, u₀, u₁, α, ρ, internalnorm,t)
  @.. out = DiffEqBase.calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm,t)
  nothing
end

@inline function DiffEqBase.calculate_residuals!(out::Array{T}, ũ::Array{T}, u₀::Array{T},
                                              u₁::Array{T}, α::T2, ρ::Real,
                                              internalnorm,t) where {T<:Number,T2<:Number}
  @tight_loop_macros for i in eachindex(out)
    @inbounds out[i] = DiffEqBase.calculate_residuals(ũ[i], u₀[i], u₁[i], α, ρ, internalnorm,t)
  end
  nothing
end

"""
    calculate_residuals!(out, u₀, u₁, α, ρ)

Save element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
in `out`.
"""
@inline function DiffEqBase.calculate_residuals!(out, u₀, u₁, α, ρ, internalnorm,t)
  @.. out = DiffEqBase.calculate_residuals(u₀, u₁, α, ρ, internalnorm,t)
end

@inline function DiffEqBase.calculate_residuals!(out::Array{T}, u₀::Array{T},
                                              u₁::Array{T}, α::T2, ρ::Real,
                                              internalnorm,t) where {T<:Number,T2<:Number}
  @tight_loop_macros for i in eachindex(out)
    @inbounds out[i] = DiffEqBase.calculate_residuals(u₀[i], u₁[i], α, ρ, internalnorm,t)
  end
end

"""
    DiffEqBase.calculate_residuals(ũ, u₀, u₁, α, ρ)

Calculate element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
"""
@inline function DiffEqBase.calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm,t)
  @.. DiffEqBase.calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm,t)
end

@inline function DiffEqBase.calculate_residuals(ũ::Array{T}, u₀::Array{T}, u₁::Array{T}, α::T2,
                                             ρ::Real, internalnorm,t) where
                                             {T<:Number,T2<:Number}
    out = similar(ũ)
    DiffEqBase.calculate_residuals!(out, ũ, u₀, u₁, α, ρ, internalnorm,t)
    out
end

"""
    calculate_residuals(u₀, u₁, α, ρ)

Calculate element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
"""
@inline function DiffEqBase.calculate_residuals(u₀, u₁, α, ρ, internalnorm,t)
  @.. DiffEqBase.calculate_residuals(u₀, u₁, α, ρ, internalnorm,t)
end

@inline function DiffEqBase.calculate_residuals(u₀::Array{T}, u₁::Array{T}, α::T2,
                                             ρ::Real, internalnorm,t) where
                                             {T<:Number,T2<:Number}
    out = similar(u₀)
    DiffEqBase.calculate_residuals!(out, u₀, u₁, α, ρ, internalnorm,t)
    out
end

macro swap!(x,y)
  quote
    local tmp = $(esc(x))
    $(esc(x)) = $(esc(y))
    $(esc(y)) = tmp
  end
end

islinear(f) = f isa DiffEqBase.AbstractDiffEqLinearOperator && is_constant(f)

macro cache(expr)
  name = expr.args[2].args[1].args[1]
  fields = expr.args[3].args[2:2:end]
  cache_vars = Expr[]
  jac_vars = Pair{Symbol,Expr}[]
  for x in fields
    if x.args[2] == :uType || x.args[2] == :rateType ||
       x.args[2] == :kType || x.args[2] == :uNoUnitsType
      push!(cache_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :JCType
      push!(cache_vars,:(c.$(x.args[1]).duals...))
    elseif x.args[2] == :GCType
      push!(cache_vars,:(c.$(x.args[1]).duals))
    elseif x.args[2] == :DiffCacheType
      push!(cache_vars,:(c.$(x.args[1]).du))
      push!(cache_vars,:(c.$(x.args[1]).dual_du))
    elseif x.args[2] == :JType || x.args[2] == :WType
      push!(jac_vars,x.args[1] => :(c.$(x.args[1])))
    end
  end
  quote
    $expr
    $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
    $(esc(:jac_iter))($(esc(:c))::$name) = tuple($(jac_vars...))
  end
end
