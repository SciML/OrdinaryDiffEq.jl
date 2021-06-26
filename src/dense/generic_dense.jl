## Integrator Dispatches

# Can get rid of an allocation here with a function
# get_tmp_arr(integrator.cache) which gives a pointer to some
# cache array which can be modified.

@inline function _searchsortedfirst(v::AbstractVector, x, lo::Integer, forward::Bool)
    u = oftype(lo, 1)
    lo = lo - u
    hi = length(v) + u
    @inbounds while lo < hi - u
        m = (lo + hi) >>> 1
        @inbounds if (forward && v[m] < x) || (!forward && v[m] > x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

@inline function _searchsortedlast(v::AbstractVector, x, lo::Integer, forward::Bool)
  u = oftype(lo, 1)
  lo = lo - u
  hi = length(v) + u
  @inbounds while lo < hi - u
      m = (lo + hi) >>> 1
      @inbounds if (forward && v[m] > x) || (!forward && v[m] < x)
          hi = m
      else
          lo = m
      end
  end
  return lo
end

@inline function _ode_addsteps!(integrator,f=integrator.f,always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if !(typeof(integrator.cache) <: CompositeCache)
    DiffEqBase.addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.u,
                  integrator.dt,f,integrator.p,integrator.cache,
                  always_calc_begin,allow_calc_end,force_calc_end)
  else
    DiffEqBase.addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.u,
                  integrator.dt,f,integrator.p,
                  @inbounds(integrator.cache.caches[integrator.cache.current]),
                  always_calc_begin,allow_calc_end,force_calc_end)
  end
end
@inline DiffEqBase.addsteps!(integrator::ODEIntegrator,args...) = _ode_addsteps!(integrator,args...)

@inline function ode_interpolant(Θ,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  DiffEqBase.addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    val = ode_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache,idxs,deriv)
  else
    val = composite_ode_interpolant(Θ,integrator,integrator.cache.caches,integrator.cache.current,idxs,deriv)
  end
  val
end

@generated function composite_ode_interpolant(Θ, integrator, caches::T, current, idxs, deriv) where {T <: Tuple}
  expr = Expr(:block)
  for i in 1:length(T.types)
    push!(expr.args, quote
      if $i == current
        return ode_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,caches[$i],idxs,deriv)
      end
    end)
  end
  push!(expr.args, quote
    throw("Cache $current is not available. There are only $(length(caches)) caches.")
  end)
  return expr
end

@inline function ode_interpolant!(val,Θ,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  DiffEqBase.addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_interpolant!(val,Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache,idxs,deriv)
  else
    ode_interpolant!(val,Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache.caches[integrator.cache.current],idxs,deriv)
  end
end

@generated function composite_ode_interpolant!(val, Θ, integrator, caches::T, current, idxs, deriv) where {T <: Tuple}
  expr = Expr(:block)
  for i in 1:length(T.types)
    push!(expr.args, quote
      if $i == current
        return ode_interpolant!(val,Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,caches[$i],idxs,deriv)
      end
    end)
  end
  push!(expr.args, quote
    throw("Cache $current is not available. There are only $(length(caches)) caches.")
  end)
  return expr
end

@inline function current_interpolant(t::Number,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  Θ = (t-integrator.tprev)/integrator.dt
  ode_interpolant(Θ,integrator,idxs,deriv)
end

@inline function current_interpolant(t,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  Θ = (t.-integrator.tprev)./integrator.dt
  [ode_interpolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_interpolant!(val,t::Number,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  Θ = (t-integrator.tprev)/integrator.dt
  ode_interpolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_interpolant!(val,t,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  Θ = (t.-integrator.tprev)./integrator.dt
  [ode_interpolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant(t::Number,integrator::DiffEqBase.DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  ode_extrapolant(Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant!(val,t::Number,integrator::DiffEqBase.DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  ode_extrapolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant(t::AbstractArray,integrator::DiffEqBase.DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [ode_extrapolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant!(val,t,integrator::DiffEqBase.DEIntegrator,idxs=nothing,deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [ode_extrapolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function ode_extrapolant!(val,Θ,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  DiffEqBase.addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_interpolant!(val,Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,integrator.cache,idxs,deriv)
  else
    composite_ode_extrapolant!(val,Θ,integrator,integrator.cache.caches,integrator.cache.current,idxs,deriv)
  end
end

@generated function composite_ode_extrapolant!(val, Θ, integrator, caches::T, current, idxs, deriv) where {T <: Tuple}
  expr = Expr(:block)
  for i in 1:length(T.types)
    push!(expr.args, quote
      if $i == current
        return ode_interpolant!(val,Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,caches[$i],idxs,deriv)
      end
    end)
  end
  push!(expr.args, quote
    throw("Cache $current is not available. There are only $(length(caches)) caches.")
  end)
  return expr
end


@inline function ode_extrapolant(Θ,integrator::DiffEqBase.DEIntegrator,idxs,deriv)
  DiffEqBase.addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_interpolant(Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,integrator.cache,idxs,deriv)
  else
    composite_ode_extrapolant(Θ,integrator,integrator.cache.caches,integrator.cache.current,idxs,deriv)
  end
end

@generated function composite_ode_extrapolant(Θ, integrator, caches::T, current, idxs, deriv) where {T <: Tuple}
  expr = Expr(:block)
  for i in 1:length(T.types)
    push!(expr.args, quote
      if $i == current
        return ode_interpolant(Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,caches[$i],idxs,deriv)
      end
    end)
  end
  push!(expr.args, quote
    throw("Cache $current is not available. There are only $(length(caches)) caches.")
  end)
  return expr
end

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation(tvals,id::I,idxs,deriv::D,p,continuity::Symbol=:left) where {I,D}
  @unpack ts,timeseries,ks,f,cache = id
  @inbounds tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)

  # start the search thinking it's ts[1]-ts[2]
  i₋ = 1
  i₊ = 2
  vals = map(idx) do j
    t = tvals[j]

    if continuity === :left
      # we have i₋ = i₊ = 1 if t = ts[1], i₊ = i₋ + 1 = lastindex(ts) if t > ts[end],
      # and otherwise i₋ and i₊ satisfy ts[i₋] < t ≤ ts[i₊]
      i₊ = min(lastindex(ts), _searchsortedfirst(ts,t,i₊,tdir > 0))
      i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
      # we have i₋ = i₊ - 1 = 1 if t < ts[1], i₊ = i₋ = lastindex(ts) if t = ts[end],
      # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ t < ts[i₊]
      i₋ = max(1, _searchsortedlast(ts,t,i₋,tdir > 0))
      i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end

    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t-ts[i₋]) / dt

    if typeof(cache) <: (FunctionMapCache) || typeof(cache) <: FunctionMapConstantCache
      return ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],0,cache,idxs,deriv)
    elseif !id.dense
      return linear_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
    elseif typeof(cache) <: CompositeCache
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache.caches[id.alg_choice[i₊]]) # update the kcurrent
      return ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache.caches[id.alg_choice[i₊]],idxs,deriv)
    else
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache) # update the kcurrent
      return ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache,idxs,deriv)
    end
  end
  invpermute!(vals, idx)
  DiffEqArray(vals, tvals)
end

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation!(vals,tvals,id::I,idxs,deriv::D,p,continuity::Symbol=:left) where {I,D}
  @unpack ts,timeseries,ks,f,cache = id
  @inbounds tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals,rev=tdir<0)

  # start the search thinking it's in ts[1]-ts[2]
  i₋ = 1
  i₊ = 2
  @inbounds for j in idx
    t = tvals[j]

    if continuity === :left
      # we have i₋ = i₊ = 1 if t = ts[1], i₊ = i₋ + 1 = lastindex(ts) if t > ts[end],
      # and otherwise i₋ and i₊ satisfy ts[i₋] < t ≤ ts[i₊]
      i₊ = min(lastindex(ts), _searchsortedfirst(ts,t,i₊,tdir > 0))
      i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
      # we have i₋ = i₊ - 1 = 1 if t < ts[1], i₊ = i₋ = lastindex(ts) if t = ts[end],
      # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ t < ts[i₊]
      i₋ = max(1, _searchsortedlast(ts,t,i₋,tdir > 0))
      i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end

    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t-ts[i₋]) / dt

    if typeof(cache) <: (FunctionMapCache) || typeof(cache) <: FunctionMapConstantCache
      if eltype(timeseries) <: AbstractArray
        ode_interpolant!(vals[j],Θ,dt,timeseries[i₋],timeseries[i₊],0,cache,idxs,deriv)
      else
        vals[j] = ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],0,cache,idxs,deriv)
      end
    elseif !id.dense
      if eltype(timeseries) <: AbstractArray
        linear_interpolant!(vals[j],Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
      else
        vals[j] = linear_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
      end
    elseif typeof(cache) <: CompositeCache
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache.caches[id.alg_choice[i₊]]) # update the kcurrent
      if eltype(timeseries) <: AbstractArray
        ode_interpolant!(vals[j],Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache.caches[id.alg_choice[i₊]],idxs,deriv)
      else
        vals[j] = ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache.caches[id.alg_choice[i₊]],idxs,deriv)
      end
    else
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache) # update the kcurrent
      if eltype(vals[j]) <: AbstractArray
        ode_interpolant!(vals[j],Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache,idxs,deriv)
      else
        vals[j] = ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache,idxs,deriv)
      end
    end
  end

  vals
end

"""
ode_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation(tval::Number,id::I,idxs,deriv::D,p,continuity::Symbol=:left) where {I,D}
  @unpack ts,timeseries,ks,f,cache = id
  @inbounds tdir = sign(ts[end]-ts[1])

  if continuity === :left
    # we have i₋ = i₊ = 1 if tval = ts[1], i₊ = i₋ + 1 = lastindex(ts) if tval > ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] < tval ≤ ts[i₊]
    i₊ = min(lastindex(ts), _searchsortedfirst(ts,tval,2,tdir > 0))
    i₋ = i₊ > 1 ? i₊ - 1 : i₊
  else
    # we have i₋ = i₊ - 1 = 1 if tval < ts[1], i₊ = i₋ = lastindex(ts) if tval = ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ tval < ts[i₊]
    i₋ = max(1, _searchsortedlast(ts,tval,1,tdir > 0))
    i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
  end

  @inbounds begin
    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval-ts[i₋]) / dt

    if typeof(cache) <: (FunctionMapCache) || typeof(cache) <: FunctionMapConstantCache
      val = ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],0,cache,idxs,deriv)
    elseif !id.dense
      val = linear_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
    elseif typeof(cache) <: CompositeCache
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache.caches[id.alg_choice[i₊]]) # update the kcurrent
      val = ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache.caches[id.alg_choice[i₊]],idxs,deriv)
    else
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache) # update the kcurrent
      val = ode_interpolant(Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache,idxs,deriv)
    end
  end

  val
end

"""
ode_interpolation!(out,tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation!(out,tval::Number,id::I,idxs,deriv::D,p,continuity::Symbol=:left) where {I,D}
  @unpack ts,timeseries,ks,f,cache = id
  @inbounds tdir = sign(ts[end]-ts[1])

  if continuity === :left
    # we have i₋ = i₊ = 1 if tval = ts[1], i₊ = i₋ + 1 = lastindex(ts) if tval > ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] < tval ≤ ts[i₊]
    i₊ = min(lastindex(ts), _searchsortedfirst(ts,tval,2,tdir > 0))
    i₋ = i₊ > 1 ? i₊ - 1 : i₊
  else
    # we have i₋ = i₊ - 1 = 1 if tval < ts[1], i₊ = i₋ = lastindex(ts) if tval = ts[end],
    # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ tval < ts[i₊]
    i₋ = max(1, _searchsortedlast(ts,tval,1,tdir > 0))
    i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
  end

  @inbounds begin
    dt = ts[i₊] - ts[i₋]
    Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval-ts[i₋]) / dt

    if typeof(cache) <: (FunctionMapCache) || typeof(cache) <: FunctionMapConstantCache
      ode_interpolant!(out,Θ,dt,timeseries[i₋],timeseries[i₊],0,cache,idxs,deriv)
    elseif !id.dense
      linear_interpolant!(out,Θ,dt,timeseries[i₋],timeseries[i₊],idxs,deriv)
    elseif typeof(cache) <: CompositeCache
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache.caches[id.alg_choice[i₊]]) # update the kcurrent
      ode_interpolant!(out,Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache.caches[id.alg_choice[i₊]],idxs,deriv)
    else
      DiffEqBase.addsteps!(ks[i₊],ts[i₋],timeseries[i₋],timeseries[i₊],dt,f,p,cache) # update the kcurrent
      ode_interpolant!(out,Θ,dt,timeseries[i₋],timeseries[i₊],ks[i₊],cache,idxs,deriv)
    end
  end

  out
end

"""
By default, Hermite interpolant so update the derivative at the two ends
"""
function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache,always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin
    if typeof(cache) <: OrdinaryDiffEqMutableCache
      rtmp = similar(u, eltype(eltype(k)))
      f(rtmp,uprev,p,t)
      copyat_or_push!(k,1,rtmp)
      f(rtmp,u,p,t+dt)
      copyat_or_push!(k,2,rtmp)
    else
      copyat_or_push!(k,1,f(uprev,p,t))
      copyat_or_push!(k,2,f(u,p,t+dt))
    end
  end
  nothing
end

"""
ode_interpolant and ode_interpolant! dispatch
"""
function ode_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{TI}}) where TI
  _ode_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T)
end

function ode_interpolant(Θ,dt,y₀,y₁,k,cache::OrdinaryDiffEqMutableCache,idxs,T::Type{Val{TI}}) where TI
  if typeof(idxs) <: Number || typeof(y₀) <: Union{Number,SArray}
    # typeof(y₀) can be these if saveidxs gives a single value
    _ode_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T)
  elseif typeof(idxs) <: Nothing
    out = oneunit(Θ) .* y₁
    _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs,T)
  else
    out = oneunit(Θ) .* y₁[idxs]
    _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs,T)
  end
end

function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{TI}}) where TI
  _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs,T)
end

##################### Hermite Interpolants

# If no dispatch found, assume Hermite
function _ode_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{TI}}) where TI
  hermite_interpolant(Θ,dt,y₀,y₁,k,Val{typeof(cache) <: OrdinaryDiffEqMutableCache},idxs,T)
end

function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{TI}}) where TI
  hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs,T)
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{false}},idxs::Nothing,T::Type{Val{0}}) # Default interpolant is Hermite
  #@.. (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
  @inbounds (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{true}},idxs::Nothing,T::Type{Val{0}}) # Default interpolant is Hermite
  #@.. (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
  @inbounds @.. (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{0}}) # Default interpolant is Hermite
  # return @.. (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+(Θ-1)*dt*k[1][idxs] + Θ*dt*k[2][idxs])
  return (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+(Θ-1)*dt*k[1][idxs] + Θ*dt*k[2][idxs])
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs::Nothing,T::Type{Val{0}}) # Default interpolant is Hermite
  @inbounds @.. out = (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
  #@inbounds for i in eachindex(out)
  #  out[i] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*k[1][i] + Θ*dt*k[2][i])
  #end
  #out
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs,T::Type{Val{0}}) # Default interpolant is Hermite
  @views @.. out = (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+(Θ-1)*dt*k[1][idxs] + Θ*dt*k[2][idxs])
  #@inbounds for (j,i) in enumerate(idxs)
  #  out[j] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*k[1][i] + Θ*dt*k[2][i])
  #end
  #out
end

"""
Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{false}},idxs::Nothing,T::Type{Val{1}}) # Default interpolant is Hermite
  #@.. k[1] + Θ*(-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(3*dt*k[1] + 3*dt*k[2] + 6*y₀ - 6*y₁) + 6*y₁)/dt
  @inbounds k[1] + Θ*(-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(3*dt*k[1] + 3*dt*k[2] + 6*y₀ - 6*y₁) + 6*y₁)/dt
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{true}},idxs::Nothing,T::Type{Val{1}}) # Default interpolant is Hermite
  @inbounds @.. k[1] + Θ*(-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(3*dt*k[1] + 3*dt*k[2] + 6*y₀ - 6*y₁) + 6*y₁)/dt
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{1}}) # Default interpolant is Hermite
  # return @.. k[1][idxs] + Θ*(-4*dt*k[1][idxs] - 2*dt*k[2][idxs] - 6*y₀[idxs] + Θ*(3*dt*k[1][idxs] + 3*dt*k[2][idxs] + 6*y₀[idxs] - 6*y₁[idxs]) + 6*y₁[idxs])/dt
  return k[1][idxs] + Θ*(-4*dt*k[1][idxs] - 2*dt*k[2][idxs] - 6*y₀[idxs] + Θ*(3*dt*k[1][idxs] + 3*dt*k[2][idxs] + 6*y₀[idxs] - 6*y₁[idxs]) + 6*y₁[idxs])/dt
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs::Nothing,T::Type{Val{1}}) # Default interpolant is Hermite
  @inbounds @.. out = k[1] + Θ*(-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(3*dt*k[1] + 3*dt*k[2] + 6*y₀ - 6*y₁) + 6*y₁)/dt
  #@inbounds for i in eachindex(out)
  #  out[i] = k[1][i] + Θ*(-4*dt*k[1][i] - 2*dt*k[2][i] - 6*y₀[i] + Θ*(3*dt*k[1][i] + 3*dt*k[2][i] + 6*y₀[i] - 6*y₁[i]) + 6*y₁[i])/dt
  #end
  #out
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs,T::Type{Val{1}}) # Default interpolant is Hermite
  @views @.. out = k[1][idxs] + Θ*(-4*dt*k[1][idxs] - 2*dt*k[2][idxs] - 6*y₀[idxs] + Θ*(3*dt*k[1][idxs] + 3*dt*k[2][idxs] + 6*y₀[idxs] - 6*y₁[idxs]) + 6*y₁[idxs])/dt
  #@inbounds for (j,i) in enumerate(idxs)
  #  out[j] = k[1][i] + Θ*(-4*dt*k[1][i] - 2*dt*k[2][i] - 6*y₀[i] + Θ*(3*dt*k[1][i] + 3*dt*k[2][i] + 6*y₀[i] - 6*y₁[i]) + 6*y₁[i])/dt
  #end
  #out
end

"""
Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{false}},idxs::Nothing,T::Type{Val{2}}) # Default interpolant is Hermite
  #@.. (-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
  @inbounds (-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{true}},idxs::Nothing,T::Type{Val{2}}) # Default interpolant is Hermite
  #@.. (-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
  @inbounds @.. (-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{2}}) # Default interpolant is Hermite
  #out = similar(y₀,axes(idxs))
  #@views @.. out = (-4*dt*k[1][idxs] - 2*dt*k[2][idxs] - 6*y₀[idxs] + Θ*(6*dt*k[1][idxs] + 6*dt*k[2][idxs] + 12*y₀[idxs] - 12*y₁[idxs]) + 6*y₁[idxs])/(dt*dt)
  @views out = (-4*dt*k[1][idxs] - 2*dt*k[2][idxs] - 6*y₀[idxs] + Θ*(6*dt*k[1][idxs] + 6*dt*k[2][idxs] + 12*y₀[idxs] - 12*y₁[idxs]) + 6*y₁[idxs])/(dt*dt)
  out
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs::Nothing,T::Type{Val{2}}) # Default interpolant is Hermite
  @inbounds @.. out = (-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
  #@inbounds for i in eachindex(out)
  #  out[i] = (-4*dt*k[1][i] - 2*dt*k[2][i] - 6*y₀[i] + Θ*(6*dt*k[1][i] + 6*dt*k[2][i] + 12*y₀[i] - 12*y₁[i]) + 6*y₁[i])/(dt*dt)
  #end
  #out
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs,T::Type{Val{2}}) # Default interpolant is Hermite
  @views @.. out = (-4*dt*k[1][idxs] - 2*dt*k[2][idxs] - 6*y₀[idxs] + Θ*(6*dt*k[1][idxs] + 6*dt*k[2][idxs] + 12*y₀[idxs] - 12*y₁[idxs]) + 6*y₁[idxs])/(dt*dt)
  #@inbounds for (j,i) in enumerate(idxs)
  #  out[j] = (-4*dt*k[1][i] - 2*dt*k[2][i] - 6*y₀[i] + Θ*(6*dt*k[1][i] + 6*dt*k[2][i] + 12*y₀[i] - 12*y₁[i]) + 6*y₁[i])/(dt*dt)
  #end
  #out
end

"""
Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{false}},idxs::Nothing,T::Type{Val{3}}) # Default interpolant is Hermite
  #@.. (6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁)/(dt*dt*dt)
  @inbounds (6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁)/(dt*dt*dt)
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,::Type{Val{true}},idxs::Nothing,T::Type{Val{3}}) # Default interpolant is Hermite
  #@.. (6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁)/(dt*dt*dt)
  @inbounds @.. (6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁)/(dt*dt*dt)
end

@muladd function hermite_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{3}}) # Default interpolant is Hermite
  #out = similar(y₀,axes(idxs))
  #@views @.. out = (6*dt*k[1][idxs] + 6*dt*k[2][idxs] + 12*y₀[idxs] - 12*y₁[idxs])/(dt*dt*dt)
  @views out = (6*dt*k[1][idxs] + 6*dt*k[2][idxs] + 12*y₀[idxs] - 12*y₁[idxs])/(dt*dt*dt)
  out
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs::Nothing,T::Type{Val{3}}) # Default interpolant is Hermite
  @inbounds @.. out = (6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁)/(dt*dt*dt)
  #for i in eachindex(out)
  #  out[i] = (6*dt*k[1][i] + 6*dt*k[2][i] + 12*y₀[i] - 12*y₁[i])/(dt*dt*dt)
  #end
  #out
end

@muladd function hermite_interpolant!(out,Θ,dt,y₀,y₁,k,idxs,T::Type{Val{3}}) # Default interpolant is Hermite
  @views @.. out = (6*dt*k[1][idxs] + 6*dt*k[2][idxs] + 12*y₀[idxs] - 12*y₁[idxs])/(dt*dt*dt)
  #for (j,i) in enumerate(idxs)
  #  out[j] = (6*dt*k[1][i] + 6*dt*k[2][i] + 12*y₀[i] - 12*y₁[i])/(dt*dt*dt)
  #end
  #out
end

######################## Linear Interpolants

@muladd @inline function linear_interpolant(Θ,dt,y₀,y₁,idxs::Nothing,T::Type{Val{0}})
  Θm1 = (1-Θ)
  @.. Θm1*y₀ + Θ*y₁
end

@muladd @inline function linear_interpolant(Θ,dt,y₀,y₁,idxs,T::Type{Val{0}})
  Θm1 = (1-Θ)
  @.. Θm1*y₀[idxs] + Θ*y₁[idxs]
end

@muladd @inline function linear_interpolant!(out,Θ,dt,y₀,y₁,idxs::Nothing,T::Type{Val{0}})
  Θm1 = (1-Θ)
  @.. out = Θm1*y₀ + Θ*y₁
  out
end

@muladd @inline function linear_interpolant!(out,Θ,dt,y₀,y₁,idxs,T::Type{Val{0}})
  Θm1 = (1-Θ)
  @views @.. out = Θm1*y₀[idxs] + Θ*y₁[idxs]
  out
end

"""
Linear Interpolation
"""
@inline function linear_interpolant(Θ,dt,y₀,y₁,idxs::Nothing,T::Type{Val{1}})
  (y₁ - y₀)/dt
end

@inline function linear_interpolant(Θ,dt,y₀,y₁,idxs,T::Type{Val{1}})
  @.. (y₁[idxs] - y₀[idxs])/dt
end

@inline function linear_interpolant!(out,Θ,dt,y₀,y₁,idxs::Nothing,T::Type{Val{1}})
  @.. out = (y₁ - y₀)/dt
  out
end

@inline function linear_interpolant!(out,Θ,dt,y₀,y₁,idxs,T::Type{Val{1}})
  @views @.. out = (y₁[idxs] - y₀[idxs])/dt
  out
end
