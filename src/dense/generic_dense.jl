## Integrator Dispatches

# Can get rid of an allocation here with a function
# get_tmp_arr(integrator.cache) which gives a pointer to some
# cache array which can be modified.

@inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(integrator,f=integrator.f,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.u,integrator.dt,f,integrator.cache,always_calc_begin,allow_calc_end,force_calc_end)
  else
    ode_addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.u,integrator.dt,f,integrator.cache.caches[integrator.cache.current],always_calc_begin,allow_calc_end,force_calc_end)
  end
end

@inline function ode_interpolant(Θ,integrator::DEIntegrator,idxs,deriv)
  ode_addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    val = ode_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache,idxs,deriv)
  else
    val = ode_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache.caches[integrator.cache.current],idxs,deriv)
  end
  val
end

@inline function ode_interpolant!(val,Θ,integrator::DEIntegrator,idxs,deriv)
  ode_addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_interpolant!(val,Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache,idxs,deriv)
  else
    ode_interpolant!(val,Θ,integrator.dt,integrator.uprev,integrator.u,integrator.k,integrator.cache.caches[integrator.cache.current],idxs,deriv)
  end
end

@inline function current_interpolant(t::Number,integrator::DEIntegrator,idxs,deriv)
  Θ = (t-integrator.tprev)/integrator.dt
  ode_interpolant(Θ,integrator,idxs,deriv)
end

@inline function current_interpolant(t,integrator::DEIntegrator,idxs,deriv)
  Θ = (t.-integrator.tprev)./integrator.dt
  [ode_interpolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_interpolant!(val,t::Number,integrator::DEIntegrator,idxs,deriv)
  Θ = (t-integrator.tprev)/integrator.dt
  ode_interpolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_interpolant!(val,t,integrator::DEIntegrator,idxs,deriv)
  Θ = (t.-integrator.tprev)./integrator.dt
  [ode_interpolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant(t::Number,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  ode_extrapolant(Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant!(val,t::Number,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t-integrator.tprev)/(integrator.t-integrator.tprev)
  ode_extrapolant!(val,Θ,integrator,idxs,deriv)
end

@inline function current_extrapolant(t::AbstractArray,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [ode_extrapolant(ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function current_extrapolant!(val,t,integrator::DEIntegrator,idxs=size(integrator.uprev),deriv=Val{0})
  Θ = (t.-integrator.tprev)./(integrator.t-integrator.tprev)
  [ode_extrapolant!(val,ϕ,integrator,idxs,deriv) for ϕ in Θ]
end

@inline function ode_extrapolant!(val,Θ,integrator::DEIntegrator,idxs,deriv)
  ode_addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_interpolant!(val,Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,integrator.cache,idxs,deriv)
  else
    ode_interpolant!(val,Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,integrator.cache.caches[integrator.cache.current],idxs,deriv)
  end
end

@inline function ode_extrapolant(Θ,integrator::DEIntegrator,idxs,deriv)
  ode_addsteps!(integrator)
  if !(typeof(integrator.cache) <: CompositeCache)
    ode_interpolant(Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,integrator.cache,idxs,deriv)
  else
    ode_interpolant(Θ,integrator.t-integrator.tprev,integrator.uprev2,integrator.uprev,integrator.k,integrator.cache.caches[integrator.cache.current],idxs,deriv)
  end
end

##

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function ode_interpolation(tvals,id,idxs,deriv)
  @unpack ts,timeseries,ks,f,notsaveat_idxs,cache = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals)
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  if idxs == nothing
    vals = Vector{eltype(timeseries)}(length(tvals))
  elseif typeof(idxs) <: Number
    vals = Vector{eltype(first(timeseries))}(length(tvals))
  else
    vals = Vector{Vector{eltype(first(timeseries))}}(length(tvals))
  end
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(ts[notsaveat_idxs[i:end]]),t,rev=tdir<0)+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[notsaveat_idxs[i]] == t
      if idxs == nothing
        vals[j] = timeseries[notsaveat_idxs[i]]
      else
        vals[j] = timeseries[notsaveat_idxs[i]][idxs]
      end
    elseif ts[notsaveat_idxs[i-1]] == t # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = timeseries[notsaveat_idxs[i-1]]
      else
        vals[j] = timeseries[notsaveat_idxs[i-1]][idxs]
      end
    else
      dt = ts[notsaveat_idxs[i]] - ts[notsaveat_idxs[i-1]]
      Θ = (t-ts[notsaveat_idxs[i-1]])/dt
      if idxs == nothing && eltype(timeseries) <: AbstractArray
        idxs_internal = size(timeseries[notsaveat_idxs[i-1]])
      else
        idxs_internal = idxs
      end
      if typeof(cache) <: CompositeCache
        ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache.caches[id.alg_choice[notsaveat_idxs[i-1]]]) # update the kcurrent
        vals[j] = ode_interpolant(Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],ks[i],cache.caches[id.alg_choice[notsaveat_idxs[i-1]]],idxs_internal,deriv)
      else
        ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache) # update the kcurrent
        (typeof(cache) <: (DiscreteCache) || typeof(cache) <: DiscreteConstantCache) ? kval = 0 : kval = ks[i]
        vals[j] = ode_interpolant(Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],kval,cache,idxs_internal,deriv)
      end
    end
  end
  vals
end

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function ode_interpolation!(vals,tvals,id,idxs,deriv)
  @unpack ts,timeseries,ks,f,notsaveat_idxs,cache = id
  tdir = sign(ts[end]-ts[1])
  idx = sortperm(tvals)
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(ts[notsaveat_idxs[i:end]]),t,rev=tdir<0)+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[notsaveat_idxs[i]] == t
      if idxs == nothing
        vals[j] = timeseries[notsaveat_idxs[i]]
      else
        vals[j] = timeseries[notsaveat_idxs[i]][idxs]
      end
    elseif ts[notsaveat_idxs[i-1]] == t # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = timeseries[notsaveat_idxs[i-1]]
      else
        vals[j] = timeseries[notsaveat_idxs[i-1]][idxs]
      end
    else
      dt = ts[notsaveat_idxs[i]] - ts[notsaveat_idxs[i-1]]
      Θ = (t-ts[notsaveat_idxs[i-1]])/dt
      if idxs == nothing && eltype(vals) <: AbstractArray
        idxs_internal = eachindex(vals[j])
      else
        idxs_internal = idxs
      end
      if typeof(cache) <: CompositeCache
        ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache.caches[id.alg_choice[notsaveat_idxs[i-1]]]) # update the kcurrent
        if eltype(timeseries) <: AbstractArray
          ode_interpolant!(vals[j],Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],ks[i],cache.caches[id.alg_choice[notsaveat_idxs[i-1]]],idxs_internal,deriv)
        else
          vals[j] = ode_interpolant(Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],ks[i],cache.caches[id.alg_choice[notsaveat_idxs[i-1]]],idxs_internal,deriv)
        end
      else
        ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache) # update the kcurrent
        if eltype(timeseries) <: AbstractArray
          ode_interpolant!(vals[j],Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],ks[i],cache,idxs_internal,deriv)
        else
          (typeof(cache) <: (DiscreteCache) || typeof(cache) <: DiscreteConstantCache) ? kval = 0 : kval = ks[i]
          vals[j] = ode_interpolant(Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],kval,cache,idxs_internal,deriv)
        end
      end
    end
  end
end

"""
ode_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function ode_interpolation(tval::Number,id,idxs,deriv)
  @unpack ts,timeseries,ks,f,notsaveat_idxs,cache = id
  tdir = sign(ts[end]-ts[1])
  @inbounds i = searchsortedfirst(@view(ts[notsaveat_idxs]),tval,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  @inbounds if ts[notsaveat_idxs[i]] == tval
    if idxs == nothing
      val = timeseries[notsaveat_idxs[i]]
    else
      val = timeseries[notsaveat_idxs[i]][idxs]
    end
  elseif ts[notsaveat_idxs[i-1]] == tval # Can happen if it's the first value!
    if idxs == nothing
      val = timeseries[notsaveat_idxs[i-1]]
    else
      val = timeseries[notsaveat_idxs[i-1]][idxs]
    end
  else
    dt = ts[notsaveat_idxs[i]] - ts[notsaveat_idxs[i-1]]
    Θ = (tval-ts[notsaveat_idxs[i-1]])/dt
    if idxs == nothing && eltype(timeseries) <: AbstractArray
      idxs_internal = size(timeseries[notsaveat_idxs[i-1]])
    else
      idxs_internal = idxs
    end
    if typeof(cache) <: CompositeCache
      ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache.caches[id.alg_choice[notsaveat_idxs[i-1]]]) # update the kcurrent
      val = ode_interpolant(Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],ks[i],cache.caches[id.alg_choice[notsaveat_idxs[i-1]]],idxs_internal,deriv)
    else
      ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache) # update the kcurrent
      (typeof(cache) <: (DiscreteCache) || typeof(cache) <: DiscreteConstantCache) ? kval = 0 : kval = ks[i]
      val = ode_interpolant(Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],kval,cache,idxs_internal,deriv)
    end
  end
  val
end

"""
ode_interpolation!(out,tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
@inline function ode_interpolation!(out,tval::Number,id,idxs,deriv)
  @unpack ts,timeseries,ks,f,notsaveat_idxs,cache = id
  tdir = sign(ts[end]-ts[1])
  @inbounds i = searchsortedfirst(@view(ts[notsaveat_idxs]),tval,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  @inbounds if ts[notsaveat_idxs[i]] == tval
    if idxs == nothing
      copy!(out,timeseries[notsaveat_idxs[i]])
    else
      copy!(out,timeseries[notsaveat_idxs[i]][idxs])
    end
  elseif ts[notsaveat_idxs[i-1]] == tval # Can happen if it's the first value!
    if idxs == nothing
      copy!(out,timeseries[notsaveat_idxs[i-1]])
    else
      copy!(out,timeseries[notsaveat_idxs[i-1]][idxs])
    end
  else
    dt = ts[notsaveat_idxs[i]] - ts[notsaveat_idxs[i-1]]
    Θ = (tval-ts[notsaveat_idxs[i-1]])/dt
    if idxs == nothing
      idxs_internal = eachindex(out)
    else
      idxs_internal = idxs
    end
    if typeof(cache) <: CompositeCache
      ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache.caches[id.alg_choice[notsaveat_idxs[i-1]]]) # update the kcurrent
      ode_interpolant!(out,Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],ks[i],cache.caches[id.alg_choice[notsaveat_idxs[i-1]]],idxs_internal,deriv)
    else
      ode_addsteps!(ks[i],ts[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],dt,f,cache) # update the kcurrent
      (typeof(cache) <: (DiscreteCache) || typeof(cache) <: DiscreteConstantCache) ? kval = 0 : kval = ks[i]
      ode_interpolant!(out,Θ,dt,timeseries[notsaveat_idxs[i-1]],timeseries[notsaveat_idxs[i]],kval,cache,idxs_internal,deriv)
    end
  end
end

"""
By default, simpledense
"""
function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<2 || calcVal
    if !(typeof(uprev)<:AbstractArray)
      copyat_or_push!(k,1,f(t,uprev))
      copyat_or_push!(k,2,f(t+dt,u))
    else
      rtmp = similar(cache.fsalfirst)
      f(t,uprev,rtmp)
      copyat_or_push!(k,1,rtmp)
      f(t+dt,u,rtmp)
      copyat_or_push!(k,2,rtmp,Val{false})
    end
  end
  nothing
end

@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::OrdinaryDiffEqMutableCache,idxs,T::Type{Val{0}})
  if typeof(idxs) <: Tuple
    out = similar(y₀,idxs)
    idxs_internal=eachindex(y₀)
  else
    out = similar(y₀,length(idxs))
    idxs_internal=idxs
  end
  ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs_internal,T)
  if typeof(idxs) <: Number
    return first(out)
  else
    return out
  end
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{0}}) # Default interpolant is Hermite
  if typeof(y₀) <: AbstractArray
    if typeof(idxs) <: Tuple
      out = similar(y₀,idxs)
    else
      out = similar(y₀,length(idxs))
    end
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*k[1][i] + Θ*dt*k[2][i])
    end
  else
    out = (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
  end
  out
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache,idxs,T::Type{Val{0}}) # Default interpolant is Hermite
  @inbounds for (j,i) in enumerate(idxs)
    out[j] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*k[1][i] + Θ*dt*k[2][i])
  end
end
