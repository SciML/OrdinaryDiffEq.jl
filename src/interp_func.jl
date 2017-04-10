immutable InterpolationData{F,uType,tType,kType,cacheType} <: Function
  f::F
  timeseries::uType
  ts::tType
  ks::kType
  notsaveat_idxs::Vector{Int}
  dense::Bool
  cache::cacheType
end

immutable CompositeInterpolationData{F,uType,tType,kType,cacheType} <: Function
  f::F
  timeseries::uType
  ts::tType
  ks::kType
  alg_choice::Vector{Int}
  notsaveat_idxs::Vector{Int}
  dense::Bool
  cache::cacheType
end

(interp::InterpolationData)(tvals,idxs,deriv) = (typeof(interp.cache) <: Union{DiscreteCache,DiscreteConstantCache}) || interp.dense ? ode_interpolation(tvals,interp,idxs,deriv) : error("Dense output must be enabled in order to interpolate. This requires save_everystep=true and dense=true.")
(interp::CompositeInterpolationData)(tvals,idxs,deriv) = (typeof(interp.cache) <: Union{DiscreteCache,DiscreteConstantCache}) || interp.dense ? ode_interpolation(tvals,interp,idxs,deriv) : error("Dense output must be enabled in order to interpolate. This requires save_everystep=true and dense=true.")
(interp::InterpolationData)(val,tvals,idxs,deriv) = (typeof(interp.cache) <: Union{DiscreteCache,DiscreteConstantCache}) || interp.dense ? ode_interpolation!(val,tvals,interp,idxs,deriv) : error("Dense output must be enabled in order to interpolate. This requires save_everystep=true and dense=true.")
(interp::CompositeInterpolationData)(val,tvals,idxs,deriv) = (typeof(interp.cache) <: Union{DiscreteCache,DiscreteConstantCache}) || interp.dense ? ode_interpolation!(val,tvals,interp,idxs,deriv) : error("Dense output must be enabled in order to interpolate. This requires save_everystep=true and dense=true.")

function InterpolationData(id::InterpolationData,f)
  InterpolationData(f,id.timeseries,
                      id.ts,
                      id.ks,
                      id.notsaveat_idxs,
                      id.dense,
                      id.cache)
end

function CompositeInterpolationData(id::CompositeInterpolationData,f)
  CompositeInterpolationData(f,id.timeseries,
                               id.ts,
                               id.ks,
                               id.alg_choice,
                               id.notsaveat_idxs,
                               id.dense,
                               id.cache)
end
