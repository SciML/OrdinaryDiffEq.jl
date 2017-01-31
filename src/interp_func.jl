immutable InterpolationData{F,uType,tType,kType,cacheType} <: Function
  f::F
  timeseries::uType
  ts::tType
  ks::kType
  notsaveat_idxs::Vector{Int}
  cache::cacheType
end

immutable CompositeInterpolationData{F,uType,tType,kType,cacheType} <: Function
  f::F
  timeseries::uType
  ts::tType
  ks::kType
  alg_choice::Vector{Int}
  notsaveat_idxs::Vector{Int}
  cache::cacheType
end

(interp::InterpolationData)(tvals) = ode_interpolation(tvals,interp)
(interp::CompositeInterpolationData)(tvals) = ode_interpolation(tvals,interp)
(interp::InterpolationData)(val,tvals) = ode_interpolation!(val,tvals,interp)
(interp::CompositeInterpolationData)(val,tvals) = ode_interpolation!(val,tvals,interp)

function InterpolationData(id::InterpolationData,f)
  InterpolationData(f,id.timeseries,
                      id.ts,
                      id.ks,
                      id.notsaveat_idxs,
                      id.cache)
end

function CompositeInterpolationData(id::CompositeInterpolationData,f)
  CompositeInterpolationData(f,id.timeseries,
                               id.ts,
                               id.ks,
                               id.alg_choice,
                               id.notsaveat_idxs,
                               id.cache)
end
