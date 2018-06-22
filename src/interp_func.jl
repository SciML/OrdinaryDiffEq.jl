abstract type OrdinaryDiffEqInterpolation{cacheType} <: AbstractDiffEqInterpolation end

struct InterpolationData{F,uType,tType,kType,cacheType} <: OrdinaryDiffEqInterpolation{cacheType}
  f::F
  timeseries::uType
  ts::tType
  ks::kType
  dense::Bool
  cache::cacheType
end

struct CompositeInterpolationData{F,uType,tType,kType,cacheType} <: OrdinaryDiffEqInterpolation{cacheType}
  f::F
  timeseries::uType
  ts::tType
  ks::kType
  alg_choice::Vector{Int}
  dense::Bool
  cache::cacheType
end

DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where {cacheType<:FunctionMapConstantCache} = "left-endpoint piecewise constant"
DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where {cacheType<:FunctionMapCache} = "left-endpoint piecewise constant"
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{DP5ConstantCache,DP5Cache,DP5ThreadedCache}
  interp.dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache,Rosenbrock23Cache,Rosenbrock32Cache}
  interp.dense ? "specialized 2nd order \"free\" stiffness-aware interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Rodas4ConstantCache,Rodas4Cache}
  interp.dense ? "specialized 3rd order \"free\" stiffness-aware interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{SSPRK22,SSPRK22ConstantCache,SSPRK33,SSPRK33ConstantCache,SSPRK432,SSPRK432ConstantCache}
  interp.dense ? "2nd order \"free\" SSP interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{OwrenZen3Cache,OwrenZen3ConstantCache}
  interp.dense ? "specialized 3rd order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{OwrenZen4Cache,OwrenZen4ConstantCache}
  interp.dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{OwrenZen5Cache,OwrenZen5ConstantCache}
  interp.dense ? "specialized 5th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Tsit5Cache,Tsit5ConstantCache}
  interp.dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{BS5ConstantCache,BS5Cache}
  interp.dense ? "specialized 5th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Vern6Cache,Vern6ConstantCache}
  interp.dense ? "specialized 6th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Vern7Cache,Vern7ConstantCache}
  interp.dense ? "specialized 7th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Vern8Cache,Vern8ConstantCache}
  interp.dense ? "specialized 8th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{Vern9Cache,Vern9ConstantCache}
  interp.dense ? "specialized 9th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType<:Union{DP8ConstantCache,DP8Cache}
  interp.dense ? "specialized 7th order interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{cacheType}) where cacheType
  interp.dense ? "3rd order Hermite" : "1st order linear"
end

(interp::InterpolationData)(tvals,idxs,deriv,p) = ode_interpolation(tvals,interp,idxs,deriv,p)
(interp::CompositeInterpolationData)(tvals,idxs,deriv,p) = ode_interpolation(tvals,interp,idxs,deriv,p)
(interp::InterpolationData)(val,tvals,idxs,deriv,p) = ode_interpolation!(val,tvals,interp,idxs,deriv,p)
(interp::CompositeInterpolationData)(val,tvals,idxs,deriv,p) = ode_interpolation!(val,tvals,interp,idxs,deriv,p)

function InterpolationData(id::InterpolationData,f)
  InterpolationData(f,id.timeseries,
                      id.ts,
                      id.ks,
                      id.dense,
                      id.cache)
end

function CompositeInterpolationData(id::CompositeInterpolationData,f)
  CompositeInterpolationData(f,id.timeseries,
                               id.ts,
                               id.ks,
                               id.alg_choice,
                               id.dense,
                               id.cache)
end
