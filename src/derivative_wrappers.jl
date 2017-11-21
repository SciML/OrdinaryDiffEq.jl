mutable struct TimeGradientWrapper{fType,uType} <: Function
  f::fType
  uprev::uType
  fx1::uType
end
(p::TimeGradientWrapper)(t) = (du2 = similar(p.uprev); p.f(t,p.uprev,du2); du2)
(p::TimeGradientWrapper)(du2,t) = p.f(t,p.uprev,du2)

mutable struct UJacobianWrapper{fType,tType,CacheType} <: Function
  f::fType
  t::tType
  x1::CacheType
  fx1::CacheType
end

(p::UJacobianWrapper)(du1,uprev) = p.f(p.t,uprev,du1)
(p::UJacobianWrapper)(uprev) = (du1 = similar(uprev); p.f(p.t,uprev,du1); du1)

mutable struct TimeDerivativeWrapper{F,uType} <: Function
  f::F
  u::uType
end
(p::TimeDerivativeWrapper)(t) = p.f(t,p.u)

mutable struct UDerivativeWrapper{F,tType} <: Function
  f::F
  t::tType
end
(p::UDerivativeWrapper)(u) = p.f(p.t,u)

function derivative!(df::AbstractArray{<:Number}, f, x::Union{Number,AbstractArray{<:Number}}, fx::AbstractArray{<:Number}, integrator::DEIntegrator)
    if alg_autodiff(integrator.alg)
        ForwardDiff.derivative!(df, f, fx, x)
    else
        RealOrComplex = eltype(integrator.u) <: Complex ? Val{:Complex} : Val{:Real}
        DiffEqDiffTools.finite_difference!(df, f, x, integrator.alg.diff_type, RealOrComplex, Val{:DiffEqDerivativeWrapper}, fx)
    end
    nothing
end

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DEIntegrator, jac_config)
    if alg_autodiff(integrator.alg)
      ForwardDiff.jacobian!(J, f, fx, x, jac_config)
    else
      RealOrComplex = eltype(integrator.u) <: Complex ? Val{:Complex} : Val{:Real}
      DiffEqDiffTools.finite_difference_jacobian!(J, f, x, integrator.alg.diff_type, RealOrComplex, Val{:JacobianWrapper}, fx)
    end
    nothing
end
