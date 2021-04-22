"""
From MATLAB ODE Suite by Shampine
"""
@def rosenbrock2332unpack begin
  if typeof(cache) <: OrdinaryDiffEqMutableCache
    d = cache.tab.d
  else
    d = cache.d
  end
end

@def rosenbrock2332pre0 begin
  @rosenbrock2332unpack
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs::Nothing,T::Type{Val{0}})
  @rosenbrock2332pre0
  @inbounds y₀ + dt*(c1*k[1] + c2*k[2])
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs::Nothing,T::Type{Val{0}})
  @rosenbrock2332pre0
  @inbounds @.. y₀ + dt*(c1*k[1] + c2*k[2])
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs,T::Type{Val{0}})
  @rosenbrock2332pre0
  @.. y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs::Nothing,T::Type{Val{0}})
  @rosenbrock2332pre0
  @inbounds @.. out = y₀ + dt*(c1*k[1] + c2*k[2])
  out
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs,T::Type{Val{0}})
  @rosenbrock2332pre0
  @views @.. out = y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
  out
end

# First Derivative of the dense output
@def rosenbrock2332pre1 begin
  @rosenbrock2332unpack
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs::Nothing,T::Type{Val{1}})
  @rosenbrock2332pre1
  @.. c1diff*k[1] + c2diff*k[2]
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs,T::Type{Val{1}})
  @rosenbrock2332pre1
  @.. c1diff*k[1][idxs] + c2diff*k[2][idxs]
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs::Nothing,T::Type{Val{1}})
  @rosenbrock2332pre1
  @.. out = c1diff*k[1] + c2diff*k[2]
  out
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock23Cache,Rosenbrock32ConstantCache,Rosenbrock32Cache},idxs,T::Type{Val{1}})
  @rosenbrock2332pre1
  @views @.. out = c1diff*k[1][idxs] + c2diff*k[2][idxs]
  out
end

"""
From MATLAB ODE Suite by Shampine
"""
@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Rodas4ConstantCache,idxs::Nothing,T::Type{Val{0}})
  Θ1 = 1 - Θ
  @inbounds Θ1*y₀ + Θ*(y₁ + Θ1*(k[1] + Θ*k[2]))
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Rodas4Cache,idxs::Nothing,T::Type{Val{0}})
  Θ1 = 1 - Θ
  @inbounds @.. Θ1*y₀ + Θ*(y₁ + Θ1*(k[1] + Θ*k[2]))
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rodas4ConstantCache,Rodas4Cache},idxs,T::Type{Val{0}})
  Θ1 = 1 - Θ
  @.. Θ1*y₀[idxs] + Θ*(y₁[idxs] + Θ1*(k[1][idxs] + Θ*k[2][idxs]))
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rodas4ConstantCache,Rodas4Cache},idxs::Nothing,T::Type{Val{0}})
  Θ1 = 1 - Θ
  @.. out = Θ1*y₀ + Θ*(y₁ + Θ1*(k[1] + Θ*k[2]))
  out
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rodas4ConstantCache,Rodas4Cache},idxs,T::Type{Val{0}})
  Θ1 = 1 - Θ
  @views @.. out = Θ1*y₀[idxs] + Θ*(y₁[idxs] + Θ1*(k[1][idxs] + Θ*k[2][idxs]))
  out
end

# First Derivative 
@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Rodas4ConstantCache,idxs::Nothing,T::Type{Val{1}})
  @inbounds (k[1] + Θ*(-2*k[1] + 2*k[2] - 3*k[2]*Θ) - y₀ + y₁)/dt
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Rodas4Cache,idxs::Nothing,T::Type{Val{1}})
  @inbounds @.. (k[1] + Θ*(-2*k[1] + 2*k[2] - 3*k[2]*Θ) - y₀ + y₁)/dt
end

@muladd function _ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rodas4ConstantCache,Rodas4Cache},idxs,T::Type{Val{1}})
  @.. (k[1][idxs] + Θ*(-2*k[1][idxs] + 2*k[2][idxs] - 3*k[2][idxs]*Θ) - y₀[idxs] + y₁[idxs])/dt
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rodas4ConstantCache,Rodas4Cache},idxs::Nothing,T::Type{Val{1}})
  @.. out = (k[1] + Θ*(-2*k[1] + 2*k[2] - 3*k[2]*Θ) - y₀ + y₁)/dt
  out
end

@muladd function _ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rodas4ConstantCache,Rodas4Cache},idxs,T::Type{Val{1}})
  @views @.. out = (k[1][idxs] + Θ*(-2*k[1][idxs] + 2*k[2][idxs] - 3*k[2][idxs]*Θ) - y₀[idxs] + y₁[idxs])/dt
  out
end
