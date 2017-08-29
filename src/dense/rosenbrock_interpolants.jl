"""
From MATLAB ODE Suite by Shampine
"""
@muladd @inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs,T::Type{Val{0}})
  d = cache.d
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
  @. y₀ + dt*(c1*k[1] + c2*k[2])
end

# First Derivative of the dense output
@muladd @inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs,T::Type{Val{1}})
  d = cache.d
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
  @. c1diff*k[1] + c2diff*k[2]
end

"""
From MATLAB ODE Suite by Shampine
"""
@muladd @inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs,T::Type{Val{0}})
  d = cache.tab.d
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
  if out == nothing
    if idxs == nothing
      return @. y₀ + dt*(c1*k[1] + c2*k[2])
    else
      return @. y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
    end
  elseif idxs == nothing
    @. out = y₀ + dt*(c1*k[1] + c2*k[2])
  else
    @views @. out = y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
  end
end

@muladd @inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs,T::Type{Val{1}})
  d = cache.tab.d
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
  if out == nothing
    if idxs == nothing
      return @. c1diff*k[1] + c2diff*k[2]
    else
      return @. c1diff*k[1][idxs] + c2diff*k[2][idxs]
    end
  elseif idxs == nothing
    @. out = c1diff*k[1] + c2diff*k[2]
  else
    @views @. out = c1diff*k[1][idxs] + c2diff*k[2][idxs]
  end
end

@muladd @inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Rodas4ConstantCache,idxs,T::Type{Val{0}})
  Θ1 = 1 - Θ
  @. Θ1*y₀ + Θ*(y₁ + Θ1*(k[1] + Θ*k[2]))
end

"""
From MATLAB ODE Suite by Shampine
"""
@muladd @inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Rodas4Cache,idxs,T::Type{Val{0}})
  Θ1 = 1 - Θ
  if out == nothing
    if idxs == nothing
      return @. Θ1*y₀ + Θ*(y₁ + Θ1*(k[1] + Θ*k[2]))
    else
      return @. Θ1*y₀[idxs] + Θ*(y₁[idxs] + Θ1*(k[1][idxs] + Θ*k[2][idxs]))
    end
  elseif idxs == nothing
    @. out = Θ1*y₀ + Θ*(y₁ + Θ1*(k[1] + Θ*k[2]))
  else
    @views @. out = Θ1*y₀[idxs] + Θ*(y₁[idxs] + Θ1*(k[1][idxs] + Θ*k[2][idxs]))
  end
end
