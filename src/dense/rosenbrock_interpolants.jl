"""
From MATLAB ODE Suite by Shampine
"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs,T::Type{Val{0}})
  d = cache.d
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
  #@. y₀ + dt*(c1*k[1] + c2*k[2])
  y₀ + dt*(c1*k[1] + c2*k[2])
end

# First Derivative of the dense output
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs,T::Type{Val{1}})
  d = cache.d
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
  #@. c1diff*k[1] + c2diff*k[2]
  c1diff*k[1] + c2diff*k[2]
end

"""
From MATLAB ODE Suite by Shampine
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs,T::Type{Val{0}})
  d = cache.tab.d
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
  if out == nothing
    return y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
  elseif idxs == nothing
    #@. out = y₀ + dt*(c1*k[1] + c2*k[2])
    @inbounds for i in eachindex(out)
      out[i] = y₀[i] + dt*(c1*k[1][i] + c2*k[2][i])
    end
  else
    #@views @. out = y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(c1*k[1][i] + c2*k[2][i])
    end
  end
end

@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs,T::Type{Val{1}})
  d = cache.tab.d
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
  if out == nothing
    return c1diff*k[1][idxs] + c2diff*k[2][idxs]
  elseif idxs == nothing
    #@. out = c1diff*k[1] + c2diff*k[2]
    @inbounds for i in eachindex(out)
      out[i] = c1diff*k[1][i] + c2diff*k[2][i]
    end
  else
    #@views @. out = c1diff*k[1][idxs] + c2diff*k[2][idxs]
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = c1diff*k[1][i] + c2diff*k[2][i]
    end
  end
end

@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Rodas4ConstantCache,idxs,T::Type{Val{0}})
  y₀*(1-Θ)+Θ*(y₁+(1-Θ)*(k[1] + Θ*k[2]))
end
