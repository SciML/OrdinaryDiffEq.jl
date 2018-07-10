# bdf_utils

function U!(k, U, inplace)
  for j = 1:k
    for r = 1:k
      if inplace
        @. U[j,r] = inv(factorial(j)) * prod([m-r for m in 0:(j-1)])
      else
        U[j,r] = inv(factorial(j)) * prod([m-r for m in 0:(j-1)]) 
      end
    end
  end  
end

function R!(k, ρ, cache)
  @unpack R = cache
  for j = 1:k
    for r = 1:k
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        @. R[j,r] = inv(factorial(j)) * prod([m-r*ρ for m in 0:(j-1)])
      else
        R[j,r] = inv(factorial(j)) * prod([m-r*ρ for m in 0:(j-1)])
      end
    end
  end
end

function bdf1_EEst!(integrator, u, uprev, cache)
  @unpack t, dt = integrator
  uprev2 = integrator.uprev2
  tprev = integrator.tprev
  dt1 = dt*(t+dt-tprev)
  dt2 = (t-tprev)*(t+dt-tprev)
  c = 7/12
  r = c*dt^2
  tmp = r*abs((u - uprev)/dt1 - (uprev - uprev2)/dt2)
  atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
  integrator.EEst = integrator.opts.internalnorm(atmp)
  @show integrator.EEst
end