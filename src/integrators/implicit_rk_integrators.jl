@inline function initialize!(integrator,cache::ImplicitEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::ImplicitEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack uf = cache
  uf.t = t
  if integrator.iter > 1 && !integrator.u_modified
    u = current_extrapolant(t+dt,integrator)
  else
    u = uprev
  end

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end

  z = u - uprev
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  b = -z + dt*f(t+dt,uprev + z)
  dz = W\b
  ndz = abs(dz)
  z = z + dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.iter > 1
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  while do_newton
    iter += 1
    b = -z + dt*f(t+dt,uprev + z)
    dz = W\b
    ndzprev = ndz
    ndz = abs(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev + z

  if integrator.opts.adaptive && integrator.iter > 1
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    DD3 = ((u - uprev)/((dt)*(t+dt-tprev)) + (uprev-uprev2)/((t-tprev)*(t+dt-tprev)))
    dEst = (dt^2)*abs(DD3/6)
    integrator.EEst = dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
  else
    integrator.EEst = 1
  end

  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end#

@inline function initialize!(integrator,cache::ImplicitEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@inline function perform_step!(integrator,cache::ImplicitEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  if integrator.iter > 1 && !integrator.u_modified
    current_extrapolant!(u,t+dt,integrator)
  else
    copy!(u,uprev)
  end

  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if cache.newton_iters == 1 && cache.ηold < 1e-3
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    if integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps()
      new_W = true
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = @muladd mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = u - uprev
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  f(t+dt,u,k)
  scale!(k,dt)
  k .-= z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz
  @. u = uprev + z

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.iter > 1
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  while do_newton
    iter += 1
    f(t+dt,u,k)
    scale!(k,dt)
    k .-= z
    if has_invW(f)
      A_mul_B!(dz,W,k) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z .+= dz
    @. u = uprev + z
  end

  cache.ηold = η
  cache.newton_iters = iter
  @. u = uprev + z

  if integrator.opts.adaptive && integrator.iter > 1
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    dt1 = (dt)*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds DD3 = (u[i] - uprev[i])/dt1 + (uprev[i]-uprev2[i])/dt2
      dEst = (dt^2)*abs(DD3)/6
      @inbounds k[i] = dEst/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(k)
  else
    integrator.EEst = 1
  end

  f(t+dt,u,integrator.fsallast)
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack uf = cache
  uf.t = t
  dto2 = dt/2

  if integrator.iter > 1 && !integrator.u_modified
    u = current_extrapolant(t+dt,integrator)
  else
    u = uprev
  end
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dto2*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dto2*J
  end
  z = u - uprev
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  b = -z + dto2*f(t+dt,uprev + z)
  dz = W\b
  ndz = abs(dz)
  z = z + dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.iter > 1
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  while do_newton
    iter += 1
    b = -z + dto2*f(t+dto2,uprev + z)
    dz = W\b
    ndzprev = ndz
    ndz = abs(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev + 2z
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast

  if integrator.opts.adaptive
    if integrator.iter > 2
      # Use 3rd divided differences a la SPICE and Shampine
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2
      DD31 = ((u - uprev)/((dt)*(t+dt-tprev)) + (uprev-uprev2)/((t-tprev)*(t+dt-tprev)))
      DD30 = ((uprev - uprev2)/((t-tprev)*(t-tprev2)) + (uprev2-uprev3)/((tprev-tprev2)*(t-tprev2)))
      dEst = (dt^3)*abs(((DD31 - DD30)/(t+dt-tprev2))/12)
      integrator.EEst = dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      if integrator.EEst <= 1
        cache.uprev3 = uprev2
        cache.tprev2 = tprev
      end
    elseif integrator.iter > 1
      integrator.EEst = 1
      cache.uprev3 = integrator.uprev2
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::TrapezoidCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@inline function perform_step!(integrator,cache::TrapezoidCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  if integrator.iter > 1 && !integrator.u_modified
    current_extrapolant!(u,t+dt,integrator)
  else
    copy!(u,uprev)
  end
  dto2 = dt/2
  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if cache.newton_iters == 1 && cache.ηold < 1e-3
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    if integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps()
      new_W = true
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = @muladd mass_matrix[i,j]-dto2*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = u - uprev
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  f(t+dto2,u,k)
  scale!(k,dto2)
  k .-= z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz
  @. u = uprev + z

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.iter > 1
    @show η*ndz,κ*tol
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  @show "thisdz",ndz

  while do_newton
    iter += 1
    f(t+dto2,u,k)
    scale!(k,dto2)
    k .-= z
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    @show "seconddz",ndz
    θ = ndz/ndzprev
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z .+= dz
    @. u = uprev + z
  end

  cache.ηold = η
  cache.newton_iters = iter
  @. u = uprev + 2*z

  if integrator.opts.adaptive
    if integrator.iter > 2
      # Use 3rd divided differences a la SPICE and Shampine
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2
      dt1 = (dt)*(t+dt-tprev)
      dt2 = ((t-tprev)*(t+dt-tprev))
      dt3 = ((t-tprev)*(t-tprev2))
      dt4 = ((tprev-tprev2)*(t-tprev2))
      @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
        @inbounds DD31 = (u[i] - uprev[i])/dt1 + (uprev[i]-uprev2[i])/dt2
        @inbounds DD30 = (uprev[i] - uprev2[i])/dt3 + (uprev2[i]-uprev3[i])/dt4
        dEst = (dt^3)*abs(((DD31 - DD30)/(t+dt-tprev2))/12)
        @inbounds k[i] = dEst/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      end
      integrator.EEst = integrator.opts.internalnorm(k)
      if integrator.EEst <= 1
        copy!(cache.uprev3,uprev2)
        cache.tprev2 = tprev
      end
    elseif integrator.iter > 1
      integrator.EEst = 1
      copy!(cache.uprev3,integrator.uprev2)
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  f(t+dt,u,integrator.fsallast)
  @pack integrator = t,dt,u
end
