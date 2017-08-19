function initialize!(integrator,cache::ImplicitEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::ImplicitEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  uf.t = t

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    u = current_extrapolant(t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    u = uprev + integrator.fsalfirst*dt
  else # :constant
    u = uprev
  end

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end

  z = @. (u - uprev)
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  b = -z .+ dt.*f(t+dt,uprev + z)
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    b = -z .+ dt.*f(t+dt,uprev + z)
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev + z

  if integrator.opts.adaptive && integrator.success_iter > 0
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    DD3 = ((u - uprev)/((dt)*(t+dt-tprev)) + (uprev-uprev2)/((t-tprev)*(t+dt-tprev)))
    dEst = (dt^2)*abs(DD3/6)
    integrator.EEst = @. dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
  else
    integrator.EEst = 1
  end

  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end#

function initialize!(integrator,cache::ImplicitEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::ImplicitEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config = cache
  mass_matrix = integrator.sol.prob.mass_matrix


  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + integrator.fsalfirst*dt
  else # :constant
    copy!(u,uprev)
  end

  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = (u - uprev)
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  f(t+dt,u,k)
  scale!(k,dt)
  if mass_matrix == I
    k .-= z
  else
    A_mul_B!(du1,mass_matrix,z)
    k .-= du1
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz
  @. u = uprev + z

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    f(t+dt,u,k)
    scale!(k,dt)
    if mass_matrix == I
      k .-= z
    else
      A_mul_B!(du1,mass_matrix,z)
      k .-= du1
    end
    if has_invW(f)
      A_mul_B!(dz,W,k) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z .+= dz
    @. u = uprev + z
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive && integrator.success_iter > 0
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
end

function initialize!(integrator,cache::ImplicitMidpointConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::ImplicitMidpointConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  uf.t = t
  dto2 = dt/2

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    u = current_extrapolant(t+dto2,integrator)
  elseif integrator.alg.extrapolant == :linear
    u = uprev + integrator.fsalfirst*dto2
  else # :constant
    u = uprev
  end

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dto2*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dto2*J
  end
  z = @. (u - uprev)
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  b = -z .+ dt.*f(t+dto2,(uprev + u)/2)
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz
  u = uprev + z

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    b = -z .+ dt.*f(t+dto2,(uprev + u)/2)
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
    u = uprev + z
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev + z
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast

  #=
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
      integrator.EEst = @. dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      if integrator.EEst <= 1
        cache.uprev3 = uprev2
        cache.tprev2 = tprev
      end
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      cache.uprev3 = integrator.uprev2
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end
  =#

  integrator.u = u
end

function initialize!(integrator,cache::ImplicitMidpointCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::ImplicitMidpointCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  dto2 = dt/2

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dto2,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + integrator.fsalfirst*dto2
  else
    copy!(u,uprev)
  end


  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dto2*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = (u - uprev)
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  f(t+dto2,u,k)
  scale!(k,dt)
  if mass_matrix == I
    k .-= z
  else
    A_mul_B!(du1,mass_matrix,z)
    k .-= du1
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz
  # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
  @. u = uprev + z/2

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    f(t+dto2,u,k)
    scale!(k,dt)
    if mass_matrix == I
      k .-= z
    else
      A_mul_B!(du1,mass_matrix,z)
      k .-= du1
    end
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z .+= dz
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    @. u = uprev + z/2
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  @. u = uprev + z

  #=
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
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      copy!(cache.uprev3,integrator.uprev2)
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end
  =#

  f(t+dt,u,integrator.fsallast)
end

function initialize!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  uf.t = t
  dto2 = dt/2

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    u = current_extrapolant(t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    u = uprev + integrator.fsalfirst*dt
  else # :constant
    u = uprev
  end

  zprev = integrator.fsalfirst*dt

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dto2*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dto2*J
  end
  z = @. (u - uprev)
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  u = uprev + (zprev + z)/2
  b = -z .+ dt.*f(t+dt,u)
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = uprev + (zprev + z)/2
    b = -z .+ dt.*f(t+dt,u)
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev + (zprev + z)/2
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
      integrator.EEst = @. dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      if integrator.EEst <= 1
        cache.uprev3 = uprev2
        cache.tprev2 = tprev
      end
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      cache.uprev3 = integrator.uprev2
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  integrator.u = u
end

function initialize!(integrator,cache::TrapezoidCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::TrapezoidCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config = cache
  mass_matrix = integrator.sol.prob.mass_matrix
  fsalfirst = integrator.fsalfirst

  dto2 = dt/2

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + integrator.fsalfirst*dt
  else
    copy!(u,uprev)
  end


  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = (u - uprev)
  iter = 0
  κ = cache.κ
  tol = cache.tol

  iter += 1
  f(t+dt,u,k)
  scale!(k,dt)
  if mass_matrix == I
    k .-= z
  else
    A_mul_B!(du1,mass_matrix,z)
    k .-= du1
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz
  @. u = uprev + (dt*fsalfirst + z)/2

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    f(t+dt,u,k)
    scale!(k,dt)
    if mass_matrix == I
      k .-= z
    else
      A_mul_B!(du1,mass_matrix,z)
      k .-= du1
    end
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z .+= dz
    @. u = uprev + (dt*fsalfirst + z)/2
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter

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
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      copy!(cache.uprev3,integrator.uprev2)
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  f(t+dt,u,integrator.fsallast)
end

function initialize!(integrator,cache::TRBDF2ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::TRBDF2ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  uf.t = t
  γ = 2 - sqrt(2)
  ω = sqrt(2)/4
  d = γ/2

  b1 = ω
  bhat1 = (1-ω)/3
  b2 = ω
  bhat2 = (3ω + 1)/3
  b3 = d
  bhat3 = d/3

  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - d*dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - d*dt*J
  end

  #FSAL
  zprev = dt*integrator.fsalfirst

  ##### Solve Trapezoid Step

  # TODO: Add extrapolation
  zᵧ = zprev
  iter = 1
  uᵧ = @. (uprev + d*zprev) + d*zᵧ
  b = dt.*f(t+γdt,uᵧ) .- zᵧ
  Δzᵧ = W\b
  ndz = integrator.opts.internalnorm(Δzᵧ)
  zᵧ = zᵧ + Δzᵧ

  uᵧ = @. (uprev + d*zprev) + d*zᵧ

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    uᵧ = @. (uprev + d*zprev) + d*zᵧ
    b = dt.*f(t+γdt,uᵧ) .- zᵧ
    Δzᵧ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(Δzᵧ)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    zᵧ = zᵧ + Δzᵧ
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  uᵧ = @. (uprev + d*zprev) + d*zᵧ

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine
  z = (1.5 + sqrt(2))*zprev + (2.5 + 2sqrt(2))*zᵧ -
      (6 + 4.5sqrt(2))*(uᵧ - uprev)

  iter = 1
  u = @. (uprev + ω*zprev + ω*zᵧ) + d*z
  b = dt.*f(t+dt,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. (uprev + ω*zprev + ω*zᵧ) + d*z
    b = dt.*f(t+dt,u) .- z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. (uprev + ω*zprev + ω*zᵧ) + d*z

  ################################### Finalize

  integrator.fsallast = z/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = (bhat1-b1)*zprev + (bhat2-b2)*zᵧ + (bhat3-b3)*z
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end

  integrator.u = u
end

function initialize!(integrator,cache::TRBDF2Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::TRBDF2Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,uᵧ,Δzᵧ,Δz,zprev,zᵧ,z,k,J,W,jac_config,tmp = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  γ = 2 - sqrt(2)
  ω = sqrt(2)/4
  d = γ/2

  b1 = ω
  bhat1 = (1-ω)/3
  b2 = ω
  bhat2 = (3ω + 1)/3
  b3 = d
  bhat3 = d/3

  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,dt*d,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
      ddt = d*dt
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-ddt*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. zprev = dt*integrator.fsalfirst

  # TODO: Add extrapolation
  @. zᵧ = zprev
  iter = 1

  @. uᵧ = (uprev + d*zprev) + d*zᵧ
  f(t+γdt,uᵧ,k)
  @. k = dt*k - zᵧ
  if has_invW(f)
    A_mul_B!(vec(Δzᵧ),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(Δzᵧ),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(Δzᵧ)
  @. zᵧ = zᵧ + Δzᵧ

  @. uᵧ = (uprev + d*zprev) + d*zᵧ

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. uᵧ = (uprev + d*zprev) + d*zᵧ
    f(t+γdt,uᵧ,k)
    @. k = dt*k - zᵧ
    if has_invW(f)
      A_mul_B!(vec(Δzᵧ),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(Δzᵧ),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(Δzᵧ)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. zᵧ = zᵧ + Δzᵧ
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. uᵧ = (uprev + d*zprev) + d*zᵧ

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine
  a1 = (1.5 + sqrt(2)); a2 = (2.5 + 2sqrt(2)); a3 = (6 + 4.5sqrt(2))
  @. z = a1*zprev + a2*zᵧ - a3*(uᵧ - uprev)
  iter = 1

  @. u = (uprev + ω*zprev + ω*zᵧ) + d*z
  f(t+dt,u,k)
  @. k = dt*k - z
  if has_invW(f)
    A_mul_B!(vec(Δz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(Δz),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(Δz)
  @. z = z + Δz

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    @. u = (uprev + ω*zprev + ω*zᵧ) + d*z
    f(t+dt,u,k)
    @. k = dt*k - z
    if has_invW(f)
      A_mul_B!(vec(Δz),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(Δz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(Δz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z = z + Δz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = (uprev + ω*zprev + ω*zᵧ) + d*z

  ################################### Finalize

  if integrator.opts.adaptive
    be1 = (bhat1-b1); be2 = (bhat2-b2); be3 = (bhat3-b3)
    @. tmp = be1*zprev + be2*zᵧ + be3*z
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::SDIRK2ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SDIRK2ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end

  # TODO: Add extrapolant initial guess
  u = uprev

  ##### Step 1

  z₁ = u - uprev
  iter = 1
  u = @. uprev + z₁
  b = dt.*f(t+dt,u) .- z₁
  dz₁ = W\b
  ndz = integrator.opts.internalnorm(dz₁)
  z₁ = z₁ + dz₁

  u = @. uprev + z₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + z₁
    b = dt.*f(t+dt,u) .- z₁
    dz₁ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ, c₂ = 0 => z₂ = α₁z₁ = 0
  z₂ = zero(u)

  iter = 1
  u = @. uprev - z₁ + z₂
  b = dt.*f(t,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev - z₁ + z₂
    b = dt.*f(t,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    z₂ = z₂ + dz₂
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + z₁/2 + z₂/2

  integrator.fsallast = f(t,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. z₁/2 - z₂/2
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end

  integrator.u = u
end

function initialize!(integrator,cache::SDIRK2Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::SDIRK2Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,z₁,z₂,k,J,W,jac_config,tmp = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+γ*dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    u .= uprev .+ integrator.fsalfirst.*γ.*dt
  else
    copy!(u,uprev)
  end

  ##### Step 1

  @. z₁ = u - uprev

  iter = 1
  @. u = uprev + z₁
  f(t+dt,u,k)
  @. k = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₁),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₁)
  @. z₁ = z₁ + dz₁
  @. u = uprev + z₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1

    @. u = uprev + z₁
    f(t+dt,u,k)
    @. k = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₁),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess
  @. z₂ = zero(u)

  iter = 1
  @.  u = uprev - z₁ + z₂
  f(t,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter


    @.  u = uprev - z₁ + z₂
    f(t,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  @. u = uprev + z₁/2 + z₂/2

  if integrator.opts.adaptive
    @. tmp = z₁/2 - z₂/2
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  f(t,u,integrator.fsallast)
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::SSPSDIRK2ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SSPSDIRK2ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  uf.t = t

  γ = eltype(u)(1//4)
  c2 = typeof(t)(3//4)

  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
  end

  # TODO: Add extrapolant initial guess
  u = uprev

  ##### Step 1

  z₁ = u - uprev
  iter = 1
  u = @. uprev + γ*z₁
  b = dt.*f(t+dt,u) .- z₁
  dz₁ = W\b
  ndz = integrator.opts.internalnorm(dz₁)
  z₁ = z₁ + dz₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁
    b = dt.*f(t+dt,u) .- z₁
    dz₁ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ
  z₂ = c2/γ

  iter = 1
  u = @. uprev + z₁/2 + z₂/4
  b = dt.*f(t,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + z₁/2 + z₂/4
    b = dt.*f(t,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    z₂ = z₂ + dz₂
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + z₁/2 + z₂/2

  integrator.fsallast = f(t,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  integrator.u = u
end

function initialize!(integrator,cache::SSPSDIRK2Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::SSPSDIRK2Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,z₁,z₂,k,J,W,jac_config,tmp = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  γ = eltype(u)(1//4)
  c2 = typeof(t)(3//4)
  γdt = γ*dt
  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+γ*dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    u .= uprev .+ integrator.fsalfirst.*γ.*dt
  else
    copy!(u,uprev)
  end

  ##### Step 1

  @. z₁ = u - uprev

  iter = 1
  @. u = uprev + γ*z₁
  f(t+dt,u,k)
  @. k = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₁),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₁)
  @. z₁ = z₁ + dz₁
  @. u = uprev + z₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(t+dt,u,k)
    @. k = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₁),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ
  @. z₂ = c2/γ

  iter = 1
  @.  u = uprev + z₁/2 + z₂/4
  f(t,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    @.  u = uprev + z₁/2 + z₂/4
    f(t,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  @. u = uprev + z₁/2 + z₂/2

  f(t,u,integrator.fsallast)
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::Union{Kvaerno3ConstantCache,KenCarp3ConstantCache},f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Union{Kvaerno3ConstantCache,KenCarp3ConstantCache},f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a31,a32,a41,a42,a43,bhat1,bhat2,bhat3,bhat4,c3,α31,α32 = cache.tab
  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  #FSAL Step 1
  z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = z₁

  iter = 1
  u = @. uprev + γ*(z₁+z₂)
  b = dt.*f(t+2γ*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*(z₁+z₂)
    b = dt.*f(t+2γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess is from Hermite derivative on z₁ and z₂
  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  if typeof(cache) <: Kvaerno3ConstantCache
    z₄ = @. a31*z₁ + a32*z₂ + γ*z₃ # use yhat as prediction
  elseif typeof(cache) <: KenCarp3ConstantCache
    @unpack α41,α42 = cache.tab
    z₄ = @. α41*z₁ + α42*z₂
  end

  iter = 1
  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  b = dt.*f(t+dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    b = dt.*f(t+dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

  integrator.fsallast = z₄/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. (bhat1-a41)*z₁ + (bhat2-a42)*z₂ + (bhat3-a43)*z₃ + (bhat4-γ)*z₄
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end

  integrator.u = u
end

function initialize!(integrator,cache::Union{Kvaerno3Cache,KenCarp3Cache},f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::Union{Kvaerno3Cache,KenCarp3Cache},f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,z₁,z₂,z₃,z₄,k,J,W,jac_config,tmp = cache
  @unpack γ,a31,a32,a41,a42,a43,bhat1,bhat2,bhat3,bhat4,c3,α31,α32 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  #FSAL Step 1
  @. z₁ = dt*integrator.fsalfirst

  ##### Step 1

  # TODO: Add extrapolation for guess
  @. z₂ = z₁

  iter = 1
  @. u = uprev + γ*(z₁+z₂)
  f(t+2γ*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*(z₁+z₂)
    f(t+2γ*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess is from Hermite derivative on z₁ and z₂
  @. z₃ = α31*z₁ + α32*z₂

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  if typeof(cache) <: Kvaerno3Cache
    @. z₄ = a31*z₁ + a32*z₂ + γ*z₃ # use yhat as prediction
  elseif typeof(cache) <: KenCarp3Cache
    @unpack α41,α42 = cache.tab
    @. z₄ = α41*z₁ + α42*z₂
  end
  iter = 1

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  f(t+dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    f(t+dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

  if integrator.opts.adaptive
    @. tmp = (bhat1-a41)*z₁ + (bhat2-a42)*z₂ + (bhat3-a43)*z₃ + (bhat4-γ)*z₄
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₄/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::Cash4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Cash4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack b1hat1, b2hat1, b3hat1, b4hat1, b1hat2, b2hat2, b3hat2, b4hat2 = cache.tab
  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  # TODO: Add extrapolation for guess
  z₁ = zero(u)

  iter = 1
  u = @. uprev + γ*z₁
  b = dt.*f(t+γ*dt,u) .- z₁
  dz₁ = W\b
  ndz = integrator.opts.internalnorm(dz₁)
  z₁ = z₁ + dz₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁
    b = dt.*f(t+γ*dt,u) .- z₁
    dz₁ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = zero(u)

  iter = 1
  u = @. uprev + a21*z₁ + γ*z₂
  b = dt.*f(t+c2*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a21*z₁ + γ*z₂
    b = dt.*f(t+c2*γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess starts from z₁
  z₃ = @. z₁

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  z₄ = @. z₃

  iter = 1
  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  b = dt.*f(t+c4*dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    b = dt.*f(t+c4*dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = @. b1hat2*z₁ + b2hat2*z₂ + b3hat2*z₃ + b4hat2*z₄

  iter = 1
  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  b = dt.*f(t+dt,u) .- z₅
  dz₅ = W\b
  ndz = integrator.opts.internalnorm(dz₅)
  z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₅ + γ*z₅
    b = dt.*f(t+dt,u) .- z₅
    dz₅ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅

  integrator.fsallast = z₅/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if integrator.alg.embedding == 3
      tmp = @. (b1hat2-a51)*z₁ + (b2hat2-a52)*z₂ + (b3hat2-a53)*z₃ + (b4hat2-a54)*z₄ - γ*z₅
    else
      tmp = @. (b1hat1-a51)*z₁ + (b2hat1-a52)*z₂ + (b3hat1-a53)*z₃ + (b4hat1-a54)*z₄ - γ*z₅
    end
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.u = u
end

function initialize!(integrator,cache::Cash4Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::Cash4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,dz₅,z₁,z₂,z₃,z₄,z₅,k,J,W,jac_config,tmp = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack b1hat1, b2hat1, b3hat1, b4hat1, b1hat2, b2hat2, b3hat2, b4hat2 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  # TODO: Add extrapolation for guess
  @. z₁ .= zero(z₁)

  iter = 1
  @. u = uprev + γ*z₁
  f(t+γ*dt,u,k)
  @. k = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₁),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₁)
  @. z₁ = z₁ + dz₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(t+γ*dt,u,k)
    @. k = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₁),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  # TODO: Add extrapolation for guess
  @. z₂ .= zero(z₂)

  iter = 1
  @. u = uprev + a21*z₁ + γ*z₂
  f(t+c2*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a21*z₁ + γ*z₂
    f(t+c2*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess starts from z₁
  @. z₃ = z₁

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = z₃

  iter = 1

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  f(t+dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    f(t+dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use constant z prediction
  @. z₅ = b1hat2*z₁ + b2hat2*z₂ + b3hat2*z₃ + b4hat2*z₄

  iter = 1
  @. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  f(t+dt,u,k)
  @. k = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₅)
  @. z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
    f(t+dt,u,k)
    @. k = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  @. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅

  if integrator.opts.adaptive
    if integrator.alg.embedding == 3
      @. tmp = (b1hat2-a51)*z₁ + (b2hat2-a52)*z₂ + (b3hat2-a53)*z₃ + (b4hat2-a54)*z₄ - γ*z₅
    else
      @. tmp = (b1hat1-a51)*z₁ + (b2hat1-a52)*z₂ + (b3hat1-a53)*z₃ + (b4hat1-a54)*z₄ - γ*z₅
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₅/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::Hairer4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Hairer4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α43 = cache.tab
  @unpack bhat1, bhat2, bhat3, bhat4 = cache.tab
  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  # TODO: Add extrapolation for guess
  z₁ = zero(u)

  iter = 1
  u = @. uprev + γ*z₁
  b = dt.*f(t+γ*dt,u) .- z₁
  dz₁ = W\b
  ndz = integrator.opts.internalnorm(dz₁)
  z₁ = z₁ + dz₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁
    b = dt.*f(t+γ*dt,u) .- z₁
    dz₁ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  z₂ = α21*z₁

  iter = 1
  u = @. uprev + a21*z₁ + γ*z₂
  b = dt.*f(t+c2*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a21*z₁ + γ*z₂
    b = dt.*f(t+c2*γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α43*z₃

  iter = 1
  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  b = dt.*f(t+c4*dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    b = dt.*f(t+c4*dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = @. bhat1*z₁ + bhat2*z₂ + bhat3*z₃ + bhat4*z₄

  iter = 1
  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  b = dt.*f(t+dt,u) .- z₅
  dz₅ = W\b
  ndz = integrator.opts.internalnorm(dz₅)
  z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₅ + γ*z₅
    b = dt.*f(t+dt,u) .- z₅
    dz₅ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅

  integrator.fsallast = z₅/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. (bhat1-a51)*z₁ + (bhat2-a52)*z₂ + (bhat3-a53)*z₃ + (bhat4-a54)*z₄ - γ*z₅
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.u = u
end

function initialize!(integrator,cache::Hairer4Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::Hairer4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,dz₅,z₁,z₂,z₃,z₄,z₅,k,J,W,jac_config,tmp = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α43 = cache.tab
  @unpack bhat1, bhat2, bhat3, bhat4 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+γ*dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + γ*dt*integrator.fsalfirst
  else
    copy!(u,uprev)
  end

  @. z₁ = u - uprev

  iter = 1
  @. u = uprev + γ*z₁
  f(t+γ*dt,u,k)
  @. k = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₁),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₁)
  @. z₁ = z₁ + dz₁

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(t+γ*dt,u,k)
    @. k = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz₁),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₁),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₁)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₁ = z₁ + dz₁
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  @. z₂ = α21*z₁

  iter = 1
  @. u = uprev + a21*z₁ + γ*z₂
  f(t+c2*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a21*z₁ + γ*z₂
    f(t+c2*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  @. z₃ = α31*z₁ + α32*z₂

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α43*z₃

  iter = 1

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  f(t+c4*dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    f(t+c4*dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat prediction
  @. z₅ = bhat1*z₁ + bhat2*z₂ + bhat3*z₃ + bhat4*z₄

  iter = 1
  #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end
  f(t+dt,u,k)
  @. k = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₅)
  @. z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
    end
    f(t+dt,u,k)
    @. k = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end

  if integrator.opts.adaptive
    #@. tmp = (bhat1-a51)*z₁ + (bhat2-a52)*z₂ + (bhat3-a53)*z₃ + (bhat4-a54)*z₄ - γ*z₅
    @tight_loop_macros for i in eachindex(u)
      tmp[i] = (bhat1-a51)*z₁[i] + (bhat2-a52)*z₂[i] + (bhat3-a53)*z₃[i] + (bhat4-a54)*z₄[i] - γ*z₅[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₅/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::Kvaerno4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end


@muladd function perform_step!(integrator,cache::Kvaerno4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α42 = cache.tab
  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  iter = 1
  u = @. uprev + γ*z₁ + γ*z₂
  b = dt.*f(t+2*γ*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁ + γ*z₂
    b = dt.*f(t+2*γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂

  iter = 1
  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  b = dt.*f(t+c4*dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    b = dt.*f(t+c4*dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = @. a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

  iter = 1
  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  b = dt.*f(t+dt,u) .- z₅
  dz₅ = W\b
  ndz = integrator.opts.internalnorm(dz₅)
  z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₅ + γ*z₅
    b = dt.*f(t+dt,u) .- z₅
    dz₅ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅

  integrator.fsallast = z₅/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. (a41-a51)*z₁ + (a42-a52)*z₂ + (a43-a53)*z₃ + (γ-a54)*z₄ - γ*z₅
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.u = u
end

function initialize!(integrator,cache::Kvaerno4Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::Kvaerno4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,dz₅,z₁,z₂,z₃,z₄,z₅,k,J,W,jac_config,tmp = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α42 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  @. z₂ = zero(u)

  iter = 1
  @. u = uprev + γ*z₁ + γ*z₂
  f(t+2γ*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁ + γ*z₂
    f(t+2γ*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  @. z₃ = α31*z₁ + α32*z₂

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂

  iter = 1

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  f(t+c4*dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    f(t+c4*dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat prediction
  @. z₅ = a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

  iter = 1
  #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end
  f(t+dt,u,k)
  @. k = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₅)
  @. z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
    end
    f(t+dt,u,k)
    @. k = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end

  if integrator.opts.adaptive
    @. tmp = (a41-a51)*z₁ + (a42-a52)*z₂ + (a43-a53)*z₃ + (γ-a54)*z₄ - γ*z₅
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₅/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::KenCarp4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::KenCarp4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c3,c4,c5 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α53,α54,α61,α62,α63,α64,α65 = cache.tab
  @unpack bhat1,bhat3,bhat4,bhat5,bhat6 = cache.tab

  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  iter = 1
  u = @. uprev + γ*z₁ + γ*z₂
  b = dt.*f(t+2*γ*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁ + γ*z₂
    b = dt.*f(t+2*γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂

  iter = 1
  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  b = dt.*f(t+c4*dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    b = dt.*f(t+c4*dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₅ = @. α51*z₁ + α52*z₂ + α53*z₃ + α54*z₄

  iter = 1
  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  b = dt.*f(t+c5*dt,u) .- z₅
  dz₅ = W\b
  ndz = integrator.opts.internalnorm(dz₅)
  z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₅ + γ*z₅
    b = dt.*f(t+c5*dt,u) .- z₅
    dz₅ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₆ = @. α61*z₁ + α62*z₂ + α63*z₃ + α64*z₄

  iter = 1
  u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  b = dt.*f(t+dt,u) .- z₆
  dz₆ = W\b
  ndz = integrator.opts.internalnorm(dz₆)
  z₆ = z₆ + dz₆

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
    b = dt.*f(t+dt,u) .- z₆
    dz₆ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₆)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₆ = z₆ + dz₆
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆

  integrator.fsallast = z₆/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. (bhat1-a61)*z₁ + (bhat3-a63)*z₃ + (bhat4-a64)*z₄ + (bhat5-a65)*z₅ + (bhat6-γ)*z₆
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.u = u
end

function initialize!(integrator,cache::KenCarp4Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::KenCarp4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,dz₅,dz₆,z₁,z₂,z₃,z₄,z₅,z₆,k,J,W,jac_config,tmp = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c3,c4,c5 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α53,α54,α61,α62,α63,α64,α65 = cache.tab
  @unpack bhat1,bhat3,bhat4,bhat5,bhat6 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  @. z₂ = zero(u)

  iter = 1
  @. u = uprev + γ*z₁ + γ*z₂
  f(t+2γ*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁ + γ*z₂
    f(t+2γ*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  @. z₃ = α31*z₁ + α32*z₂

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂

  iter = 1

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  f(t+c4*dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    f(t+c4*dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  @. z₅ = α51*z₁ + α52*z₂ + α53*z₃ + α54*z₄

  iter = 1
  #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end
  f(t+c5*dt,u,k)
  @. k = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₅)
  @. z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
    end
    f(t+c5*dt,u,k)
    @. k = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  #@. z₆ = α61*z₁ + α62*z₂ + α63*z₃ + α64*z₄ + α65*z₅
  @tight_loop_macros for i in eachindex(u)
    z₆[i] = α61*z₁[i] + α62*z₂[i] + α63*z₃[i] + α64*z₄[i] + α65*z₅[i]
  end

  iter = 1
  #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
  end
  f(t+dt,u,k)
  @. k = dt*k - z₆
  if has_invW(f)
    A_mul_B!(vec(dz₆),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₆),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₆)
  @. z₆ = z₆ + dz₆

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
    end
    f(t+dt,u,k)
    @. k = dt*k - z₆
    if has_invW(f)
      A_mul_B!(vec(dz₆),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₆),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₆)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₆ = z₆ + dz₆
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
  end

  if integrator.opts.adaptive
    #@. tmp = (bhat1-a61)*z₁ + (bhat3-a63)*z₃ + (bhat4-a64)*z₄ + (bhat5-a65)*z₅ + (bhat6-γ)*z₆
    for i in eachindex(u)
      tmp[i] = (bhat1-a61)*z₁[i] + (bhat3-a63)*z₃[i] + (bhat4-a64)*z₄[i] + (bhat5-a65)*z₅[i] + (bhat6-γ)*z₆[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₆/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::Kvaerno5ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Kvaerno5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,c3,c4,c5,c6 = cache.tab
  @unpack α31,α32,α41,α42,α43,α51,α52,α53,α61,α62,α63 = cache.tab

  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  iter = 1
  u = @. uprev + γ*z₁ + γ*z₂
  b = dt.*f(t+2*γ*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁ + γ*z₂
    b = dt.*f(t+2*γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂ + α43*z₃

  iter = 1
  u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  b = dt.*f(t+c4*dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    b = dt.*f(t+c4*dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₅ = @. α51*z₁ + α52*z₂ + α53*z₃

  iter = 1
  u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  b = dt.*f(t+c5*dt,u) .- z₅
  dz₅ = W\b
  ndz = integrator.opts.internalnorm(dz₅)
  z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₅ + γ*z₅
    b = dt.*f(t+c5*dt,u) .- z₅
    dz₅ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  z₆ = @. α61*z₁ + α62*z₂ + α63*z₃

  iter = 1
  u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  b = dt.*f(t+c6*dt,u) .- z₆
  dz₆ = W\b
  ndz = integrator.opts.internalnorm(dz₆)
  z₆ = z₆ + dz₆

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
    b = dt.*f(t+c6*dt,u) .- z₆
    dz₆ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₆)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₆ = z₆ + dz₆
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  # Prediction from embedding
  z₇ = @. a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆

  iter = 1
  u = @. uprev + a71*z₁ +  a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
  b = dt.*f(t+dt,u) .- z₇
  dz₇ = W\b
  ndz = integrator.opts.internalnorm(dz₇)
  z₇ = z₇ + dz₇

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a71*z₁ +  a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
    b = dt.*f(t+dt,u) .- z₇
    dz₇ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₇)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₇ = z₇ + dz₇
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a71*z₁ +  a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇

  integrator.fsallast = z₇/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. (a61-a71)*z₁ + (a63-a73)*z₃ + (a64-a74)*z₄ + (a65-a75)*z₅ + (γ-a76)*z₆ - γ*z₇
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.u = u
end

function initialize!(integrator,cache::Kvaerno5Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::Kvaerno5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,dz₅,dz₆,dz₇,z₁,z₂,z₃,z₄,z₅,z₆,z₇,k,J,W,jac_config,tmp = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,c3,c4,c5,c6 = cache.tab
  @unpack α31,α32,α41,α42,α43,α51,α52,α53,α61,α62,α63 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  @. z₂ = zero(u)

  iter = 1
  @. u = uprev + γ*z₁ + γ*z₂
  f(t+2γ*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁ + γ*z₂
    f(t+2γ*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  @. z₃ = α31*z₁ + α32*z₂

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂ + α43*z₃

  iter = 1

  @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
  f(t+c4*dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
    f(t+c4*dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  @. z₅ = α51*z₁ + α52*z₂ + α53*z₃

  iter = 1
  #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end
  f(t+c5*dt,u,k)
  @. k = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₅)
  @. z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a51*z₁[i] + a52*z₂[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
    end
    f(t+c5*dt,u,k)
    @. k = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  #@. z₆ = α61*z₁ + α62*z₂ + α63*z₃
  @tight_loop_macros for i in eachindex(u)
    z₆[i] = α61*z₁[i] + α62*z₂[i] + α63*z₃[i]
  end

  iter = 1
  #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
  end
  f(t+c6*dt,u,k)
  @. k = dt*k - z₆
  if has_invW(f)
    A_mul_B!(vec(dz₆),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₆),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₆)
  @. z₆ = z₆ + dz₆

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
    end
    f(t+c6*dt,u,k)
    @. k = dt*k - z₆
    if has_invW(f)
      A_mul_B!(vec(dz₆),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₆),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₆)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₆ = z₆ + dz₆
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  # Prediction is embedded method

  #@. z₇ = a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  @tight_loop_macros for i in eachindex(u)
    z₇[i] = a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
  end

  iter = 1
  #@. u = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i] + γ*z₇[i]
  end
  f(t+dt,u,k)
  @. k = dt*k - z₇
  if has_invW(f)
    A_mul_B!(vec(dz₇),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₇),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₇)
  @. z₇ = z₇ + dz₇

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i] + γ*z₇[i]
    end
    f(t+dt,u,k)
    @. k = dt*k - z₇
    if has_invW(f)
      A_mul_B!(vec(dz₇),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₇),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₇)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₇ = z₇ + dz₇
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  #@. u = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i] + γ*z₇[i]
  end

  if integrator.opts.adaptive
    #@. tmp = (a61-a71)*z₁ + (a63-a73)*z₃ + (a64-a74)*z₄ + (a65-a75)*z₅ + (γ-a76)*z₆ - γ*z₇
    for i in eachindex(u)
      tmp[i] = (a61-a71)*z₁[i] + (a63-a73)*z₃[i] + (a64-a74)*z₄[i] + (a65-a75)*z₅[i] + (γ-a76)*z₆[i] - γ*z₇[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₇/dt
  cache.ηold = η
  cache.newton_iters = iter
end

function initialize!(integrator,cache::KenCarp5ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::KenCarp5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  @unpack γ,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a84,a85,a86,a87,c3,c4,c5,c6,c7 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α61,α62,α71,α72,α73,α74,α75,α81,α82,α83,α84,α85 = cache.tab
  @unpack bhat1,bhat4,bhat5,bhat6,bhat7,bhat8 = cache.tab

  uf.t = t
  γdt = γ*dt
  κ = cache.κ
  tol = cache.tol

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*γ*J
  end

  z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  iter = 1
  u = @. uprev + γ*z₁ + γ*z₂
  b = dt.*f(t+2*γ*dt,u) .- z₂
  dz₂ = W\b
  ndz = integrator.opts.internalnorm(dz₂)
  z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁ + γ*z₂
    b = dt.*f(t+2*γ*dt,u) .- z₂
    dz₂ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
  b = dt.*f(t+c3*dt,u) .- z₃
  dz₃ = W\b
  ndz = integrator.opts.internalnorm(dz₃)
  z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a31*z₁ + a32*z₂ + γ*z₃
    b = dt.*f(t+c3*dt,u) .- z₃
    dz₃ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂

  iter = 1
  u = @. uprev + a41*z₁ + a43*z₃ + γ*z₄
  b = dt.*f(t+c4*dt,u) .- z₄
  dz₄ = W\b
  ndz = integrator.opts.internalnorm(dz₄)
  z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a41*z₁ + a43*z₃ + γ*z₄
    b = dt.*f(t+c4*dt,u) .- z₄
    dz₄ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₅ = @. α51*z₁ + α52*z₂

  iter = 1
  u = @. uprev + a51*z₁ + a53*z₃ + a54*z₄ + γ*z₅
  b = dt.*f(t+c5*dt,u) .- z₅
  dz₅ = W\b
  ndz = integrator.opts.internalnorm(dz₅)
  z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a51*z₁ + a53*z₃ + a54*z₅ + γ*z₅
    b = dt.*f(t+c5*dt,u) .- z₅
    dz₅ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  z₆ = @. α61*z₁ + α62*z₂

  iter = 1
  u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  b = dt.*f(t+c6*dt,u) .- z₆
  dz₆ = W\b
  ndz = integrator.opts.internalnorm(dz₆)
  z₆ = z₆ + dz₆

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
    b = dt.*f(t+c6*dt,u) .- z₆
    dz₆ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₆)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₆ = z₆ + dz₆
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  z₇ = @. α71*z₁ + α72*z₂ + α73*z₃ + α74*z₄ + α75*z₅

  iter = 1
  u = @. uprev + a71*z₁ +  a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
  b = dt.*f(t+c7*dt,u) .- z₇
  dz₇ = W\b
  ndz = integrator.opts.internalnorm(dz₇)
  z₇ = z₇ + dz₇

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a71*z₁ +  a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
    b = dt.*f(t+c7*dt,u) .- z₇
    dz₇ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₇)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₇ = z₇ + dz₇
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 8

  # Prediction from embedding
  z₈ = @. α81*z₁ + α82*z₂ + α83*z₃ + α84*z₄ + α85*z₅

  iter = 1
  u = @. uprev + a81*z₁ +  a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇ + γ*z₈
  b = dt.*f(t+dt,u) .- z₈
  dz₈ = W\b
  ndz = integrator.opts.internalnorm(dz₈)
  z₈ = z₈ + dz₈

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + a81*z₁ +  a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇ + γ*z₈
    b = dt.*f(t+dt,u) .- z₈
    dz₈ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₈)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z₈ = z₈ + dz₈
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  u = @. uprev + a81*z₁ +  a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇ + γ*z₈

  integrator.fsallast = z₈/dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. (bhat1-a81)*z₁ + (bhat4-a84)*z₄ + (bhat5-a85)*z₅ + (bhat6-a86)*z₆ + (bhat7-a87)*z₇ + (bhat8-γ)*z₈
    if integrator.alg.smooth_est # From Shampine
      Est = W\tmp
    else
      Est = tmp
    end
    integrator.EEst = integrator.opts.internalnorm(@. abs(Est)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.u = u
end

function initialize!(integrator,cache::KenCarp5Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator,cache::KenCarp5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz₁,dz₂,dz₃,dz₄,dz₅,dz₆,dz₇,dz₈,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,k,J,W,jac_config,tmp = cache
  @unpack γ,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a84,a85,a86,a87,c3,c4,c5,c6,c7 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α61,α62,α71,α72,α73,α74,α75,α81,α82,α83,α84,α85 = cache.tab
  @unpack bhat1,bhat4,bhat5,bhat6,bhat7,bhat8 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  uf.t = t
  κ = cache.κ
  tol = cache.tol

  if has_invW(f)
    f(Val{:invW},t,uprev,γ*dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
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
          @inbounds W[i,j] = mass_matrix[i,j]-dt*γ*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  @. z₂ = zero(u)

  iter = 1
  @. u = uprev + γ*z₁ + γ*z₂
  f(t+2γ*dt,u,k)
  @. k = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₂),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz₂)
  @. z₂ = z₂ + dz₂

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁ + γ*z₂
    f(t+2γ*dt,u,k)
    @. k = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz₂),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₂),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₂)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₂ = z₂ + dz₂
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  @. z₃ = α31*z₁ + α32*z₂

  iter = 1
  @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
  f(t+c3*dt,u,k)
  @. k = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₃)
  @. z₃ = z₃ + dz₃

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a31*z₁ + a32*z₂ + γ*z₃
    f(t+c3*dt,u,k)
    @. k = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz₃),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₃),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₃)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₃ = z₃ + dz₃
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂

  iter = 1

  @. u = uprev + a41*z₁ + a43*z₃ + γ*z₄
  f(t+c4*dt,u,k)
  @. k = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₄)
  @. z₄ = z₄ + dz₄

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + a41*z₁ + a43*z₃ + γ*z₄
    f(t+c4*dt,u,k)
    @. k = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz₄),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₄),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₄)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₄ = z₄ + dz₄
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  @. z₅ = α51*z₁ + α52*z₂

  iter = 1
  #@. u = uprev + a51*z₁ + a53*z₃ + a54*z₄ + γ*z₅
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a51*z₁[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
  end
  f(t+c5*dt,u,k)
  @. k = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₅)
  @. z₅ = z₅ + dz₅

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a51*z₁ + a53*z₃ + a54*z₄ + γ*z₅
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a51*z₁[i] + a53*z₃[i] + a54*z₄[i] + γ*z₅[i]
    end
    f(t+c5*dt,u,k)
    @. k = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz₅),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₅),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₅)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₅ = z₅ + dz₅
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  #@. z₆ = α61*z₁ + α62*z₂
  @tight_loop_macros for i in eachindex(u)
    z₆[i] = α61*z₁[i] + α62*z₂[i]
  end

  iter = 1
  #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
  end
  f(t+c6*dt,u,k)
  @. k = dt*k - z₆
  if has_invW(f)
    A_mul_B!(vec(dz₆),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₆),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₆)
  @. z₆ = z₆ + dz₆

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
    end
    f(t+dt,u,k)
    @. k = dt*k - z₆
    if has_invW(f)
      A_mul_B!(vec(dz₆),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₆),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₆)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₆ = z₆ + dz₆
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  #@. z₇ = α71*z₁ + α72*z₂ + α73*z₃ + α74*z₄ + α75*z₅
  @tight_loop_macros for i in eachindex(u)
    z₇[i] = α71*z₁[i] + α72*z₂[i] + α73*z₃[i] + α74*z₄[i] + α75*z₅[i]
  end

  iter = 1
  #@. u = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i] + γ*z₇[i]
  end
  f(t+c7*dt,u,k)
  @. k = dt*k - z₇
  if has_invW(f)
    A_mul_B!(vec(dz₇),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₇),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₇)
  @. z₇ = z₇ + dz₇

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆ + γ*z₇
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i] + γ*z₇[i]
    end
    f(t+c7*dt,u,k)
    @. k = dt*k - z₇
    if has_invW(f)
      A_mul_B!(vec(dz₇),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₇),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₇)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₇ = z₇ + dz₇
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 8

  #@. z₈ = α81*z₁ + α82*z₂ + α83*z₃ + α84*z₄ + α85*z₅
  @tight_loop_macros for i in eachindex(u)
    z₈[i] = α81*z₁[i] + α82*z₂[i] + α83*z₃[i] + α84*z₄[i] + α85*z₅[i]
  end

  iter = 1
  #@. u = uprev + a81*z₁ + a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇ + γ*z₈
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a81*z₁[i] + a84*z₄[i] + a85*z₅[i] + a86*z₆[i] + a87*z₇[i] + γ*z₈[i]
  end
  f(t+dt,u,k)
  @. k = dt*k - z₈
  if has_invW(f)
    A_mul_B!(vec(dz₈),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz₈),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(dz₈)
  @. z₈ = z₈ + dz₈

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κ*tol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    #@. u = uprev + a81*z₁ + a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇ + γ*z₈
    @tight_loop_macros for i in eachindex(u)
      u[i] = uprev[i] + a81*z₁[i] + a84*z₄[i] + a85*z₅[i] + a86*z₆[i] + a87*z₇[i] + γ*z₈[i]
    end
    f(t+dt,u,k)
    @. k = dt*k - z₈
    if has_invW(f)
      A_mul_B!(vec(dz₈),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz₈),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz₈)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    @. z₈ = z₈ + dz₈
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################### Finalize

  #@. u = uprev + a81*z₁ + a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇ + γ*z₈
  @tight_loop_macros for i in eachindex(u)
    u[i] = uprev[i] + a81*z₁[i] + a84*z₄[i] + a85*z₅[i] + a86*z₆[i] + a87*z₇[i] + γ*z₈[i]
  end

  if integrator.opts.adaptive
    #@. tmp = (bhat1-a81)*z₁ + (bhat4-a84)*z₄ + (bhat5-a85)*z₅ + (bhat6-a86)*z₆ + (bhat7-a87)*z₇ + (bhat8-γ)*z₈
    #@. tmp = (a61-a71)*z₁ + (a63-a73)*z₃ + (a64-a74)*z₄ + (a65-a75)*z₅ + (γ-a76)*z₆ - γ*z₇
    for i in eachindex(u)
      tmp[i] = (bhat1-a81)*z₁[i] + (bhat4-a84)*z₄[i] + (bhat5-a85)*z₅[i] + (bhat6-a86)*z₆[i] + (bhat7-a87)*z₇[i] + (bhat8-γ)*z₈[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      tmp[i] = abs(k[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z₈/dt
  cache.ηold = η
  cache.newton_iters = iter
end
