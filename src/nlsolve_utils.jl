# return `fail_convergence`
function diffeq_nlsolve!(integrator,
                         cache::OrdinaryDiffEqConstantCache,
                         # `z` is the initial guess
                         z, τ,
                         ::Type{Val{:newton}})
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = unwrap_alg(integrator, true)
  # precalculations
  κtol = κ*tol

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - τ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - τ*J
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + τ
  u = uprev + z
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = uprev + z
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return (u,true)
  end

  u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter
  return (u,false)
end

function diffeq_nlsolve!(integrator,
                         cache::OrdinaryDiffEqMutableCache,
                         τ,
                         ::Type{Val{:newton}}, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = unwrap_alg(integrator, true)
  # precalculations
  κtol = κ*tol
  # calculate W
  new_W = calc_W!(integrator, cache, τ, repeat_step)
  # initial step of Newton iteration
  iter = 1
  tstep = t + τ
  @. u = uprev + z
  f(k, u, p, tstep)
  if mass_matrix == I
    @. b = dt*k - z
  else
    A_mul_B!(vec(b),mass_matrix,vec(z))
    @. b = dt*k - b
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = uprev + z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. b = dt*k - z
    else
      A_mul_B!(vec(b),mass_matrix,vec(z))
      @. b = dt*k - b
    end
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return true
  end

  @. u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter
  return false
end
