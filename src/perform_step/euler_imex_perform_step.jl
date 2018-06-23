function initialize!(integrator, cache::IMEXEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::IMEXEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  f1 = integrator.f.f1
  du₁ = f1(uprev,p,t)
  du₂ = integrator.fsalfirst - du₁
  # Explicit Part
  tmp = uprev + dt*du₂
  # Implicit part
  # Precalculations
  γ = 1//1
  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  alg = unwrap_alg(integrator, true)
  if alg.extrapolant == :linear
    z = dt*du₁
  else # :constant
    z = zero(u)
  end

  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, 1, 1, Val{:newton})
  fail_convergence && return
  u = tmp + z

  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
