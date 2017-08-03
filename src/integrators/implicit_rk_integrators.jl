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
  do_newton = true
  κ = 1e-2
  tol = 1e-8

  iter += 1
  b = -z + dt*f(t+dt,uprev + z)
  dz = W\b
  ndz = abs(dz)

  while do_newton
    iter += 1
    b = -z + dt*f(t+dt,uprev + z)
    dz = W\b
    ndzprev = ndz
    ndz = abs(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    if η*ndz <= κ*tol
      do_newton=false
    end
    z = z + dz
  end

  u = uprev + z
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

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
    if has_jac(f)
      f(Val{:jac},t,uprev,J)
    else
      if alg_autodiff(integrator.alg)
        ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
      else
        Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
      end
    end
    for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]-dt*J[i,j]
    end
  end

  @. z = u - uprev
  iter = 0
  do_newton = true
  κ = 1e-2
  tol = 1e-8

  iter += 1
  f(t+dt,u,k)
  scale!(k,dt)
  k .-= z
  if has_invW(f)
    A_mul_B!(dz,W,k) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),true)
  end
  ndz = integrator.opts.internalnorm(dz)

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
    if η*ndz <= κ*tol
      do_newton=false
    end
    z .+= dz
  end

  @. u = uprev + z
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
  do_newton = true
  κ = 1e-2
  tol = 1e-8

  iter += 1
  b = -z + dto2*f(t+dt,uprev + z)
  dz = W\b
  ndz = abs(dz)

  while do_newton
    iter += 1
    b = -z + dto2*f(t+dto2,uprev + z)
    dz = W\b
    ndzprev = ndz
    ndz = abs(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    if η*ndz <= κ*tol
      do_newton=false
    end
    z = z + dz
  end

  u = uprev + 2z
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
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

  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if has_jac(f)
      f(Val{:jac},t,uprev,J)
    else
      if alg_autodiff(integrator.alg)
        ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
      else
        Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
      end
    end
    for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]-dt*J[i,j]
    end
  end

  @. z = u - uprev
  iter = 0
  do_newton = true
  κ = 1e-2
  tol = 1e-8

  iter += 1
  f(t+dt,u,k)
  scale!(k,dt)
  k .-= z
  if has_invW(f)
    A_mul_B!(dz,W,k) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),true)
  end
  ndz = integrator.opts.internalnorm(dz)

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
    if η*ndz <= κ*tol
      do_newton=false
    end
    z .+= dz
  end

  @. u = uprev + z
  f(t+dt,u,integrator.fsallast)
  @pack integrator = t,dt,u
end
