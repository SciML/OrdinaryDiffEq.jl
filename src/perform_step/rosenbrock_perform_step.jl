function initialize!(integrator, cache::Union{Rosenbrock23Cache,
                                              Rosenbrock32Cache,})
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [k₁,k₂]
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

function initialize!(integrator, cache::Union{Rosenbrock23ConstantCache,
                                              Rosenbrock32ConstantCache})
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = zero(integrator.fsalfirst)
  integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock23Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,opts = integrator
  @unpack k₁,k₂,k₃,du1,du2,f₁,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,jac_config,atmp,weight = cache
  @unpack c₃₂,d = cache.tab

  # Assignments
  sizeu  = size(u)
  mass_matrix = integrator.f.mass_matrix

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  new_W = calc_rosenbrock_differentiation!(integrator, cache, γ, γ, repeat_step, false)

  calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
                       opts.abstol, opts.reltol, opts.internalnorm, t)

  cache.linsolve(vec(k₁), W, vec(linsolve_tmp), new_W,
      Pl=DiffEqBase.ScaleVector(weight, true),
      Pr=DiffEqBase.ScaleVector(weight, false), reltol=opts.reltol)

  @.. k₁ = -k₁
  integrator.destats.nsolve += 1

  @.. u = uprev + dto2*k₁
  f(f₁,u,p,t+dto2)
  integrator.destats.nf += 1

  if mass_matrix == I
    tmp .= k₁
  else
    mul!(vec(tmp),mass_matrix,k₁)
  end

  @.. linsolve_tmp = f₁ - tmp
  cache.linsolve(vec(k₂), W, vec(linsolve_tmp))
  @.. k₂ = -k₂
  integrator.destats.nsolve += 1

  @.. k₂ += k₁
  @.. u = uprev + dt*k₂

  if integrator.opts.adaptive
    f( fsallast,  u, p, t+dt)
    integrator.destats.nf += 1

    if mass_matrix == I
      @.. linsolve_tmp = fsallast - c₃₂*(k₂-f₁) - 2(k₁-fsalfirst) + dt*dT
    else
      @.. du2 = c₃₂*k₂ + 2k₁
      mul!(du1,mass_matrix,du2)
      @.. linsolve_tmp = fsallast - du1 + c₃₂*f₁ + 2fsalfirst + dt*dT
    end


    cache.linsolve(vec(k₃), W, vec(linsolve_tmp))
    @.. k₃ = -k₃
    integrator.destats.nsolve += 1

    @.. tmp = dto6*(k₁ - 2*k₂ + k₃)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end

@muladd function perform_step!(integrator, cache::Rosenbrock32Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k₁,k₂,k₃,du1,du2,f₁,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,jac_config,atmp = cache
  @unpack c₃₂,d = cache.tab

  # Assignments
  sizeu  = size(u)
  mass_matrix = integrator.f.mass_matrix

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  calc_rosenbrock_differentiation!(integrator, cache, γ, γ, repeat_step, false)

  cache.linsolve(vec(k₁), W, vec(linsolve_tmp), !repeat_step)
  @.. k₁ = -k₁
  integrator.destats.nsolve += 1

  @.. u = uprev + dto2*k₁
  f(f₁,u,p,t+dto2)
  integrator.destats.nf += 1

  if mass_matrix == I
    tmp .= k₁
  else
    mul!(tmp,mass_matrix,k₁)
  end

  @.. linsolve_tmp = f₁ - tmp

  cache.linsolve(vec(k₂), W, vec(linsolve_tmp))
  @.. k₂ = -k₂
  integrator.destats.nsolve += 1

  @.. k₂ += k₁
  @.. tmp = uprev + dt*k₂
  f( fsallast,  tmp, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = fsallast - c₃₂*(k₂-f₁) - 2(k₁-fsalfirst) + dt*dT
  else
    @.. du2 = c₃₂*k₂ + 2k₁
    mul!(du1,mass_matrix,du2)
    @.. linsolve_tmp = fsallast - du1 + c₃₂*f₁ + 2fsalfirst + dt*dT
  end

  cache.linsolve(vec(k₃), W, vec(linsolve_tmp))
  @.. k₃ = -k₃
  integrator.destats.nsolve += 1

  @.. u = uprev + dto6*(k₁ + 4k₂ + k₃)

  if integrator.opts.adaptive
    @.. tmp = dto6*(k₁ - 2*k₂ + k₃)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end

@muladd function perform_step!(integrator, cache::Rosenbrock23ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c₃₂,d,tf,uf = cache

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  mass_matrix = integrator.f.mass_matrix

  # Time derivative
  dT = calc_tderivative(integrator, cache)

  W = calc_W(integrator, cache, γ, repeat_step)
  k₁ = _reshape(W\-_vec((integrator.fsalfirst + γ*dT)), axes(uprev))
  integrator.destats.nsolve += 1
  f₁ = f(uprev  + dto2*k₁, p, t+dto2)
  integrator.destats.nf += 1

  if mass_matrix == I
    k₂ = _reshape(W\-_vec(f₁-k₁), axes(uprev)) + k₁
  else
    k₂ = _reshape(W\-_vec(f₁-mass_matrix*k₁), axes(uprev)) + k₁
  end
  integrator.destats.nsolve += 1
  u = uprev  + dt*k₂

  if integrator.opts.adaptive
    integrator.fsallast = f(u, p, t+dt)
    integrator.destats.nf += 1

    if mass_matrix == I
      k₃ = _reshape(W\-_vec((integrator.fsallast - c₃₂*(k₂-f₁) - 2*(k₁-integrator.fsalfirst) + dt*dT)), axes(uprev))
    else
      linsolve_tmp = integrator.fsallast - mass_matrix*(c₃₂*k₂ + 2*k₁) +c₃₂*f₁ + 2*integrator.fsalfirst + dt*dT
      k₃ = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
    end
    integrator.destats.nsolve += 1

    utilde =  dto6*(k₁ - 2*k₂ + k₃)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  integrator.u = u
end

@muladd function perform_step!(integrator, cache::Rosenbrock32ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c₃₂,d,tf,uf = cache

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  mass_matrix = integrator.f.mass_matrix

  # Time derivative
  dT = calc_tderivative(integrator, cache)

  W = calc_W(integrator, cache, γ, repeat_step)

  #f₀ = f(uprev, p, t)

  k₁ = _reshape(W\-_vec((integrator.fsalfirst + γ*dT)), axes(uprev))
  integrator.destats.nsolve += 1
  f₁ = f(uprev  + dto2*k₁, p, t+dto2)
  integrator.destats.nf += 1

  if mass_matrix == I
    k₂ = _reshape(W\-_vec(f₁-k₁), axes(uprev)) + k₁
  else
    linsolve_tmp = f₁ - mass_matrix * k₁
    k₂ = _reshape(W\-_vec(linsolve_tmp), axes(uprev)) + k₁
  end

  integrator.destats.nsolve += 1
  tmp = uprev  + dt*k₂
  integrator.fsallast = f(tmp, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    k₃ = _reshape(W\-_vec((integrator.fsallast - c₃₂*(k₂-f₁) - 2(k₁-integrator.fsalfirst) + dt*dT)), axes(uprev))
  else
    linsolve_tmp = integrator.fsallast - mass_matrix*(c₃₂*k₂ + 2k₁) + c₃₂*f₁ + 2*integrator.fsalfirst + dt*dT
    k₃ = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  end
  integrator.destats.nsolve += 1
  u = uprev  + dto6*(k₁ + 4k₂ + k₃)

  if integrator.opts.adaptive
    utilde =  dto6*(k₁ - 2k₂ + k₃)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = k₁
  integrator.k[2] = k₂
  integrator.u = u
end

function initialize!(integrator, cache::Union{Rosenbrock33ConstantCache,
                                              Rosenbrock34ConstantCache,
                                              Rosenbrock4ConstantCache,
                                              Rosenbrock5ConstantCache})
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::Union{Rosenbrock33Cache,
                                              Rosenbrock34Cache,
                                              Rosenbrock4Cache,
                                              Rosenbrock5Cache})
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [fsalfirst,fsallast]
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::Rosenbrock33ConstantCache,
                               repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,C21,C31,C32,b1,b2,b3,btilde1,btilde2,btilde3,gamma,c2,c3,d1,d2,d3 = cache.tab

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtgamma = dt*gamma

  mass_matrix = integrator.f.mass_matrix

  # Time derivative
  dT = calc_tderivative(integrator, cache)

  W = calc_W(integrator, cache, dtgamma, repeat_step, true)

  linsolve_tmp =  integrator.fsalfirst + dtd1*dT

  k1 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a21*k1
  du = f(u, p, t+c2*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd2*dT + dtC21*k1
  else
    linsolve_tmp =  du + dtd2*dT + mass_matrix * (dtC21*k1)
  end

  k2 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a31*k1 + a32*k2
  du = f(u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    linsolve_tmp =  du + dtd3*dT + mass_matrix * (dtC31*k1 + dtC32*k2)
  end

  k3 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + b1*k1 + b2*k2 + b3*k3
  integrator.fsallast = f(u, p, t + dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    utilde =  btilde1*k1 + btilde2*k2 + btilde3*k3
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

@muladd function perform_step!(integrator, cache::Rosenbrock33Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack du,du1,du2,fsalfirst,fsallast,k1,k2,k3,dT,J,W,uf,tf,linsolve_tmp,jac_config,atmp = cache
  @unpack a21,a31,a32,C21,C31,C32,b1,b2,b3,btilde1,btilde2,btilde3,gamma,c2,c3,d1,d2,d3 = cache.tab

  # Assignments
  mass_matrix = integrator.f.mass_matrix
  sizeu  = size(u)
  utilde = du

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtgamma = dt*gamma

  calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

  cache.linsolve(vec(k1), W, vec(linsolve_tmp), !repeat_step)
  @.. k1 = -k1
  integrator.destats.nsolve += 1

  @.. u = uprev + a21*k1
  f( du,  u, p, t+c2*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @.. du1 = dtC21*k1
    mul!(du2,mass_matrix,du1)
    @.. linsolve_tmp = du + dtd2*dT + du2
  end

  cache.linsolve(vec(k2), W, vec(linsolve_tmp))
  @.. k2 = -k2
  integrator.destats.nsolve += 1

  @.. u = uprev + a31*k1 + a32*k2
  f( du,  u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    @.. du1 = dtC31*k1 + dtC32*k2
    mul!(du2,mass_matrix,du1)
    @.. linsolve_tmp = du + dtd3*dT + du2
  end

  cache.linsolve(vec(k3), W, vec(linsolve_tmp))
  @.. k3 = -k3
  integrator.destats.nsolve += 1

  @.. u = uprev + b1*k1 + b2*k2 + b3*k3
  f( fsallast,  u, p, t + dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    @.. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end

################################################################################

@muladd function perform_step!(integrator, cache::Rosenbrock34ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt
  dtC41 = C41/dt
  dtC42 = C42/dt
  dtC43 = C43/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtd4 = dt*d4
  dtgamma = dt*gamma

  mass_matrix = integrator.f.mass_matrix
  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf, t)

  W = calc_W(integrator, cache, dtgamma, repeat_step, true)

  linsolve_tmp =  integrator.fsalfirst + dtd1*dT

  k1 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev # +a21*k1 a21 == 0
  # du = f(u, p, t+c2*dt) c2 == 0 and a21 == 0 => du = f(uprev, p, t) == fsalfirst

  if mass_matrix == I
    linsolve_tmp =  integrator.fsalfirst + dtd2*dT + dtC21*k1
  else
    linsolve_tmp =  integrator.fsalfirst + dtd2*dT + mass_matrix * (dtC21*k1)
  end

  k2 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a31*k1 + a32*k2
  du = f(u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    linsolve_tmp =  du + dtd3*dT + mass_matrix * (dtC31*k1 + dtC32*k2)
  end

  k3 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
  else
    linsolve_tmp =  du + dtd4*dT + mass_matrix * (dtC41*k1 + dtC42*k2 + dtC43*k3)
  end

  k4 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + b1*k1 + b2*k2 + b3*k3 + b4*k4
  integrator.fsallast = f(u, p, t + dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    utilde =  btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

@muladd function perform_step!(integrator, cache::Rosenbrock34Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack du,du1,du2,fsalfirst,fsallast,k1,k2,k3,k4,dT,J,W,uf,tf,linsolve_tmp,jac_config,atmp = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Assignments
  uidx = eachindex(integrator.uprev)
  sizeu  = size(u)
  mass_matrix = integrator.f.mass_matrix
  utilde = du

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt
  dtC41 = C41/dt
  dtC42 = C42/dt
  dtC43 = C43/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtd4 = dt*d4
  dtgamma = dt*gamma

  calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

  cache.linsolve(vec(k1), W, vec(linsolve_tmp), !repeat_step)
  @.. k1 = -k1
  integrator.destats.nsolve += 1

  #=
  a21 == 0 and c2 == 0
  so du = fsalfirst!
  @.. u = uprev + a21*k1

  f(du, u, p, t+c2*dt)
  =#

  if mass_matrix == I
    @.. linsolve_tmp = fsalfirst + dtd2*dT + dtC21*k1
  else
    @.. du1 = dtC21*k1
    mul!(du2,mass_matrix,du1)
    @.. linsolve_tmp = fsalfirst + dtd2*dT + du2
  end

  cache.linsolve(vec(k2), W, vec(linsolve_tmp))
  @.. k2 = -k2
  integrator.destats.nsolve += 1

  @.. u = uprev + a31*k1 + a32*k2
  f( du,  u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    @.. du1 = dtC31*k1 + dtC32*k2
    mul!(du2,mass_matrix,du1)
    @.. linsolve_tmp = du + dtd3*dT + du2
  end

  cache.linsolve(vec(k3), W, vec(linsolve_tmp))
  @.. k3 = -k3
  integrator.destats.nsolve += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
  else
    @.. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    mul!(du2,mass_matrix,du1)
    @.. linsolve_tmp = du + dtd4*dT + du2
  end

  cache.linsolve(vec(k4), W, vec(linsolve_tmp))
  @.. k4 = -k4
  integrator.destats.nsolve += 1

  @.. u = uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4
  f( fsallast,  u, p, t + dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    @.. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end

################################################################################

#### ROS34PW type method

@ROS34PW(:init)
@ROS34PW(:performstep)

################################################################################

#### ROS4 type method

@Rosenbrock4(:performstep)

################################################################################

#### Rodas4 type method

function initialize!(integrator, cache::Rodas4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  # Avoid undefined entries if k is an array of arrays
  integrator.k[1] = zero(integrator.u)
  integrator.k[2] = zero(integrator.u)
end

@muladd function perform_step!(integrator, cache::Rodas4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt
  dtC41 = C41/dt
  dtC42 = C42/dt
  dtC43 = C43/dt
  dtC51 = C51/dt
  dtC52 = C52/dt
  dtC53 = C53/dt
  dtC54 = C54/dt
  dtC61 = C61/dt
  dtC62 = C62/dt
  dtC63 = C63/dt
  dtC64 = C64/dt
  dtC65 = C65/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtd4 = dt*d4
  dtgamma = dt*gamma

  mass_matrix = integrator.f.mass_matrix

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf, t)

  W = calc_W(integrator, cache, dtgamma, repeat_step, true)

  du = f(uprev, p, t)
  integrator.destats.nf += 1

  linsolve_tmp =  du + dtd1*dT

  k1 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a21*k1
  du = f(u, p, t+c2*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd2*dT + dtC21*k1
  else
    linsolve_tmp =  du + dtd2*dT + mass_matrix * (dtC21*k1)
  end

  k2 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a31*k1 + a32*k2
  du = f(u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd3*dT + (dtC31*k1 + dtC32*k2)
  else
    linsolve_tmp =  du + dtd3*dT + mass_matrix * (dtC31*k1 + dtC32*k2)
  end

  k3 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a41*k1 + a42*k2 + a43*k3
  du = f(u, p, t+c4*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
  else
    linsolve_tmp =  du + dtd4*dT + mass_matrix * (dtC41*k1 + dtC42*k2 + dtC43*k3)
  end

  k4 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a51*k1 + a52*k2 + a53*k3 + a54*k4
  du = f(u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  else
    linsolve_tmp =  du + mass_matrix * (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  end

  k5 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = u + k5
  du = f(u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)
  else
    linsolve_tmp =  du + mass_matrix * (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)
  end

  k6 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = u + k6

  if integrator.opts.adaptive
    atmp = calculate_residuals(k6, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  if integrator.opts.calck
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    integrator.k[1] =  h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    integrator.k[2] =  h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
  end
  integrator.u = u
end


function initialize!(integrator, cache::Rodas4Cache)
  integrator.kshortsize = 2
  @unpack dense1,dense2 = cache
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [dense1,dense2]
end

@muladd function perform_step!(integrator, cache::Rodas4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack du,du1,du2,dT,J,W,uf,tf,k1,k2,k3,k4,k5,k6,linsolve_tmp,jac_config,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

  # Assignments
  sizeu  = size(u)
  uidx = eachindex(integrator.uprev)
  mass_matrix = integrator.f.mass_matrix

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt
  dtC41 = C41/dt
  dtC42 = C42/dt
  dtC43 = C43/dt
  dtC51 = C51/dt
  dtC52 = C52/dt
  dtC53 = C53/dt
  dtC54 = C54/dt
  dtC61 = C61/dt
  dtC62 = C62/dt
  dtC63 = C63/dt
  dtC64 = C64/dt
  dtC65 = C65/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtd4 = dt*d4
  dtgamma = dt*gamma

  calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

  cache.linsolve(vec(k1), W, vec(linsolve_tmp), !repeat_step)
  @.. k1 = -k1
  integrator.destats.nsolve += 1

  @.. u = uprev + a21*k1
  f( du,  u, p, t+c2*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @.. du1 = dtC21*k1
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd2*dT + du2
  end

  cache.linsolve(vec(k2), W, vec(linsolve_tmp))
  @.. k2 = -k2
  integrator.destats.nsolve += 1

  @.. u = uprev + a31*k1 + a32*k2
  f( du,  u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd3*dT + (dtC31*k1 + dtC32*k2)
  else
    @.. du1 = dtC31*k1 + dtC32*k2
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd3*dT + du2
  end

  cache.linsolve(vec(k3), W, vec(linsolve_tmp))
  @.. k3 = -k3
  integrator.destats.nsolve += 1

  @.. u = uprev + a41*k1 + a42*k2 + a43*k3
  f( du,  u, p, t+c4*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
  else
    @.. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd4*dT + du2
  end

  cache.linsolve(vec(k4), W, vec(linsolve_tmp))
  @.. k4 = -k4
  integrator.destats.nsolve += 1

  @.. u = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
  f( du,  u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  else
    @.. du1 = dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + du2
  end

  cache.linsolve(vec(k5), W, vec(linsolve_tmp))
  @.. k5 = -k5
  integrator.destats.nsolve += 1

  u .+= k5
  f( du,  u, p, t + dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)
  else
    @.. du1 = dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + du2
  end

  cache.linsolve(vec(k6), W, vec(linsolve_tmp))
  @.. k6 = -k6
  integrator.destats.nsolve += 1

  u .+= k6

  if integrator.opts.adaptive
    calculate_residuals!(atmp, k6, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  if integrator.opts.calck
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    @.. integrator.k[1] = h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    @.. integrator.k[2] = h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
  end
end

###############################################################################

### Rodas5 Method

@muladd function perform_step!(integrator, cache::Rosenbrock5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,C71,C72,C73,C74,C75,C76,C81,C82,C83,C84,C85,C86,C87,gamma,d1,d2,d3,d4,d5,c2,c3,c4,c5 = cache.tab

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt
  dtC41 = C41/dt
  dtC42 = C42/dt
  dtC43 = C43/dt
  dtC51 = C51/dt
  dtC52 = C52/dt
  dtC53 = C53/dt
  dtC54 = C54/dt
  dtC61 = C61/dt
  dtC62 = C62/dt
  dtC63 = C63/dt
  dtC64 = C64/dt
  dtC65 = C65/dt
  dtC71 = C71/dt
  dtC72 = C72/dt
  dtC73 = C73/dt
  dtC74 = C74/dt
  dtC75 = C75/dt
  dtC76 = C76/dt
  dtC81 = C81/dt
  dtC82 = C82/dt
  dtC83 = C83/dt
  dtC84 = C84/dt
  dtC85 = C85/dt
  dtC86 = C86/dt
  dtC87 = C87/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtd4 = dt*d4
  dtd5 = dt*d5
  dtgamma = dt*gamma

  mass_matrix = integrator.f.mass_matrix

  # Time derivative
  dT = calc_tderivative(integrator, cache)

  W = calc_W(integrator, cache, dtgamma, repeat_step, true)

  du1 = f(uprev, p, t)
  integrator.destats.nf += 1

  linsolve_tmp =  du1 + dtd1*dT

  k1 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a21*k1
  du = f(u, p, t+c2*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd2*dT + dtC21*k1
  else
    linsolve_tmp = du + dtd2*dT + mass_matrix * (dtC21*k1)
  end

  k2 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a31*k1 + a32*k2
  du = f(u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd3*dT + (dtC31*k1 + dtC32*k2)
  else
    linsolve_tmp =  du + dtd3*dT + mass_matrix * (dtC31*k1 + dtC32*k2)
  end

  k3 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a41*k1 + a42*k2 + a43*k3
  du = f(u, p, t+c4*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
  else
    linsolve_tmp =  du + dtd4*dT + mass_matrix * (dtC41*k1 + dtC42*k2 + dtC43*k3)
  end

  k4 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a51*k1 + a52*k2 + a53*k3 + a54*k4
  du = f(u, p, t+c5*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp = du + dtd5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  else
    linsolve_tmp = du + dtd5*dT + mass_matrix * (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  end

  k5 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = uprev  + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5
  du = f(u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp =  du + (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)
  else
    linsolve_tmp =  du + mass_matrix * (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)
  end

  k6 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = u + k6
  du = f(u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp = du + (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
  else
    linsolve_tmp = du + mass_matrix * (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
  end

  k7 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = u + k7
  du = f(u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    linsolve_tmp = du + (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
  else
    linsolve_tmp = du + mass_matrix * (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
  end

  k8 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
  integrator.destats.nsolve += 1
  u = u + k8
  du = f(u, p, t+dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    atmp = calculate_residuals(k8, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  if integrator.opts.calck
    #=
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    integrator.k[1] =  h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    integrator.k[2] =  h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
    =#
    integrator.k[1] = du1
    integrator.k[2] = du
  end

  integrator.fsallast = du
  integrator.u = u
end

@muladd function perform_step!(integrator, cache::Rosenbrock5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack du,du1,du2,fsalfirst,fsallast,k1,k2,k3,k4,k5,k6,k7,k8,dT,J,W,uf,tf,linsolve_tmp,jac_config,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,C71,C72,C73,C74,C75,C76,C81,C82,C83,C84,C85,C86,C87,gamma,d1,d2,d3,d4,d5,c2,c3,c4,c5 = cache.tab

  # Assignments
  sizeu  = size(u)
  uidx = eachindex(integrator.uprev)
  mass_matrix = integrator.f.mass_matrix

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt
  dtC41 = C41/dt
  dtC42 = C42/dt
  dtC43 = C43/dt
  dtC51 = C51/dt
  dtC52 = C52/dt
  dtC53 = C53/dt
  dtC54 = C54/dt
  dtC61 = C61/dt
  dtC62 = C62/dt
  dtC63 = C63/dt
  dtC64 = C64/dt
  dtC65 = C65/dt
  dtC71 = C71/dt
  dtC72 = C72/dt
  dtC73 = C73/dt
  dtC74 = C74/dt
  dtC75 = C75/dt
  dtC76 = C76/dt
  dtC81 = C81/dt
  dtC82 = C82/dt
  dtC83 = C83/dt
  dtC84 = C84/dt
  dtC85 = C85/dt
  dtC86 = C86/dt
  dtC87 = C87/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtd4 = dt*d4
  dtd5 = dt*d5
  dtgamma = dt*gamma

  calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

  cache.linsolve(vec(k1), W, vec(linsolve_tmp), !repeat_step)
  @.. k1 = -k1
  integrator.destats.nsolve += 1

  @.. u = uprev + a21*k1
  f( du,  u, p, t+c2*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @.. du1 = dtC21*k1
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd2*dT + du2
  end

  cache.linsolve(vec(k2), W, vec(linsolve_tmp))
  @.. k2 = -k2
  integrator.destats.nsolve += 1

  @.. u = uprev + a31*k1 + a32*k2
  f( du,  u, p, t+c3*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd3*dT + (dtC31*k1 + dtC32*k2)
  else
    @.. du1 = dtC31*k1 + dtC32*k2
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd3*dT + du2
  end

  cache.linsolve(vec(k3), W, vec(linsolve_tmp))
  @.. k3 = -k3
  integrator.destats.nsolve += 1

  @.. u = uprev + a41*k1 + a42*k2 + a43*k3
  f( du,  u, p, t+c4*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
  else
    @.. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd4*dT + du2
  end

  cache.linsolve(vec(k4), W, vec(linsolve_tmp))
  @.. k4 = -k4
  integrator.destats.nsolve += 1

  @.. u = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
  f( du,  u, p, t+c5*dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + dtd5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  else
    @.. du1 = dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + dtd5*dT + du2
  end

  cache.linsolve(vec(k5), W, vec(linsolve_tmp))
  @.. k5 = -k5
  integrator.destats.nsolve += 1

  @.. u = uprev + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5
  f( du,  u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)
  else
    @.. du1 = dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + du2
  end

  cache.linsolve(vec(k6), W, vec(linsolve_tmp))
  @.. k6 = -k6
  integrator.destats.nsolve += 1

  u .+= k6
  f( du,  u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
  else
    @.. du1 = dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + du2
  end

  cache.linsolve(vec(k7), W, vec(linsolve_tmp))
  @.. k7 = -k7
  integrator.destats.nsolve += 1

  u .+= k7
  f( du,  u, p, t+dt)
  integrator.destats.nf += 1

  if mass_matrix == I
    @.. linsolve_tmp = du + (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
  else
    @.. du1 = dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7
    mul!(vec(du2),mass_matrix,vec(du1))
    @.. linsolve_tmp = du + du2
  end

  cache.linsolve(vec(k8), W, vec(linsolve_tmp))
  @.. k8 = -k8
  integrator.destats.nsolve += 1

  u .+= k8
  f( fsallast,  u, p, t+dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    calculate_residuals!(atmp, k8, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  if integrator.opts.calck
    #=
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    @.. integrator.k[1] = h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    @.. integrator.k[2] = h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
    =#
  end
end

@RosenbrockW6S4OS(:init)
@RosenbrockW6S4OS(:performstep)
