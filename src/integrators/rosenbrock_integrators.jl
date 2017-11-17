function initialize!(integrator, cache::Rosenbrock23Cache)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [k₁,k₂]
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock23Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack c₃₂,d = cache.tab

  # Assignments
  sizeu  = size(u)
  mass_matrix = integrator.sol.prob.mass_matrix

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad}, t, uprev, dT)
    else
      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  @. linsolve_tmp = fsalfirst + γ*dT

  # Jacobian
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW}, t, uprev, γ, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev, J)
      else

        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j] - γ*J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  recursivecopy!(k₁, reshape(vectmp, size(u)...))
  @. u = uprev + dto2*k₁
  f(t+dto2, u, f₁)

  if mass_matrix == I
    tmp .= k₁
  else
    A_mul_B!(tmp,mass_matrix,k₁)
  end

  @. linsolve_tmp = f₁ - tmp
  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  tmp2 = reshape(vectmp2, sizeu...)
  @. k₂ = tmp2 + k₁
  @. u = uprev + dt*k₂

  if integrator.opts.adaptive
    f(t+dt, u, fsallast)

    if mass_matrix == I
      @. linsolve_tmp = fsallast - c₃₂*(k₂-f₁) - 2(k₁-fsalfirst) + dt*dT
    else
      @. du2 = c₃₂*k₂ + 2k₁
      A_mul_B!(du1,mass_matrix,du2)
      @. linsolve_tmp = fsallast - du1 + c₃₂*f₁ + 2fsalfirst + dt*dT
    end


    if has_invW(f)
      A_mul_B!(vectmp3, W, linsolve_tmp_vec)
    else
      integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
    end

    k₃ = reshape(vectmp3, sizeu...)

    @. tmp = dto6*(k₁ - 2*k₂ + k₃)
    # does not work with units - additional unitless array required!
    calculate_residuals!(tmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
end

function initialize!(integrator, cache::Rosenbrock32Cache)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [k₁,k₂]
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock32Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack c₃₂,d = cache.tab

  # Assignments
  sizeu  = size(u)
  mass_matrix = integrator.sol.prob.mass_matrix

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad}, t, uprev, dT)
    else

      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  @. linsolve_tmp = fsalfirst + γ*dT

  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW}, t, uprev, γ, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev, J)
      else

        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = mass_matrix[i,j] - γ*J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  recursivecopy!(k₁, reshape(vectmp, sizeu...))
  @. u = uprev + dto2*k₁
  f(t+dto2, u, f₁)

  if mass_matrix == I
    tmp .= k₁
  else
    A_mul_B!(tmp,mass_matrix,k₁)
  end

  @. linsolve_tmp = f₁ - tmp

  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  tmp2 = reshape(vectmp2, sizeu...)
  @. k₂ = tmp2 + k₁
  @. tmp = uprev + dt*k₂
  f(t+dt, tmp, fsallast)

  if mass_matrix == I
    @. linsolve_tmp = fsallast - c₃₂*(k₂-f₁) - 2(k₁-fsalfirst) + dt*dT
  else
    @. du2 = c₃₂*k₂ + 2k₁
    A_mul_B!(du1,mass_matrix,du2)
    @. linsolve_tmp = fsallast - du1 + c₃₂*f₁ + 2fsalfirst + dt*dT
  end

  if has_invW(f)
    A_mul_B!(vectmp3, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
  end

  k₃ = reshape(vectmp3, sizeu...)
  @. u = uprev + dto6*(k₁ + 4k₂ + k₃)

  if integrator.opts.adaptive
    @. tmp = dto6*(k₁ - 2*k₂ + k₃)
    # does not work with units - additional unitless array required!
    calculate_residuals!(tmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
end

function initialize!(integrator, cache::Rosenbrock23ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = zero(integrator.fsalfirst)
  integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock23ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack c₃₂,d,tf,uf = cache

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf, t)

  # Jacobian
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf, uprev)
    W = I - γ*J
  else
    J = ForwardDiff.derivative(uf, uprev)
    W = 1 - γ*J
  end

  k₁ = W\(@. integrator.fsalfirst + γ*dT)
  f₁ = f(t+dto2, @. uprev + dto2*k₁)

  k₂ = W\(f₁-k₁) + k₁
  u = @. uprev + dt*k₂

  if integrator.opts.adaptive
    integrator.fsallast = f(t+dt, u)

    k₃ = W\(@. integrator.fsallast - c₃₂*(k₂-f₁) - 2*(k₁-integrator.fsalfirst) + dt*dT)

    utilde = @. dto6*(k₁ - 2*k₂ + k₃)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = k₁
  integrator.k[2] = k₂
  integrator.u = u
end

function initialize!(integrator, cache::Rosenbrock32ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = zero(integrator.fsalfirst)
  integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock32ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack c₃₂,d,tf,uf = cache

  # Precalculations
  γ = dt*d
  dto2 = dt/2
  dto6 = dt/6

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf,t)

  # Jacobian
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - γ*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γ*J
  end

  #f₀ = f(t, uprev)

  k₁ = W\(@. integrator.fsalfirst + γ*dT)
  f₁ = f(t+dto2, @. uprev + dto2*k₁)

  k₂ = W\(f₁-k₁) + k₁
  tmp = @. uprev + dt*k₂
  integrator.fsallast = f(t+dt, tmp)

  k₃ = W\(@. integrator.fsallast - c₃₂*(k₂-f₁) - 2(k₁-integrator.fsalfirst) + dt*dT)
  u = @. uprev + dto6*(k₁ + 4k₂ + k₃)

  if integrator.opts.adaptive
    utilde = @. dto6*(k₁ - 2k₂ + k₃)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = k₁
  integrator.k[2] = k₂
  integrator.u = u
end

function initialize!(integrator, cache::Rosenbrock33ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Rosenbrock33ConstantCache,
                               repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf,t)

  # Jacobian
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I/dtgamma - J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1/dtgamma - J
  end

  linsolve_tmp = @. integrator.fsalfirst + dtd1*dT

  k1 = W\linsolve_tmp
  u = @. uprev + a21*k1
  du = f(t+c2*dt, u)

  linsolve_tmp = @. du + dtd2*dT + dtC21*k1

  k2 = W\linsolve_tmp
  u = @. uprev + a31*k1 + a32*k2
  du = f(t+c3*dt, u)

  linsolve_tmp = @. du + dtd3*dT + dtC31*k1 + dtC32*k2

  k3 = W\linsolve_tmp
  u = @. uprev + b1*k1 + b2*k2 + b3*k3
  integrator.fsallast = f(t, u)

  if integrator.opts.adaptive
    utilde = @. btilde1*k1 + btilde2*k2 + btilde3*k3
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Rosenbrock33Cache)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [fsalfirst,fsallast]
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock33Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,C21,C31,C32,b1,b2,b3,btilde1,btilde2,btilde3,gamma,c2,c3,d1,d2,d3 = cache.tab

  # Assignments
  mass_matrix = integrator.sol.prob.mass_matrix
  sizeu  = size(u)
  utilde = du
  atmp = du # does not work with units - additional unitless array required!

  # Precalculations
  dtC21 = C21/dt
  dtC31 = C31/dt
  dtC32 = C32/dt

  dtd1 = dt*d1
  dtd2 = dt*d2
  dtd3 = dt*d3
  dtgamma = dt*gamma

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad}, t, uprev,dT)
    else

      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  @. linsolve_tmp = fsalfirst + dtd1*dT

  # Jacobian
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW_t}, t, uprev, dtgamma, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev, J)
      else

        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  k1 = reshape(vectmp, sizeu...)
  @. u = uprev + a21*k1
  f(t+c2*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @. du1 = dtC21*k1
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd2*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2, sizeu...)
  @. u = uprev + a31*k1 + a32*k2
  f(t+c3*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    @. du1 = dtC31*k1 + dtC32*k2
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd3*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp3, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3, sizeu...)
  @. u = uprev + b1*k1 + b2*k2 + b3*k3
  f(t, u, fsallast)

  if integrator.opts.adaptive
    @. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

################################################################################

function initialize!(integrator, cache::Rosenbrock34ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Rosenbrock34ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf, t)

  # Jacobian
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf, uprev)
    W = I/dtgamma - J
  else
    J = ForwardDiff.derivative(uf, uprev)
    W = 1/dtgamma - J
  end

  linsolve_tmp = @. integrator.fsalfirst + dtd1*dT

  k1 = W\linsolve_tmp
  u = uprev # +a21*k1 a21 == 0
  # du = f(t+c2*dt,u) c2 == 0 and a21 == 0 => du = f(t,uprev) == fsalfirst

  linsolve_tmp = @. integrator.fsalfirst + dtd2*dT + dtC21*k1

  k2 = W\linsolve_tmp
  u = @. uprev + a31*k1 + a32*k2
  du = f(t+c3*dt, u)

  linsolve_tmp = @. du + dtd3*dT + dtC31*k1 + dtC32*k2

  k3 = W\linsolve_tmp
  linsolve_tmp = @. du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3

  k4 = W\linsolve_tmp
  u = @. uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4
  integrator.fsallast = f(t, u)

  if integrator.opts.adaptive
    utilde = @. btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Rosenbrock34Cache)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [fsalfirst,fsallast]
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock34Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Assignments
  uidx = eachindex(integrator.uprev)
  sizeu  = size(u)
  mass_matrix = integrator.sol.prob.mass_matrix
  utilde = du
  atmp = du # does not work with units - additional unitless array required!

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

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad},t,uprev,dT)
    else

      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  @. linsolve_tmp = fsalfirst + dtd1*dT

  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW_t}, t, uprev, dtgamma, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev,J)
      else

        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  k1 = reshape(vectmp, sizeu...)

  #=
  a21 == 0 and c2 == 0
  so du = fsalfirst!
  @. u = uprev + a21*k1

  f(t+c2*dt,u,du)
  =#

  if mass_matrix == I
    @. linsolve_tmp = fsalfirst + dtd2*dT + dtC21*k1
  else
    @. du1 = dtC21*k1
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = fsalfirst + dtd2*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2, sizeu...)
  @. u = uprev + a31*k1 + a32*k2
  f(t+c3*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    @. du1 = dtC31*k1 + dtC32*k2
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd3*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp3, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3, sizeu...)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
  else
    @. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd4*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp4, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4, W, linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4, sizeu...)
  @. u = uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4
  f(t, u, fsallast)

  if integrator.opts.adaptive
    @. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end


################################################################################

#### ROS4 type method

function initialize!(integrator, cache::Rosenbrock4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Rosenbrock4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf, t)

  # Jacobian
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf, uprev)
    W = I/dtgamma - J
  else
    J = ForwardDiff.derivative(uf, uprev)
    W = 1/dtgamma - J
  end

  linsolve_tmp = @. integrator.fsalfirst + dtd1*dT

  k1 = W\linsolve_tmp
  u = @. uprev+a21*k1
  du = f(t+c2*dt, u)

  linsolve_tmp = @. du + dtd2*dT + dtC21*k1

  k2 = W\linsolve_tmp
  u = @. uprev + a31*k1 + a32*k2
  du = f(t+c3*dt, u)

  linsolve_tmp = @. du + dtd3*dT + dtC31*k1 + dtC32*k2

  k3 = W\linsolve_tmp

  linsolve_tmp = @. du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3

  k4 = W\linsolve_tmp
  u = @. uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4
  integrator.fsallast = f(t,u)

  if integrator.opts.adaptive
    utilde = @. btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Rosenbrock4Cache)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [fsalfirst,fsallast]
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Assignments
  sizeu  = size(u)
  uidx = eachindex(integrator.uprev)
  mass_matrix = integrator.sol.prob.mass_matrix
  utilde = du
  atmp = du # does not work with units - additional unitless array required!

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

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad}, t, uprev, dT)
    else

      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  @. linsolve_tmp = fsalfirst + dtd1*dT

  # Jacobian
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW_t}, t, uprev, dtgamma, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev, J)
      else

        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  k1 = reshape(vectmp, sizeu...)
  @. u = uprev + a21*k1
  f(t+c2*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @. du1 = dtC21*k1
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd2*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2, sizeu...)
  @. u = uprev + a31*k1 + a32*k2
  f(t+c3*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd3*dT + dtC31*k1 + dtC32*k2
  else
    @. du1 = dtC31*k1 + dtC32*k2
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd3*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp3, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3, sizeu...)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
  else
    @. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd4*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp4, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4, W, linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4, sizeu...)
  @. u = uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4
  f(t, u, fsallast)

  if integrator.opts.adaptive
    @. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

################################################################################

#### Rodas4 type method

function initialize!(integrator, cache::Rodas4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  # Avoid undefined entries if k is an array of arrays
  integrator.k[1] = zero(integrator.u)
  integrator.k[2] = zero(integrator.u)
end

@muladd function perform_step!(integrator, cache::Rodas4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf, t)

  # Jacobian
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf, uprev)
    W = I/dtgamma - J
  else
    J = ForwardDiff.derivative(uf, uprev)
    W = 1/dtgamma - J
  end

  du = f(t, uprev)

  linsolve_tmp = @. du + dtd1*dT

  k1 = W\linsolve_tmp
  u = @. uprev + a21*k1
  du = f(t+c2*dt, u)

  linsolve_tmp = @. du + dtd2*dT + dtC21*k1

  k2 = W\linsolve_tmp
  u = @. uprev + a31*k1 + a32*k2
  du = f(t+c3*dt, u)

  linsolve_tmp = @. du + dtd3*dT + (dtC31*k1 + dtC32*k2)

  k3 = W\linsolve_tmp
  u = @. uprev + a41*k1 + a42*k2 + a43*k3
  du = f(t+c4*dt, u)

  linsolve_tmp = @. du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

  k4 = W\linsolve_tmp
  u = @. uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
  du = f(t+dt, u)

  linsolve_tmp = @. du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)

  k5 = W\linsolve_tmp
  u = u + k5
  du = f(t+dt, u)

  linsolve_tmp = @. du + (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)

  k6 = W\linsolve_tmp
  u = u + k6

  if integrator.opts.adaptive
    atmp = calculate_residuals(k6, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if integrator.opts.calck
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    integrator.k[1] = @. h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    integrator.k[2] = @. h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
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
  @unpack t,dt,uprev,u,f = integrator
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,vectmp5,vectmp6,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

  # Assignments
  sizeu  = size(u)
  uidx = eachindex(integrator.uprev)
  mass_matrix = integrator.sol.prob.mass_matrix
  atmp = du # does not work with units - additional unitless array required!

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

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad}, t, uprev, dT)
    else
      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  f(t, uprev, du)

  @. linsolve_tmp = du + dtd1*dT

  # Jacobian
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW_t}, t, uprev, dtgamma, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev, J)
      else

        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  k1 = reshape(vectmp, sizeu...)
  @. u = uprev + a21*k1
  f(t+c2*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @. du1 = dtC21*k1
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd2*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2, sizeu...)
  @. u = uprev + a31*k1 + a32*k2
  f(t+c3*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd3*dT + (dtC31*k1 + dtC32*k2)
  else
    @. du1 = dtC31*k1 + dtC32*k2
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd3*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp3, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3, sizeu...)
  @. u = uprev + a41*k1 + a42*k2 + a43*k3
  f(t+c4*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
  else
    @. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd4*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp4, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4, W, linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4, sizeu...)
  @. u = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
  f(t+dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  else
    @. du1 = dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp5, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp5, W, linsolve_tmp_vec)
  end

  k5 = reshape(vectmp5, sizeu...)
  u .+= k5
  f(t, u, du)

  if mass_matrix == I
    # @. linsolve_tmp = du + (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)
    @tight_loop_macros for i in uidx
      @inbounds linsolve_tmp[i] = du[i] + (dtC61*k1[i] + dtC62*k2[i] + dtC65*k5[i] + dtC64*k4[i] + dtC63*k3[i])
    end
  else
    # @. linsolve_tmp = dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3
    @tight_loop_macros for i in uidx
      @inbounds du1[i] = dtC61*k1[i] + dtC62*k2[i] + dtC65*k5[i] + dtC64*k4[i] + dtC63*k3[i]
    end
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp6, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp6, W, linsolve_tmp_vec)
  end

  k6 = reshape(vectmp6, sizeu...)
  u .+= k6

  if integrator.opts.adaptive
    calculate_residuals!(atmp, k6, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if integrator.opts.calck
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    @. integrator.k[1] = h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    @. integrator.k[2] = h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
  end
end

###############################################################################

### Rodas5 Method

function initialize!(integrator, cache::Rosenbrock5ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Rosenbrock5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

  # Time derivative
  tf.u = uprev
  dT = ForwardDiff.derivative(tf,t)

  # Jacobian
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf, uprev)
    W = I/dtgamma - J
  else
    J = ForwardDiff.derivative(uf, uprev)
    W = 1/dtgamma - J
  end

  du1 = f(t, uprev)

  linsolve_tmp = @. du1 + dtd1*dT

  k1 = W\linsolve_tmp
  u = @. uprev + a21*k1
  du = f(t+c2*dt, u)

  linsolve_tmp = @. du + dtd2*dT + dtC21*k1

  k2 = W\linsolve_tmp
  u = @. uprev + a31*k1 + a32*k2
  du = f(t+c3*dt, u)

  linsolve_tmp = @. du + dtd3*dT + (dtC31*k1 + dtC32*k2)

  k3 = W\linsolve_tmp
  u = @. uprev + a41*k1 + a42*k2 + a43*k3
  du = f(t+c4*dt, u)

  linsolve_tmp = @. du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

  k4 = W\linsolve_tmp
  u = @. uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
  du = f(t+c5*dt, u)

  # linsolve_tmp = @. du + dtd5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    tmp = similar(du)
    @tight_loop_macros for i in eachindex(tmp)
      @inbounds tmp[i] = du[i] + dtd5*dT[i] + (dtC52*k2[i] + dtC54*k4[i] + dtC51*k1[i] + dtC53*k3[i])
    end
    linsolve_tmp = tmp
  else
    linsolve_tmp = du + dtd5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  end

  k5 = W\linsolve_tmp
  u = @. uprev + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5
  du = f(t+dt, u)

  linsolve_tmp = @. du + (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)

  k6 = W\linsolve_tmp
  u = u + k6
  du = f(t+dt, u)

  # linsolve_tmp = @. du + (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    tmp = similar(du)
    @tight_loop_macros for i in eachindex(tmp)
      @inbounds tmp[i] = du[i] + (dtC71*k1[i] + dtC72*k2[i] + dtC73*k3[i] + dtC74*k4[i] + dtC75*k5[i] + dtC76*k6[i])
    end
    linsolve_tmp = tmp
  else
    linsolve_tmp = du + (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
  end

  k7 = W\linsolve_tmp
  u = u + k7
  du = f(t+dt, u)

  # linsolve_tmp = @. du + (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    tmp = similar(du)
    @tight_loop_macros for i in eachindex(tmp)
      @inbounds tmp[i] = du[i] + (dtC81*k1[i] + dtC82*k2[i] + dtC83*k3[i] + dtC84*k4[i] + dtC85*k5[i] + dtC86*k6[i] + dtC87*k7[i])
    end
    linsolve_tmp = tmp
  else
    linsolve_tmp = du + (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
  end

  k8 = W\linsolve_tmp
  u = u + k8
  du = f(t+dt, u)

  if integrator.opts.adaptive
    atmp = calculate_residuals(k8, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if integrator.opts.calck
    #=
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    integrator.k[1] = @. h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    integrator.k[2] = @. h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
    =#
    integrator.k[1] = du1
    integrator.k[2] = du
  end

  integrator.fsallast = du
  integrator.u = u
end

function initialize!(integrator, cache::Rosenbrock5Cache)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast,dense1,dense2 = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [fsalfirst,fsallast]
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,vectmp5,vectmp6,vectmp7,vectmp8,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,C71,C72,C73,C74,C75,C76,C81,C82,C83,C84,C85,C86,C87,gamma,d1,d2,d3,d4,d5,c2,c3,c4,c5 = cache.tab

  # Assignments
  sizeu  = size(u)
  uidx = eachindex(integrator.uprev)
  mass_matrix = integrator.sol.prob.mass_matrix
  atmp = du # does not work with units - additional unitless array required!

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

  # Time derivative
  if !repeat_step # skip calculation if step is repeated
    if has_tgrad(f)
      f(Val{:tgrad},t,uprev,dT)
    else
      tf.uprev = uprev
      derivative!(dT, tf, t, du2, integrator)
    end
  end

  f(t, uprev, fsalfirst)

  @. linsolve_tmp = fsalfirst + dtd1*dT

  # Jacobian
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW_t}, t, uprev, dtgamma, W) # W == inverse W

    A_mul_B!(vectmp, W, linsolve_tmp_vec)
  else
    if !repeat_step # skip calculation of J and W if step is repeated
      if has_jac(f)
        f(Val{:jac}, t, uprev, J)
      else
        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
      end
    end

    # use existing factorization of W if step is repeated
    integrator.alg.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
  end

  k1 = reshape(vectmp, sizeu...)
  @. u = uprev + a21*k1
  f(t+c2*dt, u, du)



  if mass_matrix == I
    @. linsolve_tmp = du + dtd2*dT + dtC21*k1
  else
    @. du1 = dtC21*k1
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd2*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp2, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2, W, linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2, sizeu...)
  @. u = uprev + a31*k1 + a32*k2
  f(t+c3*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd3*dT + (dtC31*k1 + dtC32*k2)
  else
    @. du1 = dtC31*k1 + dtC32*k2
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd3*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp3, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3, W, linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3, sizeu...)
  @. u = uprev + a41*k1 + a42*k2 + a43*k3
  f(t+c4*dt, u, du)

  if mass_matrix == I
    @. linsolve_tmp = du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
  else
    @. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd4*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp4, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4, W, linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4, sizeu...)
  @. u = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
  f(t+c5*dt, u, du)

  if mass_matrix == I
    #  @. linsolve_tmp = du + dtd5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
      @tight_loop_macros for i in uidx
        @inbounds linsolve_tmp[i] = du[i] + dtd5*dT[i] + (dtC52*k2[i] + dtC54*k4[i] + dtC51*k1[i] + dtC53*k3[i])
      end
  else
    @. du1 = dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + dtd5*dT + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp5, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp5, W, linsolve_tmp_vec)
  end

  k5 = reshape(vectmp5, sizeu...)
  # @. u = uprev + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]
  end
  f(t+dt, u, du)

  if mass_matrix == I
    # @. linsolve_tmp = du + (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)
    @tight_loop_macros for i in uidx
      @inbounds linsolve_tmp[i] = du[i] + (dtC61*k1[i] + dtC62*k2[i] + dtC63*k3[i] + dtC64*k4[i] + dtC65*k5[i])
    end
  else
    # @. du1 = dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5
    @tight_loop_macros for i in uidx
      @inbounds du1[i] = dtC61*k1[i] + dtC62*k2[i] + dtC63*k3[i] + dtC64*k4[i] + dtC65*k5[i]
    end
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp6, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp6, W, linsolve_tmp_vec)
  end

  k6 = reshape(vectmp6, sizeu...)
  u .+= k6
  f(t+dt, u, du)

  if mass_matrix == I
    # @. linsolve_tmp = du + (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
    @tight_loop_macros for i in uidx
      @inbounds linsolve_tmp[i] = du[i] + (dtC71*k1[i] + dtC72*k2[i] + dtC73*k3[i] + dtC74*k4[i] + dtC75*k5[i] + dtC76*k6[i])
    end
  else
    # @. du1 =dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6
    @tight_loop_macros for i in uidx
      @inbounds du1[i] = dtC71*k1[i] + dtC72*k2[i] + dtC73*k3[i] + dtC74*k4[i] + dtC75*k5[i] + dtC76*k6[i]
    end
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp7, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp7, W, linsolve_tmp_vec)
  end

  k7 = reshape(vectmp7, sizeu...)
  u .+= k7
  f(t+dt, u, du)

  if mass_matrix == I
    # @. linsolve_tmp = du + (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
    @tight_loop_macros for i in uidx
      @inbounds linsolve_tmp[i] = du[i] + (dtC81*k1[i] + dtC82*k2[i] + dtC83*k3[i] + dtC84*k4[i] + dtC85*k5[i] + dtC86*k6[i] + dtC87*k7[i])
    end
  else
    # @. du1 = dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7
    @tight_loop_macros for i in uidx
      @inbounds du1[i] = dtC81*k1[i] + dtC82*k2[i] + dtC83*k3[i] + dtC84*k4[i] + dtC85*k5[i] + dtC86*k6[i] + dtC87*k7[i]
    end
    A_mul_B!(du2,mass_matrix,du1)
    @. linsolve_tmp = du + du2
  end

  if has_invW(f)
    A_mul_B!(vectmp8, W, linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp8, W, linsolve_tmp_vec)
  end

  k8 = reshape(vectmp8, sizeu...)
  u .+= k8
  f(t+dt, u, fsallast)

  if integrator.opts.adaptive
    calculate_residuals!(atmp, k8, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if integrator.opts.calck
    #=
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    @. integrator.k[1] = h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    @. integrator.k[2] = h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
    =#
  end
end
