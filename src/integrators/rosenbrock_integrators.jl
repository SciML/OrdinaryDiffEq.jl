function initialize!(integrator,cache::Rosenbrock23Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [k₁,k₂]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock23Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack c₃₂,d = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  γ = dt*d

  @. linsolve_tmp = fsalfirst + γ*dT

  if has_invW(f)
    f(Val{:invW},t,uprev,γ,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
      @inbounds W[i,j] = mass_matrix[i,j]-γ*J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  recursivecopy!(k₁,reshape(vectmp,size(u)...))
  dto2 = dt/2

  @. u= uprev+dto2*k₁
  f(t+dto2,u,f₁)
  @. linsolve_tmp = f₁-k₁

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  tmp2 = reshape(vectmp2,sizeu...)
  @. k₂ = tmp2 + k₁
  @. u = uprev + dt*k₂

  if integrator.opts.adaptive
    f(t+dt,u,fsallast)

    @. linsolve_tmp = fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT

    if has_invW(f)
      A_mul_B!(vectmp3,W,linsolve_tmp_vec)
    else
      integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
    end

    k₃ = reshape(vectmp3,sizeu...)
    #@. tmp = (dt*(k₁ - 2k₂ + k₃)/6)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds tmp[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
end

function initialize!(integrator,cache::Rosenbrock32Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [k₁,k₂]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end#

@muladd function perform_step!(integrator,cache::Rosenbrock32Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  mass_matrix = integrator.sol.prob.mass_matrix
  @unpack c₃₂,d = cache.tab
  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  γ = dt*d

  @. linsolve_tmp = fsalfirst + γ*dT

  if has_invW(f)
    f(Val{:invW},t,uprev,γ,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
      @inbounds W[i,j] = mass_matrix[i,j]-γ*J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  recursivecopy!(k₁,reshape(vectmp,sizeu...))

  dto2 = dt/2
  @. u= uprev+dto2*k₁
  f(t+dto2,u,f₁)
  @. linsolve_tmp = f₁-k₁

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  tmp2 = reshape(vectmp2,sizeu...)
  @. k₂ = tmp2 + k₁
  @. tmp = uprev + dt*k₂

  f(t+dt,tmp,fsallast)
  @. linsolve_tmp = fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k₃ = reshape(vectmp3,sizeu...)
  @. u = uprev + dt*(k₁ + 4k₂ + k₃)/6

  if integrator.opts.adaptive
    #@. tmp = (dt*(k₁ - 2k₂ + k₃)/6)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds tmp[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
end

function initialize!(integrator,cache::Rosenbrock23ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = zero(integrator.fsalfirst)
  integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock23ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c₃₂,d,tf,uf = cache
  # Time derivative
  tf.u = uprev
  uf.t = t
  dT = ForwardDiff.derivative(tf,t)
  a = -dt*d
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I + a*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 + a*J
  end
  a = -a
  k₁ = W\(@. integrator.fsalfirst + a*dT)
  dto2 = dt/2
  f₁ = f(t+dto2,@. uprev+dto2*k₁)
  k₂ = W\(f₁-k₁) + k₁
  u = @. uprev + dt*k₂
  if integrator.opts.adaptive
    integrator.fsallast = f(t+dt,u)
    k₃ = W\(@. integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
    tmp = @. (dt*(k₁ - 2k₂ + k₃)/6)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  integrator.u = u
end

function initialize!(integrator,cache::Rosenbrock32ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = zero(integrator.fsalfirst)
  integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock32ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c₃₂,d,tf,uf = cache
  tf.u = uprev
  uf.t = t
  # Time derivative
  dT = ForwardDiff.derivative(tf,t)
  a = -dt*d
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I + a*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 + a*J
  end
  #f₀ = f(t,uprev)
  a = -a
  k₁ = W\(@. integrator.fsalfirst + a*dT)
  dto2 = dt/2
  f₁ = f(t+dto2,@. uprev+dto2*k₁)
  k₂ = W\(f₁-k₁) + k₁
  tmp = @. uprev + dt*k₂
  integrator.fsallast = f(t+dt,tmp)
  k₃ = W\(@. integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
  u = @. uprev + dt*(k₁ + 4k₂ + k₃)/6
  if integrator.opts.adaptive
    tmp = @. (dt*(k₁ - 2k₂ + k₃)/6)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  integrator.u = u
end

function initialize!(integrator,cache::Rosenbrock33ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Rosenbrock33ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,C21,C31,C32,b1,b2,b3,btilde1,btilde2,btilde3,gamma,c2,c3,d1,d2,d3 = cache.tab

  # Setup Jacobian Calc
  tf.u = uprev
  uf.t = t

  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I/(dt*gamma) - J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1/(dt*gamma) - J
  end

  d1dt = dt*d1
  linsolve_tmp = @. integrator.fsalfirst + d1dt*dT

  k1 = W\linsolve_tmp
  u = @. uprev+a21*k1
  du = f(t+c2*dt,u)

  linsolve_tmp = @. du + dt*d2*dT + C21*k1/dt

  k2 = W\linsolve_tmp

  u = @. uprev + a31*k1 + a32*k2

  du = f(t+c3*dt,u)

  linsolve_tmp = @. du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  k3 = W\linsolve_tmp

  u = @. uprev + b1*k1 + b2*k2 + b3*k3

  integrator.fsallast = f(t,u)

  if integrator.opts.adaptive
    utilde = @. (b1-btilde1)*k1 + (b2-btilde2)*k2 + (b3-btilde3)*k3
    atmp = @. utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::Rosenbrock33Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [fsalfirst,fsallast]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock33Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,C21,C31,C32,b1,b2,b3,btilde1,btilde2,btilde3,gamma,c2,c3,d1,d2,d3 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  utilde = du
  atmp = du

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  d1dt = dt*d1
  @. linsolve_tmp = fsalfirst + d1dt*dT

  if has_invW(f)
    f(Val{:invW_t},t,uprev,dt*gamma,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
        @inbounds W[i,j] = mass_matrix[i,j]/(dt*gamma)-J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  k1 = reshape(vectmp,sizeu...)

  @. u = uprev+a21*k1

  f(t+c2*dt,u,du)

  @. linsolve_tmp = du + dt*d2*dT + C21*k1/dt

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2,sizeu...)

  @. u = uprev + a31*k1 + a32*k2

  f(t+c3*dt,u,du)

  @. linsolve_tmp = du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3,sizeu...)

  @. u = uprev + b1*k1 + b2*k2 + b3*k3

  f(t,u,fsallast)

  if integrator.opts.adaptive
    @. utilde = (b1-btilde1)*k1 + (b2-btilde2)*k2 + (b3-btilde3)*k3
    @. atmp = utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

################################################################################

function initialize!(integrator,cache::Rosenbrock34ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Rosenbrock34ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Setup Jacobian Calc
  tf.u = uprev
  uf.t = t

  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I/(dt*gamma) - J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1/(dt*gamma) - J
  end

  d1dt = dt*d1
  linsolve_tmp = @. integrator.fsalfirst + d1dt*dT

  k1 = W\linsolve_tmp
  u = uprev # +a21*k1 a21 == 0
  # du = f(t+c2*dt,u) c2 == 0 and a21 == 0 => du = f(t,uprev) == fsalfirst

  linsolve_tmp = @. integrator.fsalfirst + dt*d2*dT + C21*k1/dt

  k2 = W\linsolve_tmp

  u = @. uprev + a31*k1 + a32*k2

  du = f(t+c3*dt,u)

  linsolve_tmp = @. du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  k3 = W\linsolve_tmp

  # linsolve_tmp = @. du + dt*d4*dT + C41*k1/dt + C42*k2/dt + C43*k3/dt
  linsolve_tmp = @. du + dt*d4*dT + (C41*k1 + C42*k2 + C43*k3)/dt

  k4 = W\linsolve_tmp

  u = @. uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4

  integrator.fsallast = f(t,u)

  if integrator.opts.adaptive
    utilde = @. btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    atmp = @. utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::Rosenbrock34Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [fsalfirst,fsallast]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock34Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  utilde = du
  atmp = du

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  d1dt = dt*d1

  @. linsolve_tmp = fsalfirst + d1dt*dT

  if has_invW(f)
    f(Val{:invW_t},t,uprev,dt*gamma,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
        @inbounds W[i,j] = mass_matrix[i,j]/(dt*gamma)-J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  k1 = reshape(vectmp,sizeu...)

  #=
  a21 == 0 and c2 == 0
  so du = fsalfirst!
  @. u = uprev + a21*k1

  f(t+c2*dt,u,du)
  =#

  @. linsolve_tmp = fsalfirst + dt*d2*dT + C21*k1/dt

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2,sizeu...)

  @. u = uprev + a31*k1 + a32*k2

  f(t+c3*dt,u,du)

  @. linsolve_tmp = du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3,sizeu...)

  # @. linsolve_tmp = du + dt*d4*dT + C41*k1/dt + C42*k2/dt + C43*k3/dt
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + dt*d4*dT[i] + C41*k1[i]/dt + C42*k2[i]/dt + C43*k3[i]/dt
  end

  if has_invW(f)
    A_mul_B!(vectmp4,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4,W,linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4,sizeu...)

  @. u = uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4

  f(t,u,fsallast)

  if integrator.opts.adaptive
    @. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    @. atmp = utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end


################################################################################

#### ROS4 type method

function initialize!(integrator,cache::Rosenbrock4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Rosenbrock4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Setup Jacobian Calc
  tf.u = uprev
  uf.t = t

  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I/(dt*gamma) - J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1/(dt*gamma) - J
  end

  d1dt = dt*d1
  linsolve_tmp = @. integrator.fsalfirst + d1dt*dT

  k1 = W\linsolve_tmp
  u = @. uprev+a21*k1
  du = f(t+c2*dt,u)

  linsolve_tmp = @. du + dt*d2*dT + C21*k1/dt

  k2 = W\linsolve_tmp

  u = @. uprev + a31*k1 + a32*k2

  du = f(t+c3*dt,u)

  linsolve_tmp = @. du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  k3 = W\linsolve_tmp

  # linsolve_tmp = @. du + dt*d4*dT + C41*k1/dt + C42*k2/dt + C43*k3/dt
  linsolve_tmp = @. du + dt*d4*dT + (C41*k1 + C42*k2 + C43*k3)/dt

  k4 = W\linsolve_tmp

  u = @. uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4

  integrator.fsallast = f(t,u)

  if integrator.opts.adaptive
    utilde = @. btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    atmp = @. utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::Rosenbrock4Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [fsalfirst,fsallast]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  utilde = du
  atmp = du

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  d1dt = dt*d1

  @. linsolve_tmp = fsalfirst + d1dt*dT

  if has_invW(f)
    f(Val{:invW_t},t,uprev,dt*gamma,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
        @inbounds W[i,j] = mass_matrix[i,j]/(dt*gamma)-J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  k1 = reshape(vectmp,sizeu...)

  @. u = uprev + a21*k1

  f(t+c2*dt,u,du)

  @. linsolve_tmp = du + dt*d2*dT + C21*k1/dt

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2,sizeu...)

  @. u = uprev + a31*k1 + a32*k2

  f(t+c3*dt,u,du)

  @. linsolve_tmp = du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3,sizeu...)

  # @. linsolve_tmp = du + dt*d4*dT + C41*k1/dt + C42*k2/dt + C43*k3/dt
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + dt*d4*dT[i] + C41*k1[i]/dt + C42*k2[i]/dt + C43*k3[i]/dt
  end

  if has_invW(f)
    A_mul_B!(vectmp4,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4,W,linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4,sizeu...)

  @. u = uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4

  f(t,u,fsallast)

  if integrator.opts.adaptive
    @. utilde = btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    @. atmp = utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

################################################################################

#### Rodas4 type method

function initialize!(integrator,cache::Rodas4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = zero(integrator.fsalfirst)
  integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rodas4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

  # Setup Jacobian Calc
  tf.u = uprev
  uf.t = t

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

  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I/(dt*gamma) - J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1/(dt*gamma) - J
  end

  d1dt = dt*d1

  du = f(t,uprev)

  linsolve_tmp = @. du + d1dt*dT

  k1 = W\linsolve_tmp

  u = @. uprev+a21*k1

  du = f(t+c2*dt,u)

  linsolve_tmp = @. du + dt*d2*dT + dtC21*k1

  k2 = W\linsolve_tmp

  u = @. uprev + a31*k1 + a32*k2

  du = f(t+c3*dt,u)

  linsolve_tmp = @. du + dt*d3*dT + (dtC31*k1 + dtC32*k2)

  k3 = W\linsolve_tmp

  u = @. uprev + a41*k1 + a42*k2 + a43*k3

  du = f(t+c4*dt,u)

  linsolve_tmp = @. du + dt*d4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

  k4 = W\linsolve_tmp

  u = @. uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4

  du = f(t+dt,u)

  linsolve_tmp = @. du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)

  k5 = W\linsolve_tmp

  u = u + k5

  du = f(t+dt,u)

  linsolve_tmp = @. du + (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)

  k6 = W\linsolve_tmp

  u = u + k6

  if integrator.opts.adaptive
    atmp = @. k6/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if integrator.opts.calck
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    integrator.k[1] = @. h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    integrator.k[2] = @. h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
  end

  integrator.fsallast = du
  integrator.u = u
end


function initialize!(integrator,cache::Rodas4Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast,dense1,dense2 = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [dense1,dense2]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rodas4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,vectmp5,vectmp6,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  atmp = du

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

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  d1dt = dt*d1

  f(t,uprev,du)

  @. linsolve_tmp = du + d1dt*dT

  if has_invW(f)
    f(Val{:invW_t},t,uprev,dt*gamma,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
        @inbounds W[i,j] = mass_matrix[i,j]/(dt*gamma)-J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  k1 = reshape(vectmp,sizeu...)

  @. u = uprev+a21*k1

  f(t+c2*dt,u,du)

  @. linsolve_tmp = du + dt*d2*dT + dtC21*k1

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2,sizeu...)

  @. u = uprev + a31*k1 + a32*k2

  f(t+c3*dt,u,du)

  @. linsolve_tmp = du + dt*d3*dT + (dtC31*k1 + dtC32*k2)

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3,sizeu...)

  @. u = uprev + a41*k1 + a42*k2 + a43*k3

  f(t+c4*dt,u,du)

  @. linsolve_tmp = du + dt*d4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

  if has_invW(f)
    A_mul_B!(vectmp4,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4,W,linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4,sizeu...)

  @. u = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4

  f(t+dt,u,du)

  @. linsolve_tmp = du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)

  if has_invW(f)
    A_mul_B!(vectmp5,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp5,W,linsolve_tmp_vec)
  end

  k5 = reshape(vectmp5,sizeu...)

  u .+= k5

  f(t,u,du)

  # @. linsolve_tmp = du + (dtC61*k1 + dtC62*k2 + dtC65*k5 + dtC64*k4 + dtC63*k3)
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + (dtC61*k1[i] + dtC62*k2[i] + dtC65*k5[i] + dtC64*k4[i] + dtC63*k3[i])
  end

  if has_invW(f)
    A_mul_B!(vectmp6,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp6,W,linsolve_tmp_vec)
  end

  k6 = reshape(vectmp6,sizeu...)

  u .+= k6

  if integrator.opts.adaptive
    @. atmp = k6/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
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

function initialize!(integrator,cache::Rosenbrock5ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Rosenbrock5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,C71,C72,C73,C74,C75,C76,C81,C82,C83,C84,C85,C86,C87,gamma,d1,d2,d3,d4,d5,c2,c3,c4,c5 = cache.tab
  # Setup Jacobian Calc
  tf.u = uprev
  uf.t = t

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

  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I/(dt*gamma) - J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1/(dt*gamma) - J
  end

  d1dt = dt*d1

  du1 = f(t,uprev)

  linsolve_tmp = @. du1 + d1dt*dT

  k1 = W\linsolve_tmp

  u = @. uprev+a21*k1

  du = f(t+c2*dt,u)

  linsolve_tmp = @. du + dt*d2*dT + dtC21*k1

  k2 = W\linsolve_tmp

  u = @. uprev + a31*k1 + a32*k2

  du = f(t+c3*dt,u)

  linsolve_tmp = @. du + dt*d3*dT + (dtC31*k1 + dtC32*k2)

  k3 = W\linsolve_tmp

  u = @. uprev + a41*k1 + a42*k2 + a43*k3

  du = f(t+c4*dt,u)

  linsolve_tmp = @. du + dt*d4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

  k4 = W\linsolve_tmp

  u = @. uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4

  du = f(t+c5*dt,u)

  # linsolve_tmp = @. du + dt*d5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    tmp = similar(du)
    @tight_loop_macros for i in eachindex(tmp)
      @inbounds tmp[i] = du[i] + dt*d5*dT[i] + (dtC52*k2[i] + dtC54*k4[i] + dtC51*k1[i] + dtC53*k3[i])
    end
    linsolve_tmp = tmp
  else
    linsolve_tmp = du + dt*d5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  end

  k5 = W\linsolve_tmp

  u = @. uprev + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5

  du = f(t+dt,u)

  linsolve_tmp = @. du + (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)

  k6 = W\linsolve_tmp

  u = u + k6

  du = f(t+dt,u)

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

  du = f(t+dt,u)

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

  du = f(t+dt,u)

  if integrator.opts.adaptive
    atmp = @. k8/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
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

function initialize!(integrator,cache::Rosenbrock5Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast,dense1,dense2 = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [fsalfirst,fsallast]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@muladd function perform_step!(integrator,cache::Rosenbrock5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,vectmp5,vectmp6,vectmp7,vectmp8,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,C71,C72,C73,C74,C75,C76,C81,C82,C83,C84,C85,C86,C87,gamma,d1,d2,d3,d4,d5,c2,c3,c4,c5 = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  atmp = du

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

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,uprev,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  d1dt = dt*d1

  f(t,uprev,fsalfirst)

  @. linsolve_tmp = fsalfirst + d1dt*dT

  if has_invW(f)
    f(Val{:invW_t},t,uprev,dt*gamma,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
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
        @inbounds W[i,j] = mass_matrix[i,j]/(dt*gamma)-J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  k1 = reshape(vectmp,sizeu...)

  @. u = uprev + a21*k1

  f(t+c2*dt,u,du)

  @. linsolve_tmp = du + dt*d2*dT + dtC21*k1

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2,sizeu...)

  @. u = uprev + a31*k1 + a32*k2

  f(t+c3*dt,u,du)

  @. linsolve_tmp = du + dt*d3*dT + (dtC31*k1 + dtC32*k2)

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3,sizeu...)

  @. u = uprev + a41*k1 + a42*k2 + a43*k3

  f(t+c4*dt,u,du)

  @. linsolve_tmp = du + dt*d4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

  if has_invW(f)
    A_mul_B!(vectmp4,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4,W,linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4,sizeu...)

  @. u = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4

  f(t+c5*dt,u,du)

#  @. linsolve_tmp = du + dt*d5*dT + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + dt*d5*dT[i] + (dtC52*k2[i] + dtC54*k4[i] + dtC51*k1[i] + dtC53*k3[i])
  end

  if has_invW(f)
    A_mul_B!(vectmp5,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp5,W,linsolve_tmp_vec)
  end

  k5 = reshape(vectmp5,sizeu...)

  # @. u = uprev + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]
  end

  f(t+dt,u,du)

  # @. linsolve_tmp = du + (dtC61*k1 + dtC62*k2 + dtC63*k3 + dtC64*k4 + dtC65*k5)
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + (dtC61*k1[i] + dtC62*k2[i] + dtC63*k3[i] + dtC64*k4[i] + dtC65*k5[i])
  end

  if has_invW(f)
    A_mul_B!(vectmp6,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp6,W,linsolve_tmp_vec)
  end

  k6 = reshape(vectmp6,sizeu...)

  u .+= k6

  f(t+dt,u,du)

  # @. linsolve_tmp = du + (dtC71*k1 + dtC72*k2 + dtC73*k3 + dtC74*k4 + dtC75*k5 + dtC76*k6)
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + (dtC71*k1[i] + dtC72*k2[i] + dtC73*k3[i] + dtC74*k4[i] + dtC75*k5[i] + dtC76*k6[i])
  end

  if has_invW(f)
    A_mul_B!(vectmp7,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp7,W,linsolve_tmp_vec)
  end

  k7 = reshape(vectmp7,sizeu...)

  u .+= k7

  f(t+dt,u,du)

  # @. linsolve_tmp = du + (dtC81*k1 + dtC82*k2 + dtC83*k3 + dtC84*k4 + dtC85*k5 + dtC86*k6 + dtC87*k7)
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + (dtC81*k1[i] + dtC82*k2[i] + dtC83*k3[i] + dtC84*k4[i] + dtC85*k5[i] + dtC86*k6[i] + dtC87*k7[i])
  end

  if has_invW(f)
    A_mul_B!(vectmp8,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp8,W,linsolve_tmp_vec)
  end

  k8 = reshape(vectmp8,sizeu...)

  u .+= k8

  f(t+dt,u,fsallast)

  if integrator.opts.adaptive
    @. atmp = k8/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
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
