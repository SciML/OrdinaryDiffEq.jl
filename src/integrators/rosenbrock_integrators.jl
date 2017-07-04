@inline function initialize!(integrator,cache::Rosenbrock23Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [k₁,k₂]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@inline function perform_step!(integrator,cache::Rosenbrock23Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  jidx = eachindex(J)
  @unpack c₃₂,d = cache.tab
  mass_matrix = integrator.sol.prob.mass_matrix

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,u,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  γ = dt*d

  #@. linsolve_tmp = @muladd fsalfirst + γ*dT
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = @muladd fsalfirst[i] + γ*dT[i]
  end

  if has_invW(f)
    f(Val{:invW},t,u,γ,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
  else
    if has_jac(f)
      f(Val{:jac},t,u,J)
    else
      if alg_autodiff(integrator.alg)
        ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
      else
        Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
      end
    end
    for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]-γ*J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end


  recursivecopy!(k₁,reshape(vectmp,size(u)...))
  dto2 = dt/2

  #@. u= @muladd uprev+dto2*k₁
  @tight_loop_macros for i in uidx
    @inbounds u[i] = @muladd uprev[i]+dto2*k₁[i]
  end
  f(t+dto2,u,f₁)
  #@. linsolve_tmp = f₁-k₁
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = f₁[i]-k₁[i]
  end

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  tmp2 = reshape(vectmp2,sizeu...)
  #@. k₂ = tmp2 + k₁
  #@. u = @muladd uprev + dt*k₂
  @tight_loop_macros for i in uidx
    @inbounds k₂[i] = tmp2[i] + k₁[i]
    @inbounds u[i] = @muladd uprev[i] + dt*k₂[i]
  end

  if integrator.opts.adaptive
    f(t+dt,u,integrator.fsallast)

    #@. linsolve_tmp = @muladd integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT
    @tight_loop_macros for i in uidx
      @inbounds linsolve_tmp[i] = @muladd fsallast[i] - c₃₂*(k₂[i]-f₁[i])-2(k₁[i]-fsalfirst[i])+dt*dT[i]
    end

    if has_invW(f)
      A_mul_B!(vectmp3,W,linsolve_tmp_vec)
    else
      integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
    end

    k₃ = reshape(vectmp3,sizeu...)
    #@. tmp = (dt*(k₁ - 2k₂ + k₃)/6)./@muladd(integrator.opts.abstol+max(abs(uprev),abs(u)).*integrator.opts.reltol)
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds tmp[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./@muladd(atol+max(abs(uprev[i]),abs(u[i])).*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::Rosenbrock32Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [k₁,k₂]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end#

@inline function perform_step!(integrator,cache::Rosenbrock32Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  jidx = eachindex(J)
  mass_matrix = integrator.sol.prob.mass_matrix
  @unpack c₃₂,d = cache.tab
  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  if has_tgrad(f)
    f(Val{:tgrad},t,u,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  γ = dt*d

  #@. linsolve_tmp = @muladd fsalfirst + γ*dT
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = @muladd fsalfirst[i] + γ*dT[i]
  end

  if has_invW(f)
    f(Val{:invW},t,u,γ,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
  else
    if has_jac(f)
      f(Val{:jac},t,u,J)
    else
      if alg_autodiff(integrator.alg)
        ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
      else
        Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
      end
    end
    for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]-γ*J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  recursivecopy!(k₁,reshape(vectmp,sizeu...))

  dto2 = dt/2
  #@. u= @muladd uprev+dto2*k₁
  @tight_loop_macros for i in uidx
    @inbounds u[i] = @muladd uprev[i]+dto2*k₁[i]
  end
  f(t+dto2,u,f₁)
  #@. linsolve_tmp = f₁-k₁
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = f₁[i]-k₁[i]
  end

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  tmp2 = reshape(vectmp2,sizeu...)
  #@. k₂ = tmp2 + k₁
  #@. u = @muladd uprev + dt*k₂
  @tight_loop_macros for i in uidx
    @inbounds k₂[i] = tmp2[i] + k₁[i]
    @inbounds tmp[i] = @muladd uprev[i] + dt*k₂[i]
  end


  f(t+dt,tmp,integrator.fsallast)
  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = @muladd fsallast[i] - c₃₂*(k₂[i]-f₁[i])-2(k₁[i]-fsalfirst[i])+dt*dT[i]
  end

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k₃ = reshape(vectmp3,sizeu...)
  #@. u = uprev + dt*(k₁ + 4k₂ + k₃)/6
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
  end

  if integrator.opts.adaptive
    #@. tmp = (dt*(k₁ - 2k₂ + k₃)/6)./@muladd(integrator.opts.abstol+max(abs(uprev),abs(u)).*integrator.opts.reltol)
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds tmp[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./@muladd(atol+max(abs(uprev[i]),abs(u[i])).*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::Rosenbrock23ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)
end

@inline function perform_step!(integrator,cache::Rosenbrock23ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack c₃₂,d,tf,uf = cache
  # Time derivative
  tf.u = uprev
  uf.t = t
  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
  else
    J = ForwardDiff.derivative(uf,uprev)
  end
  a = -dt*d
  W = @. @muladd 1+a*J
  a = -a
  k₁ = W\@.(@muladd(integrator.fsalfirst + a*dT))
  dto2 = dt/2
  f₁ = f(t+dto2,@.(@muladd uprev+dto2*k₁))
  k₂ = W\(f₁-k₁) + k₁
  u = @. @muladd uprev + dt*k₂
  if integrator.opts.adaptive
    integrator.fsallast = f(t+dt,u)
    k₃ = W\@.(@muladd(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT))
    tmp = @. (dt*(k₁ - 2k₂ + k₃)/6)./@muladd(integrator.opts.abstol+max.(abs.(uprev),abs.(u)).*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::Rosenbrock32ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)
end

@inline function perform_step!(integrator,cache::Rosenbrock32ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack c₃₂,d,tf,uf = cache
  tf.u = uprev
  uf.t = t
  # Time derivative
  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
  else
    J = ForwardDiff.derivative(uf,uprev)
  end
  a = -dt*d
  W = @. 1+a*J
  #f₀ = f(t,uprev)
  a = -a
  k₁ = W\@.(@muladd(integrator.fsalfirst + a*dT))
  dto2 = dt/2
  f₁ = f(t+dto2,@.(@muladd(uprev+dto2*k₁)))
  k₂ = W\(f₁-k₁) + k₁
  tmp = @. @muladd uprev + dt*k₂
  integrator.fsallast = f(t+dt,tmp)
  k₃ = W\@.(@muladd(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT))
  u = uprev + dt*(k₁ + 4k₂ + k₃)/6
  if integrator.opts.adaptive
    tmp = @. (dt*(k₁ - 2k₂ + k₃)/6)./@muladd(integrator.opts.abstol+max.(abs.(uprev),abs.(u)).*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::Rosenbrock4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)
end

@inline function perform_step!(integrator,cache::Rosenbrock4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack tf,uf = cache
  @unpack a21,a31,a32,C21,C31,C32,C41,C42,C43,b1,b2,b3,b4,btilde1,btilde2,btilde3,btilde4,gamma,c2,c3,d1,d2,d3,d4 = cache.tab

  # Setup Jacobian Calc
  tf.u = uprev
  uf.t = t

  dT = ForwardDiff.derivative(tf,t)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
  else
    J = ForwardDiff.derivative(uf,uprev)
  end

  d1dt = dt*d1
  linsolve_tmp = integrator.fsalfirst + d1dt*dT
  W = @muladd 1/(dt*gamma)-J

  k1 = W\linsolve_tmp
  u = uprev+a21*k1
  du = f(t+c2*dt,u)

  linsolve_tmp = du + dt*d2*dT + C21*k1/dt

  k2 = W\linsolve_tmp

  u = uprev + a31*k1 + a32*k2

  du = f(t+c3*dt,u)

  linsolve_tmp = du + dt*d3*dT + C31*k1/dt + C32*k2/dt

  k3 = W\linsolve_tmp

  linsolve_tmp = du + dt*d4*dT + C41*k1/dt + C42*k2/dt + C43*k3/dt

  k4 = W\linsolve_tmp

  u = uprev + b1*k1 + b2*k2 + b3*k3 + b4*k4

  fsallast = f(t,u)

  if integrator.opts.adaptive
    utilde = @muladd btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4
    atmp = ((utilde)./@muladd(integrator.opts.abstol+max(abs.(uprev),abs.(u)).*integrator.opts.reltol))
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = fsallast
  integrator.fsallast = fsallast
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::Rosenbrock4Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [fsalfirst,fsallast]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@inline function perform_step!(integrator,cache::Rosenbrock4Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,fsalfirst,fsallast,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config = cache
  jidx = eachindex(J)
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
    f(Val{:tgrad},t,u,dT)
  else
    if alg_autodiff(integrator.alg)
      ForwardDiff.derivative!(dT,tf,vec(du2),t) # Should update to inplace, https://github.com/JuliaDiff/ForwardDiff.jl/pull/219
    else
      dT = Calculus.finite_difference(tf,t,integrator.alg.diff_type)
    end
  end

  d1dt = dt*d1

  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = fsalfirst[i] + d1dt*dT[i]
    #@inbounds linsolve_tmp[i] = du[i] + d1dt*dT[i]
  end

  if has_invW(f)
    f(Val{:invW_t},t,u,dt*gamma,W) # W == inverse W
    A_mul_B!(vectmp,W,linsolve_tmp_vec)
  else
    if has_jac(f)
      f(Val{:jac},t,u,J)
    else
      if alg_autodiff(integrator.alg)
        ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
      else
        Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
      end
    end
    for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]/(dt*gamma)-J[i,j]
    end
    integrator.alg.linsolve(vectmp,W,linsolve_tmp_vec,true)
  end

  recursivecopy!(k1,reshape(vectmp,size(u)...))

  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+a21*k1[i]
  end

  f(t+c2*dt,u,du)

  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + dt*d2*dT[i] + C21*k1[i]/dt
  end

  if has_invW(f)
    A_mul_B!(vectmp2,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp2,W,linsolve_tmp_vec)
  end

  k2 = reshape(vectmp2,sizeu...)

  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + a31*k1[i] + a32*k2[i]
  end

  f(t+c3*dt,u,du)

  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + dt*d3*dT[i] + C31*k1[i]/dt + C32*k2[i]/dt
  end

  if has_invW(f)
    A_mul_B!(vectmp3,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp3,W,linsolve_tmp_vec)
  end

  k3 = reshape(vectmp3,sizeu...)

  @tight_loop_macros for i in uidx
    @inbounds linsolve_tmp[i] = du[i] + dt*d4*dT[i] + C41*k1[i]/dt + C42*k2[i]/dt + C43*k3[i]/dt
  end

  if has_invW(f)
    A_mul_B!(vectmp4,W,linsolve_tmp_vec)
  else
    integrator.alg.linsolve(vectmp4,W,linsolve_tmp_vec)
  end

  k4 = reshape(vectmp4,sizeu...)

  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i]
  end

  f(t,u,fsallast)

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds utilde[i] = @muladd btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i]
      @inbounds atmp[i] = ((utilde[i])./@muladd(atol+max(abs(uprev[i]),abs(u[i])).*rtol))
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @pack integrator = t,dt,u,k
end
