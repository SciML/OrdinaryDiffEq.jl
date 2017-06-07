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
    @fastmath for i in 1:length(u)
      @simd for j in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]-γ*J[i,j]
      end
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
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./@muladd(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i])).*integrator.opts.reltol)
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
    @fastmath for i in 1:length(u)
      @simd for j in 1:length(u)
        @inbounds W[i,j] = @muladd mass_matrix[i,j]-γ*J[i,j]
      end
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
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./@muladd(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i])).*integrator.opts.reltol)
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
