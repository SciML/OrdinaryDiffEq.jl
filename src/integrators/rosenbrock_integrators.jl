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
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,uf,tf,linsolve_tmp = cache
  jidx = eachindex(J)
  @unpack c₃₂,d = cache.tab

  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  #if alg_autodiff(alg)
    ForwardDiff.derivative!(dT,tf,t) # Time derivative of each component
    ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev))
  #else
  #  Calculus.finite_difference!((t)->vecfreturn(t,uprev,du2),[t],dT)
  #  Calculus.finite_difference_jacobian!((du1,uprev)->vecf(t,uprev,du1),vec(uprev),vec(du1),J)
  #end
  a = -dt*d
  for i in 1:length(u), j in 1:length(u)
    W[i,j] = @muladd I[i,j]+a*J[i,j]
  end

  Wfact = integrator.alg.factorization(W)
  a = -a
  for i in uidx
    linsolve_tmp[i] = @muladd fsalfirst[i] + a*dT[i]
  end

  @into! vectmp = Wfact\linsolve_tmp
  recursivecopy!(k₁,reshape(vectmp,size(u)...))
  dto2 = dt/2
  for i in uidx
    u[i]= @muladd uprev[i]+dto2*k₁[i]
  end
  f(t+dto2,u,f₁)

  for i in uidx
    linsolve_tmp[i] = f₁[i]-k₁[i]
  end

  @into! vectmp2 = Wfact\linsolve_tmp
  for i in uidx
    k₂[i] = tmp[i] + k₁[i]
    u[i] = @muladd uprev[i] + dt*k₂[i]
  end
  if integrator.opts.adaptive
    f(t+dt,u,integrator.fsallast)

    for i in uidx
      linsolve_tmp[i] = integrator.fsallast[i] - c₃₂*(k₂[i]-f₁[i])-2(k₁[i]-fsalfirst[i])+dt*dT[i]
    end

    @into! vectmp3 = Wfact\linsolve_tmp
    k₃ = reshape(vectmp3,sizeu...)
    for i in uidx
      tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./@muladd(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp2)
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
end

@inline function perform_step!(integrator,cache::Rosenbrock32Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,uf,tf,linsolve_tmp = cache
  jidx = eachindex(J)
  @unpack c₃₂,d = cache.tab
  # Setup Jacobian Calc
  sizeu  = size(u)
  tf.vf.sizeu = sizeu
  tf.uprev = uprev
  uf.vfr.sizeu = sizeu
  uf.t = t

  #if alg_autodiff(alg)
    ForwardDiff.derivative!(dT,tf,t) # Time derivative of each component
    ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev))

  a = -dt*d
  for i in 1:length(u), j in 1:length(u)
    W[i,j] = @muladd I[i,j]+a*J[i,j]
  end


  Wfact = integrator.alg.factorization(W)

  a = -a
  for i in uidx
    linsolve_tmp[i] = @muladd fsalfirst[i] + a*dT[i]
  end

  @into! vectmp = Wfact\linsolve_tmp

  recursivecopy!(k₁,reshape(vectmp,sizeu...))

  dto2 = dt/2
  for i in uidx
    u[i]= @muladd uprev[i]+dto2*k₁[i]
  end
  f(t+dto2,u,f₁)
  for i in uidx
    linsolve_tmp[i] = f₁[i]-k₁[i]
  end
  @into! vectmp2 = Wfact\linsolve_tmp
  tmp = reshape(vectmp2,sizeu...)
  for i in uidx
    k₂[i] = tmp[i] + k₁[i]
  end
  for i in uidx
    tmp[i] = @muladd uprev[i] + dt*k₂[i]
  end
  f(t+dt,tmp,integrator.fsallast)
  for i in uidx
    linsolve_tmp[i] = integrator.fsallast[i] - c₃₂*(k₂[i]-f₁[i])-2(k₁[i]-fsalfirst[i])+dt*dT[i]
  end
  @into! vectmp3 = Wfact\linsolve_tmp
  k₃ = reshape(vectmp3,sizeu...)
  for i in uidx
    u[i] = uprev[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
  end
  if integrator.opts.adaptive
    for i in uidx
      tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/@muladd(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp2)
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
  J = ForwardDiff.derivative(uf,uprev)
  a = -dt*d
  W = @muladd 1+a*J
  a = -a
  k₁ = W\@muladd(integrator.fsalfirst + a*dT)
  dto2 = dt/2
  f₁ = f(t+dto2,@muladd uprev+dto2*k₁)
  k₂ = W\(f₁-k₁) + k₁
  u = @muladd uprev + dt*k₂
  if integrator.opts.adaptive
    integrator.fsallast = f(t+dt,u)
    k₃ = W\(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
    integrator.EEst = integrator.opts.internalnorm((dt*(k₁ - 2k₂ + k₃)/6)./@muladd(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
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
  J = ForwardDiff.derivative(uf,uprev)
  a = -dt*d
  W = 1+a*J
  #f₀ = f(t,uprev)
  a = -a
  k₁ = W\@muladd(integrator.fsalfirst + a*dT)
  dto2 = dt/2
  f₁ = f(t+dto2,@muladd(uprev+dto2*k₁))
  k₂ = W\(f₁-k₁) + k₁
  tmp = @muladd uprev + dt*k₂
  integrator.fsallast = f(t+dt,tmp)
  k₃ = W\(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
  u = uprev + dt*(k₁ + 4k₂ + k₃)/6
  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm((dt*(k₁ - 2k₂ + k₃)/6)./@muladd(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  @pack integrator = t,dt,u,k
end
