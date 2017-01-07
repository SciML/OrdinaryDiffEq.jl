@inline function initialize!(integrator::ODEIntegrator,cache::Rosenbrock23Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [k₁,k₂]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::Rosenbrock23Cache)
  @unpack t,dt,uprev,u,f,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,uf,tf = cache
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

  W[:] = I-dt*d*J # Can an allocation be cut here?
  @into! vectmp = W\vec(fsalfirst + dt*d*dT)
  recursivecopy!(k₁,reshape(vectmp,size(u)...))
  for i in uidx
    u[i]=uprev[i]+dt*k₁[i]/2
  end
  f(t+dt/2,u,f₁)
  @into! vectmp2 = W\vec(f₁-k₁)
  for i in uidx
    k₂[i] = tmp[i] + k₁[i]
    u[i] = uprev[i] + dt*k₂[i]
  end
  if integrator.opts.adaptive
    f(t+dt,u,integrator.fsallast)
    @into! vectmp3 = W\vec(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
    k₃ = reshape(vectmp3,sizeu...)
    for i in uidx
      tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp2)
  end
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator::ODEIntegrator,cache::Rosenbrock32Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k₁,k₂,fsalfirst,fsallast = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  integrator.k = [k₁,k₂]
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::Rosenbrock32Cache)
  @unpack t,dt,uprev,u,f,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,uf,tf = cache
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

  W[:] = I-dt*d*J # Can an allocation be cut here?
  @into! vectmp = W\vec(integrator.fsalfirst + dt*d*dT)
  recursivecopy!(k₁,reshape(vectmp,sizeu...))
  for i in uidx
    u[i]=uprev[i]+dt*k₁[i]/2
  end
  f(t+dt/2,u,f₁)
  @into! vectmp2 = W\vec(f₁-k₁)
  tmp = reshape(vectmp2,sizeu...)
  for i in uidx
    k₂[i] = tmp[i] + k₁[i]
  end
  for i in uidx
    tmp[i] = uprev[i] + dt*k₂[i]
  end
  f(t+dt,tmp,integrator.fsallast)
  @into! vectmp3 = W\vec(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
  k₃ = reshape(vectmp3,sizeu...)
  for i in uidx
    u[i] = uprev[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
  end
  if integrator.opts.adaptive
    for i in uidx
      tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp2)
  end
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator::ODEIntegrator,cache::Rosenbrock23ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::Rosenbrock23ConstantCache)
  @unpack t,dt,uprev,u,f,k = integrator
  @unpack c₃₂,d,tf,uf = cache
  # Time derivative
  tf.u = uprev
  uf.t = t
  dT = ForwardDiff.derivative(tf,t)
  J = ForwardDiff.derivative(uf,uprev)
  W = 1-dt*d*J
  k₁ = W\(integrator.fsalfirst + dt*d*dT)
  f₁ = f(t+dt/2,uprev+dt*k₁/2)
  k₂ = W\(f₁-k₁) + k₁
  u = uprev + dt*k₂
  if integrator.opts.adaptive
    integrator.fsallast = f(t+dt,u)
    k₃ = W\(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
    integrator.EEst = integrator.opts.internalnorm((dt*(k₁ - 2k₂ + k₃)/6)./(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator::ODEIntegrator,cache::Rosenbrock32ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  k = eltype(integrator.sol.k)(2)
  integrator.k = k
  integrator.fsalfirst = f(integrator.t,integrator.uprev)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::Rosenbrock32ConstantCache)
  @unpack t,dt,uprev,u,f,k = integrator
  @unpack c₃₂,d,tf,uf = cache
  tf.u = uprev
  uf.t = t
  # Time derivative
  dT = ForwardDiff.derivative(tf,t)
  J = ForwardDiff.derivative(uf,uprev)
  W = 1-dt*d*J
  #f₀ = f(t,uprev)
  k₁ = W\(integrator.fsalfirst + dt*d*dT)
  f₁ = f(t+dt/2,uprev+dt*k₁/2)
  k₂ = W\(f₁-k₁) + k₁
  tmp = uprev + dt*k₂
  integrator.fsallast = f(t+dt,tmp)
  k₃ = W\(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
  u = uprev + dt*(k₁ + 4k₂ + k₃)/6
  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm((dt*(k₁ - 2k₂ + k₃)/6)./(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.k[1] = k₁
  integrator.k[2] = k₂
  @pack integrator = t,dt,u,k
end
