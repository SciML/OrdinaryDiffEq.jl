function initialize!(integrator, cache::LawsonEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Pre-start fsal
  lin = integrator.f.f1(integrator.uprev,integrator.p,integrator.t)
  nl = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  integrator.fsalfirst = ExpRKFsal(lin, nl)

  # Avoid undefined entries if k is an array of arrays
  rate_prototype = lin
  integrator.fsallast = ExpRKFsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::LawsonEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack lin,nl = integrator.fsalfirst
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  integrator.k[1] = lin + nl

  if alg.krylov
    @muladd u = _expmv(dt, f.f1, uprev + dt*nl; m=min(alg.m, size(f.f1,1)), norm=normbound)
  else
    @muladd u = cache.exphA*(uprev + dt*nl)
  end

  # Push the fsal at t+dt
  lin = f.f1(u,p,t+dt)
  nl = f.f2(u,p,t+dt)
  integrator.k[2] = lin + nl
  @pack integrator.fsallast = lin, nl
  integrator.u = u
end

function initialize!(integrator, cache::LawsonEulerCache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  rate_prototype = cache.rtmp

  # Pre-start fsal
  integrator.fsalfirst = ExpRKFsal(rate_prototype)
  @unpack lin,nl = integrator.fsalfirst
  integrator.f.f1(lin,integrator.uprev,integrator.p,integrator.t)
  integrator.f.f2(nl,integrator.uprev,integrator.p,integrator.t)

  integrator.fsallast = ExpRKFsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::LawsonEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack lin,nl = integrator.fsalfirst
  @unpack tmp,exphA,Ks,KsCache = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  @. integrator.k[1] = lin + nl

  @muladd @. tmp = uprev + dt*nl
  if alg.krylov
    arnoldi!(Ks,f.f1,tmp; m=min(alg.m, size(f.f1,1)), norm=normbound, cache=u)
    _expmv!(u,dt,Ks; cache=KsCache)
  else
    A_mul_B!(u,exphA,tmp)
  end

  # Push the fsal at t+dt
  @unpack lin,nl = integrator.fsallast
  f.f1(lin,u,p,t+dt)
  f.f2(nl,u,p,t+dt)
  @. integrator.k[2] = lin + nl
end

function initialize!(integrator, cache::NorsettEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Pre-start fsal
  lin = integrator.f.f1(integrator.uprev,integrator.p,integrator.t)
  nl = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  integrator.fsalfirst = ExpRKFsal(lin, nl)

  # Avoid undefined entries if k is an array of arrays
  rate_prototype = lin
  integrator.fsallast = ExpRKFsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::NorsettEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack lin,nl = integrator.fsalfirst
  @unpack exphA,phihA = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  integrator.k[1] = lin + nl

  if alg.krylov
    u = phimv(dt,f.f1,nl,uprev; tol=integrator.opts.reltol, m=min(alg.m, size(f.f1,1)), norm=normbound)
  else
    u = exphA*uprev + dt*(phihA*nl)
  end

  # Push the fsal at t+dt
  lin = f.f1(u,p,t+dt)
  nl = f.f2(u,p,t+dt)
  integrator.k[2] = lin + nl
  @pack integrator.fsallast = lin, nl
  integrator.u = u
end

function initialize!(integrator, cache::NorsettEulerCache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  rate_prototype = cache.rtmp

  # Pre-start fsal
  integrator.fsalfirst = ExpRKFsal(rate_prototype)
  @unpack lin,nl = integrator.fsalfirst
  integrator.f.f1(lin,integrator.uprev,integrator.p,integrator.t)
  integrator.f.f2(nl,integrator.uprev,integrator.p,integrator.t)

  integrator.fsallast = ExpRKFsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::NorsettEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack lin,nl = integrator.fsalfirst
  @unpack tmp,rtmp,exphA,phihA = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  @. integrator.k[1] = lin + nl

  if alg.krylov
    phimv!(u,dt,f.f1,nl,uprev; tol=integrator.opts.reltol, m=min(alg.m, size(f.f1,1)), norm=normbound)
  else
    A_mul_B!(tmp,exphA,uprev)
    A_mul_B!(rtmp,phihA,nl)
    @muladd @. u = tmp + dt*rtmp
  end

  # Push the fsal at t+dt
  @unpack lin,nl = integrator.fsallast
  f.f1(lin,u,p,t+dt)
  f.f2(nl,u,p,t+dt)
  @. integrator.k[2] = lin + nl
end

function initialize!(integrator,cache::ETD2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Pre-start fsal
  lin = integrator.f.f1(integrator.uprev,integrator.p,integrator.t)
  nl = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  nlprev = zero(nl) # to be computed in the first iteration via ETD1
  integrator.fsalfirst = ETD2Fsal(lin, nl, nlprev)
    
  # Avoid undefined entries if k is an array of arrays
  rate_prototype = lin
  integrator.fsallast = ETD2Fsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator,cache::ETD2ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack lin,nl,nlprev = integrator.fsalfirst
  @unpack exphA,phihA,B1,B0 = cache
  integrator.k[1] = lin + nl

  if integrator.iter == 1 # ETD1 for initial step
    u = exphA*uprev + dt*(phihA*nl)
  else
    u = exphA*uprev + dt*(B1*nl + B0*nlprev)
  end
  integrator.u = u

  # Push the fsal at t+dt
  nlprev = nl
  lin = f.f1(u,p,t+dt)
  nl = f.f2(u,p,t+dt)
  integrator.k[2] = lin + nl
  @pack integrator.fsallast = lin, nl, nlprev
end

function initialize!(integrator, cache::ETD2Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  rate_prototype = cache.rtmp1
  
  # Pre-start fsal
  integrator.fsalfirst = ETD2Fsal(rate_prototype)
  @unpack lin,nl = integrator.fsalfirst
  integrator.f.f1(lin,integrator.uprev,integrator.p,integrator.t)
  integrator.f.f2(nl,integrator.uprev,integrator.p,integrator.t)
    
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = ETD2Fsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::ETD2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack lin,nl,nlprev = integrator.fsalfirst
  @unpack utmp,rtmp1,rtmp2,exphA,phihA,B1,B0 = cache
  @. integrator.k[1] = lin + nl

  if integrator.iter == 1 # ETD1 for initial step
    A_mul_B!(utmp, exphA, uprev)
    A_mul_B!(rtmp1, phihA, nl)
    @muladd @. u = utmp + dt*rtmp1
  else
    A_mul_B!(utmp, exphA, uprev)
    A_mul_B!(rtmp1, B1, nl)
    A_mul_B!(rtmp2, B0, nlprev)
    @muladd @. u = utmp + dt*(rtmp1 + rtmp2)
  end

  # Push the fsal at t+dt
  fsallast = integrator.fsallast
  fsallast.nlprev .= nl
  f.f1(fsallast.lin,u,p,t+dt)
  f.f2(fsallast.nl,u,p,t+dt)
  @. integrator.k[2] = fsallast.lin + fsallast.nl
end

function initialize!(integrator, cache::ETDRK4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Pre-start fsal
  lin = integrator.f.f1(integrator.uprev,integrator.p,integrator.t)
  nl = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  integrator.fsalfirst = ExpRKFsal(lin, nl)
    
  # Avoid undefined entries if k is an array of arrays
  rate_prototype = lin
  integrator.fsallast = ExpRKFsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::ETDRK4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack lin,nl = integrator.fsalfirst
  @unpack E,E2,a,b,c,Q = cache
  integrator.k[1] = lin + nl

  tmp = E2*uprev
  k1 = nl # k1 is fsaled
  s1 = tmp + Q*k1;
  k2 = integrator.f.f2(s1,p,t+dt/2)
  s2 = tmp + Q*k2;
  k3 = integrator.f.f2(s2,p,t+dt/2)
  s3 = E2*s1 + Q*(2*k3-k1);
  k4 = integrator.f.f2(s3,p,t+dt)
  u = E*uprev + a*k1 + 2b*(k2+k3) + c*k4;

  # Push the fsal at t+dt
  lin = f.f1(u,p,t+dt)
  nl = f.f2(u,p,t+dt)
  integrator.k[2] = lin + nl
  @pack integrator.fsallast = lin, nl
  integrator.u = u
end

function initialize!(integrator, cache::ETDRK4Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  rate_prototype = cache.tmp2

  # Pre-start fsal
  integrator.fsalfirst = ExpRKFsal(rate_prototype)
  @unpack lin,nl = integrator.fsalfirst
  integrator.f.f1(lin,integrator.uprev,integrator.p,integrator.t)
  integrator.f.f2(nl,integrator.uprev,integrator.p,integrator.t)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = ExpRKFsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::ETDRK4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack lin,nl = integrator.fsalfirst
  @unpack tmp2,tmp = cache
  @unpack E,E2,a,b,c,Q = cache
  @unpack k2,k3,k4,s1 = cache
  @. integrator.k[1] = lin + nl

  # Substep 1
  k1 = nl # k1 is fsaled
  A_mul_B!(tmp,E2,uprev)
  A_mul_B!(tmp2,Q,k1)
  @. s1 = tmp + tmp2

  # Substep 2
  integrator.f.f2(k2,s1,p,t+dt/2)
  A_mul_B!(tmp2,Q,k2)
  # tmp is still E2*uprev
  @. tmp2 = tmp + tmp2

  # Substep 3
  integrator.f.f2(k3,tmp2,p,t+dt/2)
  @. tmp = 2.0*k3 - k1
  A_mul_B!(tmp2,Q,tmp)
  A_mul_B!(tmp,E2,s1)
  @. tmp2 = tmp + tmp2

  # Substep 4
  integrator.f.f2(k4,tmp2,p,t+dt)

  # Update
  @. tmp2 = k2+k3
  A_mul_B!(tmp,b,tmp2)
  A_mul_B!(s1,E,uprev)
  A_mul_B!(k2,a,k1)
  A_mul_B!(k3,c,k4)
  @. u = s1 + k2 + 2tmp + k3

  # Push the fsal at t+dt
  @unpack lin,nl = integrator.fsallast
  f.f1(lin,u,p,t+dt)
  f.f2(nl,u,p,t+dt)
  @. integrator.k[2] = lin + nl
end
