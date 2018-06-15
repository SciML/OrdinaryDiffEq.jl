# Helper function to compute the G_nj factors for the classical ExpRK methods
@inline _compute_nl(f::SplitFunction, u, p, t, A) = f.f2(u, p, t)
@inline _compute_nl(f::DiffEqFunction, u, p, t, A) = f(u, p, t) - A * u
@inline _compute_nl!(G, f::SplitFunction, u, p, t, A, Au_cache) = f.f2(G, u, p, t)
@inline function _compute_nl!(G, f::DiffEqFunction, u, p, t, A, Au_cache)
  f(G, u, p, t)
  A_mul_B!(Au_cache, A, u)
  G .-= Au_cache
end

##########################################
# Common initializers for ExpRK integrators
function initialize!(integrator, cache::ExpRKConstantCache)
  # Pre-start fsal
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.fsallast = zero(integrator.fsalfirst)

  # Initialize interpolation derivatives
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end
function initialize!(integrator, cache::ExpRKCache)
  # Pre-start fsal
  integrator.fsalfirst = zero(cache.rtmp)
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.fsallast = zero(integrator.fsalfirst)

  # Initialize interpolation derivatives
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

###########################################
function perform_step!(integrator, cache::LawsonEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  nl = _compute_nl(f, uprev, p, t, A)
  @muladd v = uprev + dt * nl
  if alg.krylov
    u = expv(dt, A, v; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
  else
    u = cache.exphA * v
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::LawsonEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,G,Jcache,exphA,Ks,KsCache = cache
  if isa(f, SplitFunction)
    A = f.f1
  else
    f.jac(Jcache, uprev, p, t)
    A = Jcache
  end
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  _compute_nl!(G, f, uprev, p, t, A, rtmp)
  @muladd @. tmp = uprev + dt*G
  if alg.krylov
    arnoldi!(Ks, f.f1, tmp; m=min(alg.m, size(f.f1,1)), norm=integrator.opts.internalnorm, 
      cache=u, iop=alg.iop)
    expv!(u,dt,Ks; cache=KsCache)
  else
    A_mul_B!(u,exphA,tmp)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::NorsettEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    w = phiv(dt, A, integrator.fsalfirst, 1; m=min(alg.m, size(A,1)), 
      norm=integrator.opts.internalnorm, iop=alg.iop)
    u = uprev + dt * w[:,2]
  else
    u = uprev + dt * (cache.phihA * integrator.fsalfirst)
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::NorsettEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack rtmp,Jcache,Ks,KsCache = cache
  if isa(f, SplitFunction)
    A = f.f1
  else
    f.jac(Jcache, uprev, p, t)
    A = Jcache
  end
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    w = KsCache[1]
    arnoldi!(Ks, A, integrator.fsalfirst; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, 
      cache=u, iop=alg.iop)
    phiv!(w, dt, Ks, 1; caches=KsCache[2:end])
    @muladd @. u = uprev + dt * @view(w[:, 2])
  else
    A_mul_B!(rtmp, cache.phihA, integrator.fsalfirst)
    @muladd @. u = uprev + dt*rtmp
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
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

function perform_step!(integrator, cache::ETDRK4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack E,E2,a,b,c,Q = cache
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator

  tmp = E2*uprev
  k1 = _compute_nl(f, uprev, p, t, A)
  s1 = tmp + Q*k1;
  k2 = _compute_nl(f, s1, p, t + dt/2, A)
  s2 = tmp + Q*k2;
  k3 = _compute_nl(f, s2, p, t + dt/2, A)
  s3 = E2*s1 + Q*(2*k3-k1);
  k4 = _compute_nl(f, s3, p, t + dt, A)
  u = E*uprev + a*k1 + 2b*(k2+k3) + c*k4;

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp2,tmp,rtmp,Jcache = cache
  @unpack E,E2,a,b,c,Q = cache
  @unpack k1,k2,k3,k4,s1 = cache
  if isa(f, SplitFunction)
    A = f.f1
  else
    f.jac(Jcache, uprev, p, t)
    A = Jcache
  end

  # Substep 1
  _compute_nl!(k1, f, uprev, p, t, A, rtmp)
  A_mul_B!(tmp,E2,uprev)
  A_mul_B!(tmp2,Q,k1)
  @. s1 = tmp + tmp2

  # Substep 2
  _compute_nl!(k2, f, s1, p, t + dt/2, A, rtmp)
  A_mul_B!(tmp2,Q,k2)
  # tmp is still E2*uprev
  @. tmp2 = tmp + tmp2

  # Substep 3
  _compute_nl!(k3, f, tmp2, p, t + dt/2, A, rtmp)
  @. tmp = 2.0*k3 - k1
  A_mul_B!(tmp2,Q,tmp)
  A_mul_B!(tmp,E2,s1)
  @. tmp2 = tmp + tmp2

  # Substep 4
  _compute_nl!(k4, f, tmp2, p, t + dt, A, rtmp)

  # Update
  @. tmp2 = k2+k3
  A_mul_B!(tmp,b,tmp2)
  A_mul_B!(s1,E,uprev)
  A_mul_B!(k2,a,k1)
  A_mul_B!(k3,c,k4)
  @. u = s1 + k2 + 2tmp + k3

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end
