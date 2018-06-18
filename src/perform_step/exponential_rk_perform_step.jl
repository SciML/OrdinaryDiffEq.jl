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
# Classical ExpRK integrators
function perform_step!(integrator, cache::LawsonEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  nl = _compute_nl(f, uprev, p, t, A)
  @muladd v = uprev + dt * nl
  if alg.krylov
    u = expv(dt, A, v; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
  else
    exphA = cache.ops
    u = exphA * v
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
    phihA = cache.ops
    u = uprev + dt * (phihA * integrator.fsalfirst)
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

function perform_step!(integrator, cache::ETDRK2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    F1 = integrator.fsalfirst
    w1 = phiv(dt, A, F1, 2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U2 = uprev + dt * w1[:, 2]
    F2 = _compute_nl(f, U2, p, t + dt, A) + A * uprev
    w2 = phiv(dt, A, F2, 2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    u = uprev + dt * (w1[:, 2] - w1[:, 3] + w2[:, 3])
  else
    phi1, phi2 = cache.ops
    # The caching version uses a special formula to save computation
    G1 = f.f2(uprev, p, t)
    F1 = integrator.fsalfirst
    U2 = uprev + dt * (phi1 * F1)
    G2 = f.f2(U2, p, t + dt)
    u = U2 + dt * (phi2 * (G2 - G1))
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,F2,Jcache,Ks,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    F1 = integrator.fsalfirst
    w1, w2, phiv_caches = KsCache
    # Krylov for F1
    arnoldi!(Ks, A, F1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w1, dt, Ks, 2; caches=phiv_caches)
    # Krylov for F2
    @muladd @. tmp = uprev + dt * @view(w1[:, 2])
    _compute_nl!(F2, f, tmp, p, t + dt, A, rtmp)
    F2 .+= A_mul_B!(rtmp, A, uprev)
    arnoldi!(Ks, A, F2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w2, dt, Ks, 2; caches=phiv_caches)
    # Update u
    u .= uprev
    Base.axpy!( dt, @view(w1[:, 2]), u)
    Base.axpy!(-dt, @view(w1[:, 3]), u)
    Base.axpy!( dt, @view(w2[:, 3]), u)
  else
    phi1, phi2 = cache.ops
    F1 = integrator.fsalfirst
    # The caching version uses a special formula to save computation
    # Compute U2
    A_mul_B!(rtmp, phi1, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    # Compute G2 - G1, storing result in the cache F2
    f.f2(rtmp, uprev, p, t)
    f.f2(F2, tmp, p, t + dt)
    F2 .-= rtmp # "F2" is G2 - G1
    # Update u
    u .= tmp
    Base.axpy!(dt, A_mul_B!(rtmp, phi2, F2), u)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK3ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  Au = A * uprev
  F1 = integrator.fsalfirst
  if alg.krylov
    # Krylov on F1 (first column)
    # TODO: reuse Krylov subspace for w1_half
    w1_half = phiv(dt/2, A, F1, 1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    w1 = phiv(dt, A, F1, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U2 = uprev + dt/2 * w1_half[:, 2]
    F2 = _compute_nl(f, U2, p, t + dt/2, A) + Au
    # Krylov on F2 (second column)
    w2 = phiv(dt, A, F2, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U3 = uprev + dt * (2w2[:, 2] - w1[:, 2])
    F3 = _compute_nl(f, U3, p, t + dt, A) + Au
    # Krylov on F3 (third column)
    w3 = phiv(dt, A, F3, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    u = uprev + dt * (4w1[:,4] - 3w1[:,3] + w1[:,2]
                      -8w2[:,4] + 4w2[:,3]
                      +4w3[:,4] - w3[:,3])
  else
    A21, A3, B1, B2, B3 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    U2 = uprev + dt * (A21 * F1)
    F2 = f.f2(U2, p, t + dt/2) + Au
    # stage 3
    U3 = uprev + dt * (A3 * (2F2 - F1))
    F3 = f.f2(U3, p, t + dt) + Au
    # update u
    u = uprev + dt * (B1 * F1 + B2 * F2 + B3 * F3)
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,Au,F2,F3,Jcache,Ks,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  F1 = integrator.fsalfirst
  A_mul_B!(Au, A, uprev)
  halfdt = dt/2
  if alg.krylov
    w1_half, w1, w2, w3, phiv_caches = KsCache
    # Krylov for F1 (first column)
    arnoldi!(Ks, A, F1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w1_half, halfdt, Ks, 1; caches=phiv_caches)
    phiv!(w1, dt, Ks, 3; caches=phiv_caches)
    @muladd @. @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
    _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp); F2 .+= Au
    # Krylov for F2 (second column)
    arnoldi!(Ks, A, F2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w2, dt, Ks, 3; caches=phiv_caches)
    @muladd @. @views tmp = uprev + dt * (2*w2[:, 2] - w1[:, 2]) # tmp is U3
    _compute_nl!(F3, f, tmp, p, t + dt, A, rtmp); F3 .+= Au
    # Krylov for F3 (third column)
    arnoldi!(Ks, A, F3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w3, dt, Ks, 3; caches=phiv_caches)
    # Update u
    @views @. rtmp = 4w1[:,4] - 3w1[:,3] + w1[:,2] - 8w2[:,4] + 4w2[:,3] + 4w3[:,4] - w3[:,3]
    @muladd @. u = uprev + dt * rtmp
  else
    A21, A3, B1, B2, B3 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    A_mul_B!(rtmp, A21, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    f.f2(F2, tmp, p, t + halfdt); F2 .+= Au
    # stage 3
    @muladd @. F3 = 2 * F2 - F1 # use F3 temporarily as cache
    A_mul_B!(rtmp, A3, F3)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U3
    f.f2(F3, tmp, p, t + dt); F3 .+= Au
    # update u
    u .= uprev
    Base.axpy!(dt, A_mul_B!(rtmp, B1, F1), u)
    Base.axpy!(dt, A_mul_B!(rtmp, B2, F2), u)
    Base.axpy!(dt, A_mul_B!(rtmp, B3, F3), u)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  Au = A * uprev
  F1 = integrator.fsalfirst
  halfdt = dt/2
  if alg.krylov # TODO: reuse Krylov subspace for halfdt
    # Krylov on F1 (first column)
    w1_half = phiv(halfdt, A, F1, 1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    w1 = phiv(dt, A, F1, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U2 = uprev + halfdt * w1_half[:, 2]
    F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
    # Krylov on F2 (second column)
    w2_half = phiv(halfdt, A, F2, 1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    w2 = phiv(dt, A, F2, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U3 = uprev + halfdt * w2_half[:, 2]
    F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
    # Krylov on F3 (third column)
    w3_half = phiv(halfdt, A, F3, 1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    w3 = phiv(dt, A, F3, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    # Extra Krylov for computing F4
    rtmp = w1_half[:, 1] - F1 # (exp(hA/2) - I)F1
    wtmp = phiv(halfdt, A, rtmp, 1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U4 = uprev + halfdt * wtmp[:, 2] + dt * w3_half[:, 2]
    F4 = _compute_nl(f, U4, p, t + dt, A) + Au
    # Krylov on F4 (fourth column)
    w4 = phiv(dt, A, F4, 3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    # update u
    u = uprev + dt * (w1[:,2] - 3w1[:,3] + 4w1[:,4] + 2w2[:,3] - 4w2[:4] + 
                      2w3[:,3] - 4w3[:,4] + 4w4[:,4] - w4[:,3])
  else
    A21, A41, A43, B1, B2, B4 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    U2 = uprev + dt * (A21 * F1)
    F2 = f.f2(U2, p, t + halfdt) + Au
    # stage 3
    U3 = uprev + dt * (A21 * F2) # A32 = A21
    F3 = f.f2(U3, p, t + halfdt) + Au
    # stage 4
    U4 = uprev + dt * (A41 * F1 + A43 * F3)
    F4 = f.f2(U4, p, t) + Au
    # update u
    u = uprev + dt * (B1 * F1 + B2 * (F2 + F3) + B4 * F4) # B3 = B2
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,Au,F2,F3,F4,Jcache,Ks,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  F1 = integrator.fsalfirst
  A_mul_B!(Au, A, uprev)
  halfdt = dt/2
  if alg.krylov
    w1_half, w2_half, w1, w2, w3, w4, phiv_caches = KsCache
    # Krylov for F1 (first column)
    arnoldi!(Ks, A, F1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w1_half, halfdt, Ks, 1; caches=phiv_caches)
    phiv!(w1, dt, Ks, 3; caches=phiv_caches)
    @muladd @. @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
    _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp); F2 .+= Au
    # Krylov for F2 (second column)
    arnoldi!(Ks, A, F2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w2_half, halfdt, Ks, 1; caches=phiv_caches)
    phiv!(w2, dt, Ks, 3; caches=phiv_caches)
    @muladd @. @views tmp = uprev + halfdt * w2_half[:, 2] # tmp is U3
    _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp); F3 .+= Au
    # Krylov for F3 (third column)
    w3_half = w2_half # w2_half is no longer used
    arnoldi!(Ks, A, F3; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w3_half, halfdt, Ks, 1; caches=phiv_caches)
    phiv!(w3, dt, Ks, 3; caches=phiv_caches)
    # Extra Krylov for computing F4
    @. @views rtmp = w1_half[:, 1] - F1 # rtmp is (exp(hA/2) - I)F1
    arnoldi!(Ks, A, rtmp; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w1_half, dt, Ks, 1; caches=phiv_caches)
    @views @. rtmp = 0.5w1_half[:, 2] + w3_half[:, 2]
    @muladd @. tmp = uprev + dt * rtmp # tmp is U4
    _compute_nl!(F4, f, tmp, p, t + dt, A, rtmp); F4 .+= Au
    # Krylov for F4 (fourth column)
    arnoldi!(Ks, A, F4; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w4, dt, Ks, 3; caches=phiv_caches)
    # update u
    @views @. rtmp = w1[:,2] - 3w1[:,3] + 4w1[:,4] + 2w2[:,3] - 4w2[:4] + 
                     2w3[:,3] - 4w3[:,4] + 4w4[:,4] - w4[:,3]
    @muladd @. u = uprev + dt * rtmp
  else
    A21, A41, A43, B1, B2, B4 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    A_mul_B!(rtmp, A21, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    f.f2(F2, tmp, p, t + halfdt); F2 .+= Au
    # stage 3
    A_mul_B!(rtmp, A21, F2) # A32 = A21
    @muladd @. tmp = uprev + dt * rtmp # tmp is U3
    f.f2(F3, tmp, p, t + halfdt); F3 .+= Au
    # stage 4
    @. tmp = uprev
    Base.axpy!(dt, A_mul_B!(rtmp, A41, F1), tmp)
    Base.axpy!(dt, A_mul_B!(rtmp, A43, F3), tmp) # tmp is U4
    f.f2(F4, tmp, p, t); F4 .+= Au
    # update u
    u .= uprev
    Base.axpy!(dt, A_mul_B!(rtmp, B1, F1), u)
    F2 .+= F3; Base.axpy!(dt, A_mul_B!(rtmp, B2, F2), u) # B3 = B2
    Base.axpy!(dt, A_mul_B!(rtmp, B4, F4), u)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

######################################################
# Multistep exponential integrators
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
