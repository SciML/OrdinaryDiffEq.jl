function initialize!(integrator, cache::LawsonEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  rtmp = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  integrator.fsalfirst = rtmp # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = zero(integrator.fsalfirst)
end

function perform_step!(integrator, cache::LawsonEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  rtmp = integrator.fsalfirst
  A = f.f1
  if integrator.alg.krylov
    @muladd u = expmv(dt, A, uprev + dt*rtmp; tol=integrator.opts.reltol, m=min(integrator.alg.m, size(A,1)), norm=normbound)
  else
    @muladd u = expm(dt*A)*(uprev + dt*rtmp)
  end
  rtmp = f.f2(u,p,t+dt)
  k = A*u + rtmp # For the interpolation, needs k at the updated point
  integrator.fsallast = rtmp
  integrator.k[1] = integrator.fsalfirst # this is wrong, since it's just rtmp. Should fsal this value though
  integrator.k[2] = k
  integrator.u = u
end

function initialize!(integrator, cache::LawsonEulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst,rtmp = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = rtmp
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = fsalfirst # this is wrong, since it's just rtmp. Should fsal this value though
  integrator.k[2] = k
  integrator.f.f1(k,integrator.u,integrator.p,integrator.t)
  integrator.f.f2(rtmp,integrator.uprev,integrator.p,integrator.t) # For the interpolation, needs k at the updated point
  @. integrator.fsalfirst = k + rtmp
end

function perform_step!(integrator, cache::LawsonEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,rtmp,tmp = cache
  A = f.f1
  @muladd @. tmp = uprev + dt*integrator.fsalfirst
  if integrator.alg.krylov
    expmv!(u,dt,A,tmp; tol=integrator.opts.reltol, m=min(integrator.alg.m, size(A,1)), norm=normbound)
  else
    A_mul_B!(u,cache.expA,tmp)
  end
  A_mul_B!(tmp,A,u)
  f.f2(rtmp,u,p,t+dt)
  @. k = tmp + rtmp
end

function initialize!(integrator, cache::NorsettEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  rtmp = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  integrator.fsalfirst = rtmp # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = zero(integrator.fsalfirst)
end

function perform_step!(integrator, cache::NorsettEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  rtmp = integrator.fsalfirst
  A = f.f1
  if integrator.alg.krylov
    u = phimv(dt,A,rtmp,uprev; tol=integrator.opts.reltol, m=min(integrator.alg.m, size(A,1)), norm=normbound)
  else
    u = uprev + ((expm(dt*A)-I)/A)*(A*uprev + rtmp)
  end
  rtmp = f.f2(u,p,t+dt)
  k = A*u + rtmp # For the interpolation, needs k at the updated point
  integrator.fsallast = rtmp
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = k
  integrator.u = u
end

function initialize!(integrator, cache::NorsettEulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst,rtmp = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = rtmp
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = fsalfirst
  integrator.k[2] = k
  integrator.f.f1(k,integrator.u,integrator.p,integrator.t)
  integrator.f.f2(rtmp,integrator.uprev,integrator.p,integrator.t) # For the interpolation, needs k at the updated point
  @. integrator.fsalfirst = k + rtmp
end

function perform_step!(integrator, cache::NorsettEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,rtmp,tmp = cache
  A = f.f1

  if integrator.alg.krylov
    phimv!(u,dt,A,rtmp,uprev; tol=integrator.opts.reltol, m=min(integrator.alg.m, size(A,1)), norm=normbound)
  else
    A_mul_B!(tmp,A,uprev)
    tmp .+= rtmp
    A_mul_B!(rtmp,cache.phi1,tmp)
    @. u = uprev + rtmp
  end
  A_mul_B!(tmp,A,u)
  f.f2(rtmp,u,p,t+dt)
  @. k = tmp +  rtmp
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
  integrator.fsallast = ETD2Fsal(zero(lin), zero(nl), zero(nlprev))
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(lin) + zero(nl)
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

function initialize!(integrator, cache::ETDRK4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  rtmp = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  integrator.fsalfirst = rtmp # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = zero(integrator.fsalfirst)
end

function perform_step!(integrator, cache::ETDRK4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack E,E2,a,b,c,Q = cache
  rtmp = integrator.fsalfirst
  A = f.f1

  tmp = E2*uprev

  k1 = integrator.f.f2(uprev,p,t)
  s1 = tmp + Q*k1;
  k2 = integrator.f.f2(s1,p,t+dt/2)
  s2 = tmp + Q*k2;
  k3 = integrator.f.f2(s2,p,t+dt/2)
  s3 = E2*s1 + Q*(2*k3-k1);
  k4 = integrator.f.f2(s3,p,t+dt)
  u = E*uprev + a*k1 + 2b*(k2+k3) + c*k4;


  integrator.fsallast = integrator.f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::ETDRK4Cache)
  integrator.kshortsize = 2
  @unpack tmp,fsalfirst,tmp2 = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = tmp2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = fsalfirst
  integrator.k[2] = tmp2
  integrator.f.f1(tmp,integrator.u,integrator.p,integrator.t)
  integrator.f.f2(tmp2,integrator.uprev,integrator.p,integrator.t) # For the interpolation, needs k at the updated point
  @. integrator.fsalfirst = tmp + tmp2
end

function perform_step!(integrator, cache::ETDRK4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp2,tmp = cache
  @unpack E,E2,a,b,c,Q = cache
  @unpack k1,k2,k3,k4,s1 = cache
  A = f.f1

  # Substep 1
  integrator.f.f2(k1,uprev,p,t) # TODO: Erase for faslfirst
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

  integrator.f(tmp2, u, p, t+dt)
end
