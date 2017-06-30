@inline function initialize!(integrator,cache::DP5ThreadedCache,f=integrator.f)
  integrator.kshortsize = 4
  integrator.k = [cache.update,cache.bspl,cache.dense_tmp3,cache.dense_tmp4]
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k7
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

@inline function perform_step!(integrator,cache::DP5ThreadedCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache.tab
  dp5threaded_loop1(dt,tmp,uprev,a21,k1,uidx)
  f(@muladd(t+c1*dt),tmp,k2)
  dp5threaded_loop2(dt,tmp,uprev,a31,k1,a32,k2,uidx)
  f(@muladd(t+c2*dt),tmp,k3)
  dp5threaded_loop3(dt,tmp,uprev,a41,k1,a42,k2,a43,k3,uidx)
  f(@muladd(t+c3*dt),tmp,k4)
  dp5threaded_loop4(dt,tmp,uprev,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
  f(@muladd(t+c4*dt),tmp,k5)
  dp5threaded_loop5(dt,tmp,uprev,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
  f(t+dt,tmp,k6)
  dp5threaded_loop6(dt,u,uprev,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,update,uidx)
  f(t+dt,u,integrator.fsallast)
  if integrator.opts.adaptive
    dp5threaded_adaptiveloop(dt,utilde,uprev,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,atmp,u,integrator.opts.abstol,integrator.opts.reltol,uidx)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  dp5threaded_denseloop(bspl,update,k1,k3,k4,k5,k6,k7,integrator.k,d1,d3,d4,d5,d6,d7,uidx)
  @pack integrator = t,dt,u,k
end

@noinline function dp5threaded_denseloop(bspl,update,k1,k3,k4,k5,k6,k7,k,d1,d3,d4,d5,d6,d7,uidx)
  Threads.@threads for i in uidx
    bspl[i] = k1[i] - update[i]
    k[3][i] = update[i] - k7[i] - bspl[i]
    k[4][i] = (d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i])
  end
end

@noinline function dp5threaded_loop1(dt,tmp,uprev,a21,k1,uidx)
  a = dt*a21
  Threads.@threads for i in uidx
    tmp[i] = @muladd uprev[i]+a*k1[i]
  end
end

@noinline function dp5threaded_loop2(dt,tmp,uprev,a31,k1,a32,k2,uidx)
  Threads.@threads for i in uidx
    tmp[i] = @muladd uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
end

@noinline function dp5threaded_loop3(dt,tmp,uprev,a41,k1,a42,k2,a43,k3,uidx)
  Threads.@threads for i in uidx
    tmp[i] = @muladd uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
end

@noinline function dp5threaded_loop4(dt,tmp,uprev,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
  Threads.@threads for i in uidx
    tmp[i] = @muladd uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
end

@noinline function dp5threaded_loop5(dt,tmp,uprev,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
  Threads.@threads for i in uidx
    tmp[i] = @muladd uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
end

@noinline function dp5threaded_loop6(dt,u,uprev,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,update,uidx)
  Threads.@threads for i in uidx
    update[i] = @muladd a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
    u[i] = @muladd uprev[i]+dt*update[i]
  end
end

@noinline function dp5threaded_adaptiveloop(dt,utilde,uprev,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,atmp,u,abstol,reltol,uidx)
  #Threads.@threads for (i,atol,rtol) in zip(uidx,Iterators.cycle(abstol),Iterators.cycle(reltol)) # warnings due to iterators
  for (i,atol,rtol) in zip(uidx,Iterators.cycle(abstol),Iterators.cycle(reltol))
    utilde[i] = @muladd uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
    atmp[i] = ((utilde[i]-u[i])./@muladd(atol+max(abs(uprev[i]),abs(u[i]))*rtol))
  end
end
