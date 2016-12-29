function ode_solve{uType<:AbstractArray,tType,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{DP5Threaded,uType,tType,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k2::rateType = similar(rate_prototype)
  k3::rateType = similar(rate_prototype)
  k4::rateType = similar(rate_prototype)
  k5::rateType = similar(rate_prototype)
  k6::rateType = similar(rate_prototype)
  update::rateType = similar(rate_prototype)
  bspl::rateType = similar(rate_prototype)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx::Base.OneTo{Int64} = eachindex(u)
  integrator.kshortsize = 4
  if integrator.opts.calck
    d1,d3,d4,d5,d6,d7 = DP5_dense_ds(uEltypeNoUnits)
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,similar(rate_prototype))
    end

    # Setup k pointers
    k[1] = update
    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
    end
  end
  k1 = fsalfirst; k7 = fsallast
  f(t,u,fsalfirst);  # Pre-start fsal
  if integrator.custom_callback
    if integrator.opts.calck
      cache = (u,k...,k1,k2,k3,k4,k5,k6,k7,tmp,utmp,atmp,utilde,bspl,integrator.uprev,integrator.kprev...)
    else
      cache = (u,k1,k2,k3,k4,k5,k6,k7,tmp,utmp,atmp,utilde,bspl,integrator.uprev,integrator.kprev...)
    end
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      dp5threaded_loop1(dt,tmp,u,a21,k1,uidx)
      f(t+c1*dt,tmp,k2)
      dp5threaded_loop2(dt,tmp,u,a31,k1,a32,k2,uidx)
      f(t+c2*dt,tmp,k3)
      dp5threaded_loop3(dt,tmp,u,a41,k1,a42,k2,a43,k3,uidx)
      f(t+c3*dt,tmp,k4)
      dp5threaded_loop4(dt,tmp,u,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
      f(t+c4*dt,tmp,k5)
      dp5threaded_loop5(dt,tmp,u,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
      f(t+dt,tmp,k6)
      dp5threaded_loop6(dt,utmp,u,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,update,uidx)
      f(t+dt,utmp,fsallast)
      if integrator.opts.adaptive
        dp5threaded_adaptiveloop(dt,utilde,u,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,atmp,utmp,integrator.opts.abstol,integrator.opts.reltol,uidx)
        EEst = integrator.opts.internalnorm(atmp)
      else
        recursivecopy!(u, utmp)
      end
      if integrator.opts.calck
        dp5threaded_denseloop(bspl,update,k1,k3,k4,k5,k6,k7,k,d1,d3,d4,d5,d6,d7,uidx)
      end
      @ode_loopfooter
    end
  end
  @ode_postamble
end

@noinline function dp5threaded_denseloop(bspl,update,k1,k3,k4,k5,k6,k7,k,d1,d3,d4,d5,d6,d7,uidx)
  Threads.@threads for i in uidx
    bspl[i] = k1[i] - update[i]
    k[2][i] = bspl[i]
    k[3][i] = update[i] - k7[i] - bspl[i]
    k[4][i] = (d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i])
  end
end

@noinline function dp5threaded_loop1(dt,tmp,u,a21,k1,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+dt*(a21*k1[i])
  end
end

@noinline function dp5threaded_loop2(dt,tmp,u,a31,k1,a32,k2,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+dt*(a31*k1[i]+a32*k2[i])
  end
end

@noinline function dp5threaded_loop3(dt,tmp,u,a41,k1,a42,k2,a43,k3,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
end

@noinline function dp5threaded_loop4(dt,tmp,u,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
  Threads.@threads for i in uidx
    tmp[i] =u[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
end

@noinline function dp5threaded_loop5(dt,tmp,u,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
end

@noinline function dp5threaded_loop6(dt,utmp,u,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,update,uidx)
  Threads.@threads for i in uidx
    update[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
    utmp[i] = u[i]+dt*update[i]
  end
end

@noinline function dp5threaded_adaptiveloop(dt,utilde,u,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,atmp,utmp,abstol,reltol,uidx)
  Threads.@threads for i in uidx
    utilde[i] = u[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
    atmp[i] = ((utilde[i]-utmp[i])/(abstol+max(abs(u[i]),abs(utmp[i]))*reltol))
  end
end
