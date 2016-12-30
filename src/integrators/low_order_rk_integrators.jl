function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{BS3,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  local k1::rateType
  local k2::rateType
  local k3::rateType
  local k4::rateType
  local utilde::uType
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      k2 = f(t+c1*dt,u+dt*a21*k1)
      k3 = f(t+c2*dt,u+dt*a32*k2)
      utmp = u+dt*(a41*k1+a42*k2+a43*k3)
      k4 = f(t+dt,utmp); fsallast = k4
      if integrator.opts.adaptive
        utilde = u + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
        EEst = abs( ((utilde-utmp)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      else
        u = utmp
      end
      if integrator.opts.calck
        k = fsallast
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{BS3,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  uidx = eachindex(u)

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.tableau,integrator.uprev,integrator.kprev)
  @unpack k1,k2,k3,k4,utilde,tmp,atmp = cache

  k = fsallast
  k1 = fsalfirst # done by pointers, no copying
  k4 = fsallast

  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        tmp[i] = u[i]+dt*a21*k1[i]
      end
      f(t+c1*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*a32*k2[i]
      end
      f(t+c2*dt,tmp,k3)
      for i in uidx
        utmp[i] = u[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
      end
      f(t+dt,utmp,k4)
      if integrator.opts.adaptive
        for i in uidx
          utilde[i] = u[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i])
          atmp[i] = ((utilde[i]-utmp[i])/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      else
        recursivecopy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{BS5,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = constructBS5(uEltypeNoUnits)
  local k1::rateType
  local k2::rateType
  local k3::rateType
  local k4::rateType
  local k5::rateType
  local k6::rateType
  local k7::rateType
  local k8::rateType
  local utilde::uType
  local EEst2::uEltypeNoUnits
  integrator.kshortsize = 8
  if integrator.opts.calck
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,zero(rateType))
    end

    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
      for i in 1:3 # Make it full-sized
        push!(integrator.kprev,zero(rateType))
      end
    end
  end
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      k2 = f(t+c1*dt,u+dt*a21*k1)
      k3 = f(t+c2*dt,u+dt*(a31*k1+a32*k2))
      k4 = f(t+c3*dt,u+dt*(a41*k1+a42*k2+a43*k3))
      k5 = f(t+c4*dt,u+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
      k6 = f(t+c5*dt,u+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
      k7 = f(t+dt,u+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6))
      utmp = u+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
      fsallast = f(t+dt,utmp); k8 = fsallast
      if integrator.opts.adaptive
        uhat   = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
        utilde = u + dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
        EEst1 = abs( sum(((uhat)./(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))))
        EEst2 = abs( sum(((utilde-utmp)./(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))))
        EEst = max(EEst1,EEst2)
      else
        u = utmp
      end
      if integrator.opts.calck
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{BS5,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = constructBS5(uEltypeNoUnits)
  integrator.kshortsize = 8
  local EEst2::uEltypeNoUnits
  uidx = eachindex(u)

  if integrator.opts.calck
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,similar(rate_prototype))
    end

    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
      for i in 1:3 # Make it full-sized
        push!(integrator.kprev,similar(rate_prototype))
      end
    end
  end

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.tableau,integrator.uprev,integrator.kprev)
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde = cache

  if integrator.opts.calck
    k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8
  end
  fsalfirst = k1; fsallast = k8  # setup pointers
  f(t,u,k1) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        tmp[i] = u[i]+dt*a21*k1[i]
      end
      f(t+c1*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a31*k1[i]+a32*k2[i])
      end
      f(t+c2*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
      end
      f(t+c3*dt,tmp,k4)
      for i in uidx
        tmp[i] = u[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
      end
      f(t+c4*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
      end
      f(t+c5*dt,tmp,k6)
      for i in uidx
        tmp[i] = u[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
      end
      f(t+dt,tmp,k7)
      for i in uidx
        utmp[i] = u[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
      end
      f(t+dt,utmp,k8)
      if integrator.opts.adaptive
        for i in uidx
          uhat[i]   = dt*(bhat1*k1[i] + bhat3*k3[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i])
          utilde[i] = u[i] + dt*(btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i])
          atmp[i] = ((uhat[i])./(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
          atmptilde[i] = ((utilde[i]-utmp[i])./(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst1 = integrator.opts.internalnorm(atmp)
        EEst2 = integrator.opts.internalnorm(atmptilde)
        EEst = max(EEst1,EEst2)
      else
        recursivecopy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{Tsit5,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = constructTsit5(uEltypeNoUnits)
  local k1::rateType
  local k2::rateType
  local k3::rateType
  local k4::rateType
  local k5::rateType
  local k6::rateType
  local k7::rateType
  local utilde::uType
  integrator.kshortsize = 7
  if integrator.opts.calck
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,zero(rateType))
    end

    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
    end
  end
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      k2 = f(t+c1*dt,u+dt*(a21*k1))
      k3 = f(t+c2*dt,u+dt*(a31*k1+a32*k2))
      k4 = f(t+c3*dt,u+dt*(a41*k1+a42*k2+a43*k3))
      k5 = f(t+c4*dt,u+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
      k6 = f(t+dt,u+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
      utmp = u+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
      fsallast = f(t+dt,utmp); k7 = fsallast
      if integrator.opts.adaptive
        utilde = u + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
        EEst = abs(((utilde-utmp)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      else
        u = utmp
      end
      if integrator.opts.calck
        k[1] = k1
        k[2] = k2
        k[3] = k3
        k[4] = k4
        k[5] = k5
        k[6] = k6
        k[7] = k7
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{Tsit5,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = constructTsit5(uEltypeNoUnits)

  integrator.kshortsize = 7
  uidx = eachindex(u)

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.tableau,integrator.uprev,integrator.kprev)
  @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache

  k1 = fsalfirst; k7 = fsallast # setup pointers

  if integrator.opts.calck
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,similar(rate_prototype))
    end

    # Setup k pointers
    k[1] = k1
    k[2] = k2
    k[3] = k3
    k[4] = k4
    k[5] = k5
    k[6] = k6
    k[7] = k7
    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
    end
  end
  f(t,u,k1) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        tmp[i] = u[i]+dt*(a21*k1[i])
      end
      f(t+c1*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a31*k1[i]+a32*k2[i])
      end
      f(t+c2*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
      end
      f(t+c3*dt,tmp,k4)
      for i in uidx
        tmp[i] = u[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
      end
      f(t+c4*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
      end
      f(t+dt,tmp,k6)
      for i in uidx
        utmp[i] = u[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
      end
      f(t+dt,utmp,k7)
      if integrator.opts.adaptive
        for i in uidx
          utilde[i] = u[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
          atmp[i] = ((utilde[i]-utmp[i])./(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      else
        recursivecopy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{DP5,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  local k1::rateType
  local k2::rateType
  local k3::rateType
  local k4::rateType
  local k5::rateType
  local k6::rateType
  local k7::rateType
  local update::rateType
  local bspl::rateType
  integrator.kshortsize = 4
  if integrator.opts.calck
    d1,d3,d4,d5,d6,d7 = DP5_dense_ds(uEltypeNoUnits)
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,zero(rateType))
    end

    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
    end
  end
  local utilde::uType
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t<T
      @ode_loopheader
      k1 = fsalfirst
      k2 = f(t+c1*dt,u+dt*(a21*k1))
      k3 = f(t+c2*dt,u+dt*(a31*k1+a32*k2))
      k4 = f(t+c3*dt,u+dt*(a41*k1+a42*k2+a43*k3))
      k5 = f(t+c4*dt,u+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
      k6 = f(t+dt,u+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
      update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
      utmp = u+dt*update
      fsallast = f(t+dt,utmp); k7 = fsallast

      if integrator.opts.adaptive
        utilde = u + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
        EEst = abs( ((utilde-utmp)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      else
        u = utmp
      end
      if integrator.opts.calck
        k[1] = update
        bspl = k1 - update
        k[2] = bspl
        k[3] = update - k7 - bspl
        k[4] = (d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end


function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{DP5,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  update = similar(rate_prototype)
  bspl = similar(rate_prototype)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx = eachindex(u)
  integrator.kshortsize = 4
  if integrator.opts.calck
    d1,d3,d4,d5,d6,d7 = DP5_dense_ds(uEltypeNoUnits)
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,similar(rate_prototype))
    end

    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
    end
    # Setup pointers
    k[1] = update
  end
  k1 = fsalfirst; k7 = fsallast
  if !(typeof(integrator.opts.callback)<:Void)
    if integrator.opts.calck
      cache = (u,k...,k1,k2,k3,k4,k5,k6,k7,tmp,utmp,atmp,utilde,bspl,integrator.uprev,integrator.kprev...)
    else
      cache = (u,k1,k2,k3,k4,k5,k6,k7,tmp,utmp,atmp,utilde,update,bspl,integrator.uprev)
    end
  end
  f(t,u,fsalfirst);  # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        tmp[i] = u[i]+dt*(a21*k1[i])
      end
      f(t+c1*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a31*k1[i]+a32*k2[i])
      end
      f(t+c2*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
      end
      f(t+c3*dt,tmp,k4)
      for i in uidx
        tmp[i] =u[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
      end
      f(t+c4*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
      end
      f(t+dt,tmp,k6)
      for i in uidx
        update[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
        utmp[i] = u[i]+dt*update[i]
      end
      f(t+dt,utmp,k7);
      if integrator.opts.adaptive
        for i in uidx
          utilde[i] = u[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
          atmp[i] = ((utilde[i]-utmp[i])/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      else
        recursivecopy!(u, utmp)
      end
      if integrator.opts.calck
        for i in uidx
          bspl[i] = k1[i] - update[i]
          k[2][i] = bspl[i]
          k[3][i] = update[i] - k7[i] - bspl[i]
          k[4][i] = (d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i])
        end
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end
