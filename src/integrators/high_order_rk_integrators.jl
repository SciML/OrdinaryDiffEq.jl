function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{TanYam7,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructTanYam7(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType; local k10::rateType;
  local utilde::uType;
  if integrator.opts.calck
    pop!(integrator.sol.k) # Get rid of the one it starts with
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k = f(t,u)
      k1 = k
      k2 = f(t+c1*dt,u+dt*(a21*k1))
      k3 = f(t+c2*dt,u+dt*(a31*k1+a32*k2))
      k4 = f(t+c3*dt,u+dt*(a41*k1       +a43*k3))
      k5 = f(t+c4*dt,u+dt*(a51*k1       +a53*k3+a54*k4))
      k6 = f(t+c5*dt,u+dt*(a61*k1       +a63*k3+a64*k4+a65*k5))
      k7 = f(t+c6*dt,u+dt*(a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6))
      k8 = f(t+c7*dt,u+dt*(a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7))
      k9 = f(t+dt,u+dt*(a91*k1       +a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8))
      k10= f(t+dt,u+dt*(a101*k1      +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8))
      utmp = u + dt*(k1*b1+k4*b4+k5*b5+k6*b6+k7*b7+k8*b8+k9*b9)
      if integrator.opts.adaptive
        utilde = u + dt*(k1*bhat1+k4*bhat4+k5*bhat5+k6*bhat6+k7*bhat7+k8*bhat8+k10*bhat10)
        EEst = abs( ((utilde-utmp)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  if integrator.opts.calck
    k = f(t,u)
    push!(integrator.sol.k,k)
  end
  ode_postamble!(integrator)
  nothing
end


function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{TanYam7,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructTanYam7(uEltypeNoUnits)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype) ; k3 = similar(rate_prototype); k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6 = similar(rate_prototype) ; k7 = similar(rate_prototype); k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10= similar(rate_prototype) ;
  k = similar(rate_prototype)
  if integrator.calcprevs && integrator.opts.calck
    integrator.kprev = similar(rate_prototype)
  end
  utilde = similar(u); uidx = eachindex(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits)

  if integrator.opts.calck
    pop!(integrator.sol.k) # Get rid of the one it starts with
  end
  k = k1
  if integrator.custom_callback
    if integrator.opts.calck
      cache = (u,k,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,integrator.uprev,integrator.kprev,utmp,tmp,atmp)
    else
      cache = (u,k,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,integrator.uprev,utmp,tmp,atmp)
    end
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      for i in uidx
        tmp[i] = u[i]+dt*(a21*k1[i])
      end
      f(t+c1*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a31*k1[i]+a32*k2[i])
      end
      f(t+c2*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a41*k1[i]+a43*k3[i])
      end
      f(t+c3*dt,tmp,k4)
      for i in uidx
        tmp[i] = u[i]+dt*(a51*k1[i]+a53*k3[i]+a54*k4[i])
      end
      f(t+c4*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
      end
      f(t+c5*dt,tmp,k6)
      for i in uidx
        tmp[i] = u[i]+dt*(a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
      end
      f(t+c6*dt,tmp,k7)
      for i in uidx
        tmp[i] = u[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
      end
      f(t+c7*dt,tmp,k8)
      for i in uidx
        tmp[i] = u[i]+dt*(a91*k1[i]+a93*k3[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i])
      end
      f(t+dt,tmp,k9)
      for i in uidx
        tmp[i] = u[i]+dt*(a101*k1[i]+a103*k3[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i])
      end
      f(t+dt,tmp,k10)
      for i in uidx
        utmp[i] = u[i] + dt*(k1[i]*b1+k4[i]*b4+k5[i]*b5+k6[i]*b6+k7[i]*b7+k8[i]*b8+k9[i]*b9)
      end
      if integrator.opts.adaptive
        for i in uidx
          utilde[i] = u[i] + dt*(k1[i]*bhat1+k4[i]*bhat4+k5[i]*bhat5+k6[i]*bhat6+k7[i]*bhat7+k8[i]*bhat8+k10[i]*bhat10)
          atmp[i] = ((utilde[i]-utmp[i])/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      else
        recursivecopy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  if integrator.opts.calck
    f(t,u,k)
    push!(integrator.sol.k,deepcopy(k))
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{DP8,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = constructDP8(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType; local k10::rateType; local k11::rateType; local k12::rateType;
  local k13::rateType; local utilde::uType; local udiff::rateType; local bspl::rateType
  integrator.kshortsize = 7
  if integrator.opts.calck
    c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = DP8Interp(uEltypeNoUnits)
    d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = DP8Interp_polyweights(uEltypeNoUnits)
    fsal = true
    k = ksEltype()
    for i in 1:integrator.kshortsize
      push!(k,zero(rateType))
    end

    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
    end
  else
    fsal = false
  end
  if fsal # Pre-start fsal
    k1 = f(t,u)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      if fsal
        k1 = fsalfirst
      else
        k1 = f(t,u)
      end
      k2 = f(t+c2*dt,u+dt*(a0201*k1))
      k3 = f(t+c3*dt,u+dt*(a0301*k1+a0302*k2))
      k4 = f(t+c4*dt,u+dt*(a0401*k1       +a0403*k3))
      k5 = f(t+c5*dt,u+dt*(a0501*k1       +a0503*k3+a0504*k4))
      k6 = f(t+c6*dt,u+dt*(a0601*k1                +a0604*k4+a0605*k5))
      k7 = f(t+c7*dt,u+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6))
      k8 = f(t+c8*dt,u+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7))
      k9 = f(t+c9*dt,u+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8))
      k10 =f(t+c10*dt,u+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9))
      k11= f(t+c11*dt,u+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10))
      k12= f(t+dt,u+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11))
      kupdate= b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
      update = dt*kupdate
      utmp = u + update
      if integrator.opts.adaptive
        err5 = abs(dt*(k1*er1 + k6*er6 + k7*er7 + k8*er8 + k9*er9 + k10*er10 + k11*er11 + k12*er12)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)) # Order 5
        err3 = abs((update - dt*(bhh1*k1 + bhh2*k9 + bhh3*k12))/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      if integrator.opts.calck
        k13 = f(t+dt,utmp)
        fsallast = k13
        k14 = f(t+c14*dt,u+dt*(a1401*k1         +a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13))
        k15 = f(t+c15*dt,u+dt*(a1501*k1+a1506*k6+a1507*k7+a1508*k8                   +a1511*k11+a1512*k12+a1513*k13+a1514*k14))
        k16 = f(t+c16*dt,u+dt*(a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9                              +a1613*k13+a1614*k14+a1615*k15))
        udiff = kupdate
        k[1] = udiff
        bspl = k1 - udiff
        k[2] = bspl
        k[3] = udiff - k13 - bspl
        k[4] = (d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16)
        k[5] = (d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16)
        k[6] = (d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16)
        k[7] = (d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{DP8,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = constructDP8(uEltypeNoUnits)
  k1 = similar(rate_prototype); k2  = similar(rate_prototype); k3  = similar(rate_prototype);  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6  = similar(rate_prototype); k7  = similar(rate_prototype);  k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype); k12 = similar(rate_prototype)
  kupdate = similar(rate_prototype); utilde = similar(u);
  #err5 = similar(u); err3 = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u); atmp2 = similar(u,uEltypeNoUnits); update = similar(u)
  local k13::rateType; local k14::rateType; local k15::rateType; local k16::rateType;
  local udiff::rateType; local bspl::rateType
  integrator.kshortsize = 7
  if integrator.opts.calck
    c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = DP8Interp(uEltypeNoUnits)
    d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = DP8Interp_polyweights(uEltypeNoUnits)
    fsal = true
    if integrator.opts.calck
      k = ksEltype()
      for i in 1:integrator.kshortsize
        push!(k,similar(rate_prototype))
      end

      if integrator.calcprevs
        integrator.kprev = deepcopy(k)
      end
    end
    k13 = similar(rate_prototype)
    k14 = similar(rate_prototype)
    k15 = similar(rate_prototype)
    k16 = similar(rate_prototype)
    udiff = similar(rate_prototype)
    bspl = similar(rate_prototype)
    fsalfirst = k1
    fsallast = k13
  else

    fsal = false
  end

  if integrator.custom_callback
    if integrator.opts.calck
      cache = (u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k...,integrator.kprev...,integrator.uprev,udiff,bspl,utilde,update,utmp,tmp,atmp,atmp2,kupdate)
    else
      cache = (u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,integrator.uprev,utilde,update,utmp,tmp,atmp,atmp2,kupdate)
    end
  end

  if fsal # Pre-start FSAL
    f(t,u,k1)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      if !fsal
        f(t,u,k1)
      end
      for i in uidx
        tmp[i] = u[i]+dt*(a0201*k1[i])
      end
      f(t+c2*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a0301*k1[i]+a0302*k2[i])
      end
      f(t+c3*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a0401*k1[i]+a0403*k3[i])
      end
      f(t+c4*dt,tmp,k4)
      for i in uidx
        tmp[i] = u[i]+dt*(a0501*k1[i]+a0503*k3[i]+a0504*k4[i])
      end
      f(t+c5*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a0601*k1[i]+a0604*k4[i]+a0605*k5[i])
      end
      f(t+c6*dt,tmp,k6)
      for i in uidx
        tmp[i] = u[i]+dt*(a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i])
      end
      f(t+c7*dt,tmp,k7)
      for i in uidx
        tmp[i] = u[i]+dt*(a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i])
      end
      f(t+c8*dt,tmp,k8)
      for i in uidx
        tmp[i] = u[i]+dt*(a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i])
      end
      f(t+c9*dt,tmp,k9)
      for i in uidx
        tmp[i] = u[i]+dt*(a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i])
      end
      f(t+c10*dt,tmp,k10)
      for i in uidx
        tmp[i] = u[i]+dt*(a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i])
      end
      f(t+c11*dt,tmp,k11)
      for i in uidx
        tmp[i] = u[i]+dt*(a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i])
      end
      f(t+dt,tmp,k12)
      for i in uidx
        kupdate[i] = b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i]
        update[i] = dt*kupdate[i]
        utmp[i] = u[i] + update[i]
      end
      if integrator.opts.adaptive
        for i in uidx
          atmp[i] = (dt*(k1[i]*er1 + k6[i]*er6 + k7[i]*er7 + k8[i]*er8 + k9[i]*er9 + k10[i]*er10 + k11[i]*er11 + k12[i]*er12)/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
          atmp2[i]= ((update[i] - dt*(bhh1*k1[i] + bhh2*k9[i] + bhh3*k12[i]))/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        err5 = integrator.opts.internalnorm(atmp) # Order 5
        err3 = integrator.opts.internalnorm(atmp2) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        recursivecopy!(u, utmp)
      end
      if integrator.opts.calck
        f(t+dt,utmp,k13)
        for i in uidx
          tmp[i] = u[i]+dt*(a1401*k1[i]+a1407*k7[i]+a1408*k8[i]+a1409*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*k13[i])
        end
        f(t+c14*dt,tmp,k14)
        for i in uidx
          tmp[i] = u[i]+dt*(a1501*k1[i]+a1506*k6[i]+a1507*k7[i]+a1508*k8[i]+a1511*k11[i]+a1512*k12[i]+a1513*k13[i]+a1514*k14[i])
        end
        f(t+c15*dt,tmp,k15)
        for i in uidx
          tmp[i] = u[i]+dt*(a1601*k1[i]+a1606*k6[i]+a1607*k7[i]+a1608*k8[i]+a1609*k9[i]+a1613*k13[i]+a1614*k14[i]+a1615*k15[i])
        end
        f(t+c16*dt,tmp,k16)
        for i in uidx
          udiff[i]= kupdate[i]
          k[1][i] = udiff[i]
          bspl[i] = k1[i] - udiff[i]
          k[2][i] = bspl[i]
          k[3][i] = udiff[i] - k13[i] - bspl[i]
          k[4][i] = (d401*k1[i]+d406*k6[i]+d407*k7[i]+d408*k8[i]+d409*k9[i]+d410*k10[i]+d411*k11[i]+d412*k12[i]+d413*k13[i]+d414*k14[i]+d415*k15[i]+d416*k16[i])
          k[5][i] = (d501*k1[i]+d506*k6[i]+d507*k7[i]+d508*k8[i]+d509*k9[i]+d510*k10[i]+d511*k11[i]+d512*k12[i]+d513*k13[i]+d514*k14[i]+d515*k15[i]+d516*k16[i])
          k[6][i] = (d601*k1[i]+d606*k6[i]+d607*k7[i]+d608*k8[i]+d609*k9[i]+d610*k10[i]+d611*k11[i]+d612*k12[i]+d613*k13[i]+d614*k14[i]+d615*k15[i]+d616*k16[i])
          k[7][i] = (d701*k1[i]+d706*k6[i]+d707*k7[i]+d708*k8[i]+d709*k9[i]+d710*k10[i]+d711*k11[i]+d712*k12[i]+d713*k13[i]+d714*k14[i]+d715*k15[i]+d716*k16[i])
        end
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{TsitPap8,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType; local k10::rateType; local k11::rateType; local k12::rateType;
  local k13::rateType; local utilde::uType;
  if integrator.opts.calck
    pop!(integrator.sol.k) # Take out the initial
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k = f(t,u)
      k1 = k
      k2 = f(t+c1*dt,u+dt*(a0201*k1))
      k3 = f(t+c2*dt,u+dt*(a0301*k1+a0302*k2))
      k4 = f(t+c3*dt,u+dt*(a0401*k1       +a0403*k3))
      k5 = f(t+c4*dt,u+dt*(a0501*k1       +a0503*k3+a0504*k4))
      k6 = f(t+c5*dt,u+dt*(a0601*k1                +a0604*k4+a0605*k5))
      k7 = f(t+c6*dt,u+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6))
      k8 = f(t+c7*dt,u+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7))
      k9 = f(t+c8*dt,u+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8))
      k10 =f(t+c9*dt,u+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9))
      k11= f(t+c10*dt,u+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10))
      k12= f(t+dt,u+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11))
      k13= f(t+dt,u+dt*(a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10))
      update = dt*(b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12)
      utmp = u + update
      if integrator.opts.adaptive
        EEst = abs((update - dt*(k1*bhat1 + k6*bhat6 + k7*bhat7 + k8*bhat8 + k9*bhat9 + k10*bhat10 + k13*bhat13))/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  if integrator.opts.calck
    k = f(t,u)
    push!(integrator.sol.k,k)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,TabType,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{TsitPap8,uType,tType,ksEltype,TabType,SolType,rateType,F,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype); k12 = similar(rate_prototype)
  k13 = similar(rate_prototype); update = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
  utilde = similar(u);
  k = similar(rate_prototype)
  if integrator.calcprevs
    integrator.kprev = similar(rate_prototype)
  end
  if integrator.opts.calck
    pop!(integrator.sol.k)
  end
  k = k1
  if integrator.custom_callback
    cache = (u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,integrator.uprev,update,tmp,utmp, atmp,utilde,integrator.kprev)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      for i in uidx
        tmp[i] = u[i]+dt*(a0201*k1[i])
      end
      f(t+c1*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a0301*k1[i]+a0302*k2[i])
      end
      f(t+c2*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a0401*k1[i]+a0403*k3[i])
      end
      f(t+c3*dt,tmp,k4)
      for i in uidx
        tmp[i] = u[i]+dt*(a0501*k1[i]+a0503*k3[i]+a0504*k4[i])
      end
      f(t+c4*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a0601*k1[i]+a0604*k4[i]+a0605*k5[i])
      end
      f(t+c5*dt,tmp,k6)
      for i in uidx
        tmp[i] = u[i]+dt*(a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i])
      end
      f(t+c6*dt,tmp,k7)
      for i in uidx
        tmp[i] = u[i]+dt*(a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i])
      end
      f(t+c7*dt,tmp,k8)
      for i in uidx
        tmp[i] = u[i]+dt*(a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i])
      end
      f(t+c8*dt,tmp,k9)
      for i in uidx
        tmp[i] = u[i]+dt*(a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i])
      end
      f(t+c9*dt,tmp,k10)
      for i in uidx
        tmp[i] = u[i]+dt*(a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i])
      end
      f(t+c10*dt,tmp,k11)
      for i in uidx
        tmp[i] = u[i]+dt*(a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i])
      end
      f(t+dt,tmp,k12)
      for i in uidx
        tmp[i] = u[i]+dt*(a1301*k1[i]+a1304*k4[i]+a1305*k5[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i])
      end
      f(t+dt,tmp,k13)
      for i in uidx
        update[i] = dt*(b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i])
        utmp[i] = u[i] + update[i]
      end
      if integrator.opts.adaptive
        for i in uidx
          atmp[i] = ((update[i] - dt*(k1[i]*bhat1 + k6[i]*bhat6 + k7[i]*bhat7 + k8[i]*bhat8 + k9[i]*bhat9 + k10[i]*bhat10 + k13[i]*bhat13))/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      else
        recursivecopy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  if integrator.opts.calck
    f(t,u,k)
    push!(integrator.sol.k,deepcopy(k))
  end
  ode_postamble!(integrator)
  nothing
end
