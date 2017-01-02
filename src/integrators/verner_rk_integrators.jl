function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern6,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType;
  local utilde::uType; integrator.fsalfirst = f(t,u) # Pre-start fsal
  integrator.kshortsize = 9
  k = ksEltype(integrator.kshortsize)
  integrator.k = k

  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k1 = integrator.fsalfirst
      k2 = f(t+c1*dt,u+dt*(a21*k1))
      k3 = f(t+c2*dt,u+dt*(a31*k1+a32*k2))
      k4 = f(t+c3*dt,u+dt*(a41*k1       +a43*k3))
      k5 = f(t+c4*dt,u+dt*(a51*k1       +a53*k3+a54*k4))
      k6 = f(t+c5*dt,u+dt*(a61*k1       +a63*k3+a64*k4+a65*k5))
      k7 = f(t+c6*dt,u+dt*(a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6))
      k8 = f(t+dt,u+dt*(a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7))
      utmp =    u+dt*(a91*k1              +a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
      integrator.fsallast = f(t+dt,utmp); k9 = integrator.fsallast
      if integrator.opts.adaptive
        utilde = u + dt*(b1*k1 + b4*k4 + b5*k5 + b6*k6 + b7*k7 + b8*k8 + b9*k9)
        EEst = abs( ((utilde-utmp)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      end
      if integrator.opts.calck
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end



function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern6,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)


  @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,atmp = integrator.cache

  integrator.kshortsize = 9
  uidx = eachindex(u)
  integrator.fsalfirst = k1 ; integrator.fsallast = k9
  if integrator.opts.calck
    k = ksEltype(integrator.kshortsize)
    k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9 # Set the pointers
    integrator.k = k
  end
  f(t,u,k1) # Pre-start fsal
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
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
      f(t+dt,tmp,k8)
      for i in uidx
        utmp[i]=u[i]+dt*(a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i])
      end
      f(t+dt,utmp,k9)
      if integrator.opts.adaptive
        for i in uidx
          utilde[i] = u[i] + dt*(b1*k1[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i])
          atmp[i] = ((utilde[i]-utmp[i])/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern7,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType; local k10::rateType; local utilde::uType; local update::uType
  integrator.kshortsize = 10
  k = ksEltype(integrator.kshortsize)
  integrator.k = k
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k1 = f(t,u)
      k2 = f(t+c2*dt,u+dt*(a021*k1))
      k3 = f(t+c3*dt,u+dt*(a031*k1+a032*k2))
      k4 = f(t+c4*dt,u+dt*(a041*k1       +a043*k3))
      k5 = f(t+c5*dt,u+dt*(a051*k1       +a053*k3+a054*k4))
      k6 = f(t+c6*dt,u+dt*(a061*k1       +a063*k3+a064*k4+a065*k5))
      k7 = f(t+c7*dt,u+dt*(a071*k1       +a073*k3+a074*k4+a075*k5+a076*k6))
      k8 = f(t+c8*dt,u+dt*(a081*k1       +a083*k3+a084*k4+a085*k5+a086*k6+a087*k7))
      k9 = f(t+dt,u+dt*(a091*k1          +a093*k3+a094*k4+a095*k5+a096*k6+a097*k7+a098*k8))
      k10= f(t+dt,u+dt*(a101*k1          +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7))
      update = dt*(k1*b1 + k4*b4 + k5*b5 + k6*b6 + k7*b7 + k8*b8 + k9*b9)
      utmp = u + update
      if integrator.opts.adaptive
        EEst = abs( ((update - dt*(bhat1*k1 + bhat4*k4 + bhat5*k5 + bhat6*k6 + bhat7*k7 + bhat10*k10))/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      end
      if integrator.opts.calck
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern7,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnits)
  uidx = eachindex(u);


  @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,update,utilde,tmp,atmp = integrator.cache

  integrator.kshortsize = 10
  if integrator.opts.calck
    k = ksEltype(integrator.kshortsize)
    k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10 # Setup pointers
    integrator.k = k
  end

  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      f(t,u,k1)
      for i in uidx
        tmp[i] = u[i]+dt*(a021*k1[i])
      end
      f(t+c2*dt,tmp,k2)
      for i in uidx
        tmp[i] = u[i]+dt*(a031*k1[i]+a032*k2[i])
      end
      f(t+c3*dt,tmp,k3)
      for i in uidx
        tmp[i] = u[i]+dt*(a041*k1[i]+a043*k3[i])
      end
      f(t+c4*dt,tmp,k4)
      for i in uidx
        tmp[i] = u[i]+dt*(a051*k1[i]+a053*k3[i]+a054*k4[i])
      end
      f(t+c5*dt,tmp,k5)
      for i in uidx
        tmp[i] = u[i]+dt*(a061*k1[i]+a063*k3[i]+a064*k4[i]+a065*k5[i])
      end
      f(t+c6*dt,tmp,k6)
      for i in uidx
        tmp[i] = u[i]+dt*(a071*k1[i]+a073*k3[i]+a074*k4[i]+a075*k5[i]+a076*k6[i])
      end
      f(t+c7*dt,tmp,k7)
      for i in uidx
        tmp[i] = u[i]+dt*(a081*k1[i]+a083*k3[i]+a084*k4[i]+a085*k5[i]+a086*k6[i]+a087*k7[i])
      end
      f(t+c8*dt,tmp,k8)
      for i in uidx
        tmp[i] = u[i]+dt*(a091*k1[i]+a093*k3[i]+a094*k4[i]+a095*k5[i]+a096*k6[i]+a097*k7[i]+a098*k8[i])
      end
      f(t+dt,tmp,k9)
      for i in uidx
        tmp[i] = u[i]+dt*(a101*k1[i]+a103*k3[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i])
      end
      f(t+dt,tmp,k10)
      for i in uidx
        update[i] = dt*(k1[i]*b1 + k4[i]*b4 + k5[i]*b5 + k6[i]*b6 + k7[i]*b7 + k8[i]*b8 + k9[i]*b9)
        utmp[i] = u[i] + update[i]
      end
      if integrator.opts.adaptive
        for i in uidx
          atmp[i] = ((update[i] - dt*(bhat1*k1[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i] + bhat7*k7[i] + bhat10*k10[i]))/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern8,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType; local k10::rateType; local utilde::uType; local update::uType
  integrator.kshortsize = 13
  k = ksEltype(integrator.kshortsize)
  integrator.k = k
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k1 = f(t,u)
      k2 = f(t+c2*dt ,u+dt*(a0201*k1))
      k3 = f(t+c3*dt ,u+dt*(a0301*k1+a0302*k2))
      k4 = f(t+c4*dt ,u+dt*(a0401*k1       +a0403*k3))
      k5 = f(t+c5*dt ,u+dt*(a0501*k1       +a0503*k3+a0504*k4))
      k6 = f(t+c6*dt ,u+dt*(a0601*k1                +a0604*k4+a0605*k5))
      k7 = f(t+c7*dt ,u+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6))
      k8 = f(t+c8*dt ,u+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7))
      k9 = f(t+c9*dt ,u+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8))
      k10= f(t+c10*dt,u+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9))
      k11= f(t+c11*dt,u+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10))
      k12= f(t+    dt,u+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11))
      k13= f(t+    dt,u+dt*(a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10))
      update = dt*(k1*b1 + k6*b6 + k7*b7 + k8*b8 + k9*b9 + k10*b10 + k11*b11 + k12*b12)
      utmp = u + update
      if integrator.opts.adaptive
        EEst = abs( ((update - dt*(bhat1*k1 + bhat6*k6 + bhat7*k7 + bhat8*k8 + bhat9*k9 + bhat10*k10 + bhat13*k13))/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol)))
      end
      if integrator.opts.calck
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern8,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
  uidx = eachindex(u);


  @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp = integrator.cache

  integrator.kshortsize = 13
  if integrator.opts.calck
    k = ksEltype(integrator.kshortsize)
    k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13 # Setup pointers
    integrator.k = k
  end

  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      f(t,u,k1)
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
        tmp[i] = u[i]+dt*(a1301*k1[i]+a1304*k4[i]+a1305*k5[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i])
      end
      f(t+dt,tmp,k13)
      for i in uidx
        update[i] = dt*(k1[i]*b1 + k6[i]*b6 + k7[i]*b7 + k8[i]*b8 + k9[i]*b9 + k10[i]*b10 + k11[i]*b11 + k12[i]*b12)
        utmp[i] = u[i] + update[i]
      end
      if integrator.opts.adaptive
        for i in uidx
          atmp[i] = ((update[i] - dt*(bhat1*k1[i] + bhat6*k6[i] + bhat7*k7[i] + bhat8*k8[i] + bhat9*k9[i] + bhat10*k10[i] + bhat13*k13[i]))/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern9,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  local k1::rateType; local k2::rateType; local k3::rateType; local k4::rateType;
  local k5::rateType; local k6::rateType; local k7::rateType; local k8::rateType;
  local k9::rateType; local k10::rateType; local k11::rateType; local k12::rateType;
  local k13::rateType; local k14::rateType; local k15::rateType; local k16::rateType;
  local utilde::uType; local update::uType
  integrator.kshortsize = 16
  if integrator.opts.calck
    k = ksEltype(integrator.kshortsize)
    integrator.k = k
  end
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k1 = f(t,u)
      k2 = f(t+c1*dt,u+dt*(a0201*k1))
      k3 = f(t+c2*dt,u+dt*(a0301*k1+a0302*k2))
      k4 = f(t+c3*dt,u+dt*(a0401*k1       +a0403*k3))
      k5 = f(t+c4*dt,u+dt*(a0501*k1       +a0503*k3+a0504*k4))
      k6 = f(t+c5*dt,u+dt*(a0601*k1                +a0604*k4+a0605*k5))
      k7 = f(t+c6*dt,u+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6))
      k8 = f(t+c7*dt,u+dt*(a0801*k1                                  +a0806*k6+a0807*k7))
      k9 = f(t+c8*dt,u+dt*(a0901*k1                                  +a0906*k6+a0907*k7+a0908*k8))
      k10 =f(t+c9*dt,u+dt*(a1001*k1                                  +a1006*k6+a1007*k7+a1008*k8+a1009*k9))
      k11= f(t+c10*dt,u+dt*(a1101*k1                                  +a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10))
      k12= f(t+c11*dt,u+dt*(a1201*k1                                  +a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11))
      k13= f(t+c12*dt,u+dt*(a1301*k1                                  +a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10+a1311*k11+a1312*k12))
      k14= f(t+c13*dt,u+dt*(a1401*k1                                  +a1406*k6+a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13))
      k15= f(t+dt,u+dt*(a1501*k1                                  +a1506*k6+a1507*k7+a1508*k8+a1509*k9+a1510*k10+a1511*k11+a1512*k12+a1513*k13+a1514*k14))
      k16= f(t+dt,u+dt*(a1601*k1                                  +a1606*k6+a1607*k7+a1608*k8+a1609*k9+a1610*k10+a1611*k11+a1612*k12+a1613*k13))
      update = dt*(k1*b1+k8*b8+k9*b9+k10*b10+k11*b11+k12*b12+k13*b13+k14*b14+k15*b15)
      utmp = u + update
      if integrator.opts.adaptive
        EEst = abs((update - dt*(k1*bhat1 + k8*bhat8 + k9*bhat9 + k10*bhat10 + k11*bhat11 + k12*bhat12 + k13*bhat13 + k16*bhat16))/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))
      end
      if integrator.opts.calck
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13;k[14]=k14;k[15]=k15;k[16]=k16
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Vern9,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)

  uidx = eachindex(u)

  @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,update,tmp,atmp = integrator.cache
  integrator.kshortsize = 16
  if integrator.opts.calck
    k = ksEltype(integrator.kshortsize)
    k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13;k[14]=k14;k[15]=k15;k[16]=k16 # Setup pointers
    integrator.k = k
  end
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
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
        tmp[i] = u[i]+dt*(a0801*k1[i]+a0806*k6[i]+a0807*k7[i])
      end
      f(t+c7*dt,tmp,k8)
      for i in uidx
        tmp[i] = u[i]+dt*(a0901*k1[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i])
      end
      f(t+c8*dt,tmp,k9)
      for i in uidx
        tmp[i] = u[i]+dt*(a1001*k1[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i])
      end
      f(t+c9*dt,tmp,k10)
      for i in uidx
        tmp[i] = u[i]+dt*(a1101*k1[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i])
      end
      f(t+c10*dt,tmp,k11)
      for i in uidx
        tmp[i] = u[i]+dt*(a1201*k1[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i])
      end
      f(t+c11*dt,tmp,k12)
      for i in uidx
        tmp[i] = u[i]+dt*(a1301*k1[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i]+a1311*k11[i]+a1312*k12[i])
      end
      f(t+c12*dt,tmp,k13)
      for i in uidx
        tmp[i] = u[i]+dt*(a1401*k1[i]+a1406*k6[i]+a1407*k7[i]+a1408*k8[i]+a1409*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*k13[i])
      end
      f(t+c13*dt,tmp,k14)
      for i in uidx
        tmp[i] = u[i]+dt*(a1501*k1[i]+a1506*k6[i]+a1507*k7[i]+a1508*k8[i]+a1509*k9[i]+a1510*k10[i]+a1511*k11[i]+a1512*k12[i]+a1513*k13[i]+a1514*k14[i])
      end
      f(t+dt,tmp,k15)
      for i in uidx
        tmp[i] = u[i]+dt*(a1601*k1[i]+a1606*k6[i]+a1607*k7[i]+a1608*k8[i]+a1609*k9[i]+a1610*k10[i]+a1611*k11[i]+a1612*k12[i]+a1613*k13[i])
      end
      f(t+dt,tmp,k16)
      for i in uidx
        update[i] = dt*(k1[i]*b1+k8[i]*b8+k9[i]*b9+k10[i]*b10+k11[i]*b11+k12[i]*b12+k13[i]*b13+k14[i]*b14+k15[i]*b15)
        utmp[i] = u[i] + update[i]
      end
      if integrator.opts.adaptive
        for i in uidx
          atmp[i] = ((update[i] - dt*(k1[i]*bhat1 + k8[i]*bhat8 + k9[i]*bhat9 + k10[i]*bhat10 + k11[i]*bhat11 + k12[i]*bhat12 + k13[i]*bhat13 + k16[i]*bhat16))/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end
