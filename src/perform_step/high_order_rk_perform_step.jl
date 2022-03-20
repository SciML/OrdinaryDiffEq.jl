function initialize!(integrator,cache::TanYam7ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TanYam7ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,btilde1,btilde4,btilde5,btilde6,btilde7,btilde8,btilde9,btilde10 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(uprev+a*k1, p, t + c1*dt)
  k3 = f(uprev+dt*(a31*k1+a32*k2), p, t + c2*dt)
  k4 = f(uprev+dt*(a41*k1       +a43*k3), p, t + c3*dt)
  k5 = f(uprev+dt*(a51*k1       +a53*k3+a54*k4), p, t + c4*dt)
  k6 = f(uprev+dt*(a61*k1       +a63*k3+a64*k4+a65*k5), p, t + c5*dt)
  k7 = f(uprev+dt*(a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6), p, t + c6*dt)
  k8 = f(uprev+dt*(a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7), p, t + c7*dt)
  k9 = f(uprev+dt*(a91*k1       +a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8), p, t+dt)
  k10= f(uprev+dt*(a101*k1      +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8), p, t+dt)
  integrator.destats.nf += 9
  u = uprev + dt*(b1*k1+b4*k4+b5*k5+b6*k6+b7*k7+b8*k8+b9*k9)
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1+btilde4*k4+btilde5*k5+btilde6*k6+btilde7*k7+btilde8*k8+btilde9*k9+btilde10*k10)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.fsallast = f(u, p, t+dt) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::TanYam7Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::TanYam7Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack fsalfirst,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k = cache
  @unpack c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,btilde1,btilde4,btilde5,btilde6,btilde7,btilde8,btilde9,btilde10 = cache.tab
  k1 = fsalfirst
  f(k1, uprev, p, t)
  a = dt*a21
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t + c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t + c2*dt)
  @.. tmp = uprev+dt*(a41*k1+a43*k3)
  f(k4, tmp, p, t + c3*dt)
  @.. tmp = uprev+dt*(a51*k1+a53*k3+a54*k4)
  f(k5, tmp, p, t + c4*dt)
  @.. tmp = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t + c5*dt)
  @.. tmp = uprev+dt*(a71*k1+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, tmp, p, t + c6*dt)
  @.. tmp = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  f(k8, tmp, p, t + c7*dt)
  @.. tmp = uprev+dt*(a91*k1+a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
  f(k9, tmp, p, t+dt)
  @.. tmp = uprev+dt*(a101*k1+a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8)
  f(k10, tmp, p, t+dt)
  @.. u = uprev + dt*(b1*k1+b4*k4+b5*k5+b6*k6+b7*k7+b8*k8+b9*k9)
  integrator.destats.nf += 10
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1+btilde4*k4+btilde5*k5+btilde6*k6+btilde7*k7+btilde8*k8+btilde9*k9+btilde10*k10)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
  integrator.destats.nf += 1
  return nothing
end

#=
@muladd function perform_step!(integrator, cache::TanYam7Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack fsalfirst,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k = cache
  @unpack c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,btilde1,btilde4,btilde5,btilde6,btilde7,btilde8,btilde9,btilde10 = cache.tab
  k1 = fsalfirst
  f(k1, uprev, p, t)
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t + c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(k3, tmp, p, t + c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a43*k3[i])
  end
  f(k4, tmp, p, t + c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a53*k3[i]+a54*k4[i])
  end
  f(k5, tmp, p, t + c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(k6, tmp, p, t + c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(k7, tmp, p, t + c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
  end
  f(k8, tmp, p, t + c7*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a91*k1[i]+a93*k3[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i])
  end
  f(k9, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a101*k1[i]+a103*k3[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i])
  end
  f(k10, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + dt*(b1*k1[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i])
  end
  integrator.destats.nf += 10
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i]+btilde4*k4[i]+btilde5*k5[i]+btilde6*k6[i]+btilde7*k7[i]+btilde8*k8[i]+btilde9*k9[i]+btilde10*k10[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end
=#

function initialize!(integrator, cache::DP8ConstantCache)
  integrator.kshortsize = 7
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  @inbounds for i in eachindex(integrator.k)
    integrator.k[i] = zero(integrator.fsalfirst)
  end
end

@muladd function perform_step!(integrator, cache::DP8ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = cache
  k1 = integrator.fsalfirst
  a = dt*a0201
  k2 = f(uprev+a*k1, p, t + c2*dt)
  k3 = f(uprev+dt*(a0301*k1+a0302*k2), p, t + c3*dt)
  k4 = f(uprev+dt*(a0401*k1       +a0403*k3), p, t + c4*dt)
  k5 = f(uprev+dt*(a0501*k1       +a0503*k3+a0504*k4), p, t + c5*dt)
  k6 = f(uprev+dt*(a0601*k1                +a0604*k4+a0605*k5), p, t + c6*dt)
  k7 = f(uprev+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6), p, t + c7*dt)
  k8 = f(uprev+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7), p, t + c8*dt)
  k9 = f(uprev+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8), p, t + c9*dt)
  k10 =f(uprev+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9), p, t + c10*dt)
  k11= f(uprev+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10), p, t + c11*dt)
  k12= f(uprev+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11), p, t+dt)
  integrator.destats.nf += 11
  kupdate= b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
  u = uprev + dt*kupdate
  if integrator.opts.adaptive
    utilde = dt*(k1*er1 + k6*er6 + k7*er7 + k8*er8 + k9*er9 + k10*er10 + k11*er11 + k12*er12)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    err5 = integrator.opts.internalnorm(atmp,t) # Order 5
    utilde = dt*(btilde1*k1 + btilde6*k6 + btilde7*k7 + btilde8*k8 + btilde9*k9 + btilde10*k10 + btilde11*k11 + btilde12*k12)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    err3 = integrator.opts.internalnorm(atmp,t) # Order 3
    err52 = err5*err5
    if err5 ≈ 0 && err3 ≈ 0
      integrator.EEst = zero(integrator.EEst)
    else
      integrator.EEst = err52/sqrt(err52 + 0.01*err3*err3)
    end
  end
  k13 = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.fsallast = k13
  if integrator.opts.calck
    @unpack c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = cache
    @unpack d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = cache
    k14 = f(uprev+dt*(a1401*k1         +a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13), p, t + c14*dt)
    k15 = f(uprev+dt*(a1501*k1+a1506*k6+a1507*k7+a1508*k8                   +a1511*k11+a1512*k12+a1513*k13+a1514*k14), p, t + c15*dt)
    k16 = f(uprev+dt*(a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9                              +a1613*k13+a1614*k14+a1615*k15), p, t + c16*dt)
    integrator.destats.nf += 3
    udiff = kupdate
    integrator.k[1] = udiff
    bspl = k1 - udiff
    integrator.k[2] = bspl
    integrator.k[3] = udiff - k13 - bspl
    integrator.k[4] = d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16
    integrator.k[5] = d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16
    integrator.k[6] = d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16
    integrator.k[7] = d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16
  end
  integrator.u = u
end

function initialize!(integrator, cache::DP8Cache)
  integrator.kshortsize = 7
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [cache.udiff,cache.bspl,cache.dense_tmp3,cache.dense_tmp4,cache.dense_tmp5,cache.dense_tmp6,cache.dense_tmp7]
  integrator.fsalfirst = cache.k1
  integrator.fsallast = cache.k13
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::DP8Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,kupdate,utilde,tmp,atmp = cache
  f(k1, uprev, p, t)
  a = dt*a0201
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t + c2*dt)
  @.. tmp = uprev+dt*(a0301*k1+a0302*k2)
  f(k3, tmp, p, t + c3*dt)
  @.. tmp = uprev+dt*(a0401*k1+a0403*k3)
  f(k4, tmp, p, t + c4*dt)
  @.. tmp = uprev+dt*(a0501*k1+a0503*k3+a0504*k4)
  f(k5, tmp, p, t + c5*dt)
  @.. tmp = uprev+dt*(a0601*k1+a0604*k4+a0605*k5)
  f(k6, tmp, p, t + c6*dt)
  @.. tmp = uprev+dt*(a0701*k1+a0704*k4+a0705*k5+a0706*k6)
  f(k7, tmp, p, t + c7*dt)
  @.. tmp = uprev+dt*(a0801*k1+a0804*k4+a0805*k5+a0806*k6+a0807*k7)
  f(k8, tmp, p, t + c8*dt)
  @.. tmp = uprev+dt*(a0901*k1+a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8)
  f(k9, tmp, p, t + c9*dt)
  @.. tmp = uprev+dt*(a1001*k1+a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9)
  f(k10, tmp, p, t + c10*dt)
  @.. tmp = uprev+dt*(a1101*k1+a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
  f(k11, tmp, p, t + c11*dt)
  @.. tmp = uprev+dt*(a1201*k1+a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
  f(k12, tmp, p, t+dt)
  @.. kupdate = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
  @.. u = uprev + dt*kupdate
  integrator.destats.nf += 12
  if integrator.opts.adaptive
    @.. utilde = dt*(k1*er1 + k6*er6 + k7*er7 + k8*er8 + k9*er9 + k10*er10 + k11*er11 + k12*er12)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    err5 = integrator.opts.internalnorm(atmp,t) # Order 5
    @.. utilde= dt*(btilde1*k1 + btilde6*k6 + btilde7*k7 + btilde8*k8 + btilde9*k9 + btilde10*k10 + btilde11*k11 + btilde12*k12)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    err3 = integrator.opts.internalnorm(atmp,t) # Order 3
    err52 = err5*err5
    if err5 ≈ 0 && err3 ≈ 0
      integrator.EEst = zero(integrator.EEst)
    else
      integrator.EEst = err52/sqrt(err52 + 0.01*err3*err3)
    end
  end
  f(k13, u, p, t+dt)
  integrator.destats.nf += 1
  if integrator.opts.calck
    @unpack c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = cache.tab
    @unpack d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = cache.tab
    @.. tmp = uprev+dt*(a1401*k1+a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13)
    f(k14, tmp, p, t + c14*dt)
    @.. tmp = uprev+dt*(a1501*k1+a1506*k6+a1507*k7+a1508*k8+a1511*k11+a1512*k12+a1513*k13+a1514*k14)
    f(k15, tmp, p, t + c15*dt)
    @.. tmp = uprev+dt*(a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9+a1613*k13+a1614*k14+a1615*k15)
    f(k16, tmp, p, t + c16*dt)
    integrator.destats.nf += 3
    @.. udiff= kupdate
    @.. bspl = k1 - udiff
    @.. integrator.k[3] = udiff - k13 - bspl
    @.. integrator.k[4] = d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16
    @.. integrator.k[5] = d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16
    @.. integrator.k[6] = d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16
    @.. integrator.k[7] = d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16
  end
end

#=
@muladd function perform_step!(integrator, cache::DP8Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,kupdate,utilde,tmp,atmp = cache
  f(k1, uprev, p, t)
  a = dt*a0201
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t + c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0301*k1[i]+a0302*k2[i])
  end
  f(k3, tmp, p, t + c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0401*k1[i]+a0403*k3[i])
  end
  f(k4, tmp, p, t + c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0501*k1[i]+a0503*k3[i]+a0504*k4[i])
  end
  f(k5, tmp, p, t + c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0601*k1[i]+a0604*k4[i]+a0605*k5[i])
  end
  f(k6, tmp, p, t + c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i])
  end
  f(k7, tmp, p, t + c7*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i])
  end
  f(k8, tmp, p, t + c8*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i])
  end
  f(k9, tmp, p, t + c9*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i])
  end
  f(k10, tmp, p, t + c10*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i])
  end
  f(k11, tmp, p, t + c11*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i])
  end
  f(k12, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds kupdate[i] = b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i]
    @inbounds u[i] = uprev[i] + dt*kupdate[i]
  end
  integrator.destats.nf += 12
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(k1[i]*er1 + k6[i]*er6 + k7[i]*er7 + k8[i]*er8 + k9[i]*er9 + k10[i]*er10 + k11[i]*er11 + k12[i]*er12)
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    err5 = integrator.opts.internalnorm(atmp,t) # Order 5
    @tight_loop_macros for i in uidx
      @inbounds utilde[i]= dt*(btilde1*k1[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i] + btilde9*k9[i] + btilde10*k10[i] + btilde11*k11[i] + btilde12*k12[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    err3 = integrator.opts.internalnorm(atmp,t) # Order 3
    err52 = err5*err5
    if err5 ≈ 0 && err3 ≈ 0
      integrator.EEst = zero(integrator.EEst)
    else
      integrator.EEst = err52/sqrt(err52 + 0.01*err3*err3)
    end
  end
  f(k13, u, p, t+dt)
  integrator.destats.nf += 1
  if integrator.opts.calck
    @unpack c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = cache.tab
    @unpack d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = cache.tab
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1401*k1[i]+a1407*k7[i]+a1408*k8[i]+a1409*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*k13[i])
    end
    f(k14, tmp, p, t + c14*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1501*k1[i]+a1506*k6[i]+a1507*k7[i]+a1508*k8[i]+a1511*k11[i]+a1512*k12[i]+a1513*k13[i]+a1514*k14[i])
    end
    f(k15, tmp, p, t + c15*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1601*k1[i]+a1606*k6[i]+a1607*k7[i]+a1608*k8[i]+a1609*k9[i]+a1613*k13[i]+a1614*k14[i]+a1615*k15[i])
    end
    f(k16, tmp, p, t + c16*dt)
    integrator.destats.nf += 3
    @tight_loop_macros for i in uidx
      @inbounds udiff[i]= kupdate[i]
      @inbounds bspl[i] = k1[i] - udiff[i]
      @inbounds integrator.k[3][i] = udiff[i] - k13[i] - bspl[i]
      @inbounds integrator.k[4][i] = d401*k1[i]+d406*k6[i]+d407*k7[i]+d408*k8[i]+d409*k9[i]+d410*k10[i]+d411*k11[i]+d412*k12[i]+d413*k13[i]+d414*k14[i]+d415*k15[i]+d416*k16[i]
      @inbounds integrator.k[5][i] = d501*k1[i]+d506*k6[i]+d507*k7[i]+d508*k8[i]+d509*k9[i]+d510*k10[i]+d511*k11[i]+d512*k12[i]+d513*k13[i]+d514*k14[i]+d515*k15[i]+d516*k16[i]
      @inbounds integrator.k[6][i] = d601*k1[i]+d606*k6[i]+d607*k7[i]+d608*k8[i]+d609*k9[i]+d610*k10[i]+d611*k11[i]+d612*k12[i]+d613*k13[i]+d614*k14[i]+d615*k15[i]+d616*k16[i]
      @inbounds integrator.k[7][i] = d701*k1[i]+d706*k6[i]+d707*k7[i]+d708*k8[i]+d709*k9[i]+d710*k10[i]+d711*k11[i]+d712*k12[i]+d713*k13[i]+d714*k14[i]+d715*k15[i]+d716*k16[i]
    end
  end
end
=#

function initialize!(integrator, cache::TsitPap8ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

#=
@muladd function perform_step!(integrator, cache::TsitPap8ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,btilde13 = cache
  k1 = integrator.fsalfirst
  a = dt*a0201
  k2 = f(t + c1*dt, @.. uprev+a*k1)
  k3 = f(t + c2*dt, @.. uprev+dt*(a0301*k1+a0302*k2))
  k4 = f(t + c3*dt, @.. uprev+dt*(a0401*k1       +a0403*k3))
  k5 = f(t + c4*dt, @.. uprev+dt*(a0501*k1       +a0503*k3+a0504*k4))
  k6 = f(t + c5*dt, @.. uprev+dt*(a0601*k1                +a0604*k4+a0605*k5))
  k7 = f(t + c6*dt, @.. uprev+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6))
  k8 = f(t + c7*dt, @.. uprev+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7))
  k9 = f(t + c8*dt, @.. uprev+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8))
  k10 =f(t + c9*dt, @.. uprev+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9))
  k11= f(t + c10*dt, @.. uprev+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10))
  k12= f(t+dt, @.. uprev+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11))
  k13= f(t+dt, @.. uprev+dt*(a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10))
  u = @.. uprev + dt*(b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12)
  if integrator.opts.adaptive
    utilde = @.. dt*(btilde1*k1 + btilde6*k6 + btilde7*k7 + btilde8*k8 + btilde9*k9 + btilde10*k10 + btilde11*k11 + btilde12*k12 + btilde13*k13)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
=#

@muladd function perform_step!(integrator, cache::TsitPap8ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,btilde13 = cache
  k1 = integrator.fsalfirst
  a = dt*a0201
  k2 = f(uprev+a*k1, p, t + c1*dt)
  k3 = f(uprev+dt*(a0301*k1+a0302*k2), p, t + c2*dt)
  k4 = f(uprev+dt*(a0401*k1       +a0403*k3), p, t + c3*dt)
  k5 = f(uprev+dt*(a0501*k1       +a0503*k3+a0504*k4), p, t + c4*dt)
  k6 = f(uprev+dt*(a0601*k1                +a0604*k4+a0605*k5), p, t + c5*dt)
  k7 = f(uprev+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6), p, t + c6*dt)
  k8 = f(uprev+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7), p, t + c7*dt)
  k9 = f(uprev+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8), p, t + c8*dt)
  k10 =f(uprev+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9), p, t + c9*dt)
  k11= f(uprev+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10), p, t + c10*dt)
  k12= f(uprev+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11), p, t+dt)
  k13= f(uprev+dt*(a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10), p, t+dt)
  integrator.destats.nf += 12
  u = uprev + dt*(b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12)
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde6*k6 + btilde7*k7 + btilde8*k8 + btilde9*k9 + btilde10*k10 + btilde11*k11 + btilde12*k12 + btilde13*k13)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end


function initialize!(integrator, cache::TsitPap8Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::TsitPap8Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,btilde13 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,atmp,k = cache
  k1 = cache.fsalfirst
  f(k1, uprev, p, t)
  a = dt*a0201
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t + c1*dt)
  @.. tmp = uprev+dt*(a0301*k1+a0302*k2)
  f(k3, tmp, p, t + c2*dt)
  @.. tmp = uprev+dt*(a0401*k1+a0403*k3)
  f(k4, tmp, p, t + c3*dt)
  @.. tmp = uprev+dt*(a0501*k1+a0503*k3+a0504*k4)
  f(k5, tmp, p, t + c4*dt)
  @.. tmp = uprev+dt*(a0601*k1+a0604*k4+a0605*k5)
  f(k6, tmp, p, t + c5*dt)
  @.. tmp = uprev+dt*(a0701*k1+a0704*k4+a0705*k5+a0706*k6)
  f(k7, tmp, p, t + c6*dt)
  @.. tmp = uprev+dt*(a0801*k1+a0804*k4+a0805*k5+a0806*k6+a0807*k7)
  f(k8, tmp, p, t + c7*dt)
  @.. tmp = uprev+dt*(a0901*k1+a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8)
  f(k9, tmp, p, t + c8*dt)
  @.. tmp = uprev+dt*(a1001*k1+a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9)
  f(k10, tmp, p, t + c9*dt)
  @.. tmp = uprev+dt*(a1101*k1+a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
  f(k11, tmp, p, t + c10*dt)
  @.. tmp = uprev+dt*(a1201*k1+a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
  f(k12, tmp, p, t+dt)
  @.. tmp = uprev+dt*(a1301*k1+a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10)
  f(k13, tmp, p, t+dt)
  @.. u = uprev + dt*(b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12)
  integrator.destats.nf += 13
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde6*k6 + btilde7*k7 + btilde8*k8 + btilde9*k9 + btilde10*k10 + btilde11*k11 + btilde12*k12 + btilde13*k13)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
  integrator.destats.nf += 1
  return nothing
end

#=
@muladd function perform_step!(integrator, cache::TsitPap8Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,btilde1,btilde6,btilde7,btilde8,btilde9,btilde10,btilde11,btilde12,btilde13 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,atmp,k = cache
  k1 = cache.fsalfirst
  f(k1, uprev, p, t)
  a = dt*a0201
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t + c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0301*k1[i]+a0302*k2[i])
  end
  f(k3, tmp, p, t + c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0401*k1[i]+a0403*k3[i])
  end
  f(k4, tmp, p, t + c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0501*k1[i]+a0503*k3[i]+a0504*k4[i])
  end
  f(k5, tmp, p, t + c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0601*k1[i]+a0604*k4[i]+a0605*k5[i])
  end
  f(k6, tmp, p, t + c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i])
  end
  f(k7, tmp, p, t + c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i])
  end
  f(k8, tmp, p, t + c7*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i])
  end
  f(k9, tmp, p, t + c8*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i])
  end
  f(k10, tmp, p, t + c9*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i])
  end
  f(k11, tmp, p, t + c10*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i])
  end
  f(k12, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a1301*k1[i]+a1304*k4[i]+a1305*k5[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i])
  end
  f(k13, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + dt*(b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i])
  end
  integrator.destats.nf += 13
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i] + btilde9*k9[i] + btilde10*k10[i] + btilde11*k11[i] + btilde12*k12[i] + btilde13*k13[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end
=#

function initialize!(integrator, cache::PFRK87ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::PFRK87ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α0201, α0301, α0401, α0501, α0601, α0701, α0302, α0403, α0503, α0504, α0604, α0704, α0605, α0705, α0706, α0908, α1008, α1108, α1208, α1308, α1009, α1109, α1209, α1309, α1110, α1210, α1310, α1211, α1311, β1, β6, β7, β8, β9, β10, β11, β12, β13, β1tilde, β6tilde, β7tilde, β8tilde, β9tilde, β10tilde, β11tilde, β12tilde, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 = cache
  alg = unwrap_alg(integrator, false)
  ν = alg.omega*dt
  νsq = ν^2

  c_ν = -0.19781108078634084 - (0.164050909125528499*νsq) + (0.042578310088756321*(νsq^2)) - (0.002300513610963998*(νsq^3)) + (0.000033467244551879287*(νsq^4)) - (7.8661142036921924*(1/100000000)*(νsq^5))
  d_ν = 1-(0.296457092123567400*νsq) + (0.0015793885907465726*(νsq^2)) - (0.00018913011771688527*(νsq^3)) + (0.000017089234650765179*(νsq^4)) - (1.2705211682518626*(1/10000000)*(νsq^5))

  α0807 = c_ν/d_ν
  α0801 = 0.026876256 + 0.0576576*α0807
  α0804 = 0.22464336 - 0.944944*α0807
  α0805 = 0.000369024 - 0.2061696*α0807
  α0806 = 0.21311136 + 0.093456*α0807
  α0901 = 0.07239997637512857 + 0.01913119863380767*α0807
  α0904 = -0.688400520601143 - 0.3135390887207368*α0807
  α0905 = -0.17301267570583073 - 0.06840852844816077*α0807
  α0906 = 0.1440060555560846 + 0.031009360422930017*α0807
  α0907 = 0.9982362892760762 + 0.33180705811215994*α0807
  α1001 = 0.16261514523236525 - 0.12125171966747463*α0807
  α1004 = -2.1255544052061124 + 1.9871809612169453*α0807
  α1005 = -0.216403903283323 + 0.43356675517460624*α0807
  α1006 = -0.060417230254934076 - 0.1965343807796979*α0807
  α1007 = 2.4846281621788395 - 2.102961615944379*α0807
  α1101 = -1.0320124180911034 + 1.061943768952537*α0807
  α1104 = 13.666683232895137 - 17.40407843561103*α0807
  α1105 = 0.25990355211486116 - 3.797253476860588*α0807
  α1106 = -5.759316475814002 + 1.7212824826428488*α0807
  α1107 = -12.822511612651839 + 18.41810566087623*α0807
  α1201 = 0.2478349764611783 - 0.06383934946543009*α0807
  α1204 = -4.593782880309185 + 1.046256005127882*α0807
  α1205 = -0.39566692537411896 + 0.22827403748244698*α0807
  α1206 = -3.0673550479691665 - 0.10347586863902129*α0807
  α1207 = 5.386688702227177 - 1.1072148245058775*α0807
  α1301 = 0.7332242174431163 - 0.5164807626867616*α0807
  α1304 = -10.196728938160977 + 8.464545832921925*α0807
  α1305 = -0.43865244706547707 + 1.846809999910238*α0807
  α1306 = 0.5693856884667226 - 0.8371528845746959*α0807
  α1307 = 10.52865228002416 - 8.957722185570706*α0807

  k1 = integrator.fsalfirst
  k2 = f(uprev+dt*α0201*k1, p, t + c2*dt)
  k3 = f(uprev+dt*(α0301*k1+α0302*k2), p, t + c3*dt)
  k4 = f(uprev+dt*(α0401*k1       +α0403*k3), p, t + c4*dt)
  k5 = f(uprev+dt*(α0501*k1       +α0503*k3+α0504*k4), p, t + c5*dt)
  k6 = f(uprev+dt*(α0601*k1                +α0604*k4+α0605*k5), p, t + c6*dt)
  k7 = f(uprev+dt*(α0701*k1                +α0704*k4+α0705*k5+α0706*k6), p, t + c7*dt)
  k8 = f(uprev+dt*(α0801*k1                +α0804*k4+α0805*k5+α0806*k6+α0807*k7), p, t + c8*dt)
  k9 = f(uprev+dt*(α0901*k1                +α0904*k4+α0905*k5+α0906*k6+α0907*k7+α0908*k8), p, t + c9*dt)
  k10 =f(uprev+dt*(α1001*k1                +α1004*k4+α1005*k5+α1006*k6+α1007*k7+α1008*k8+α1009*k9), p, t + c10*dt)
  k11= f(uprev+dt*(α1101*k1                +α1104*k4+α1105*k5+α1106*k6+α1107*k7+α1108*k8+α1109*k9+α1110*k10), p, t + c11*dt)
  k12= f(uprev+dt*(α1201*k1                +α1204*k4+α1205*k5+α1206*k6+α1207*k7+α1208*k8+α1209*k9+α1210*k10+α1211*k11), p, t+c12*dt)
  k13= f(uprev+dt*(α1301*k1                +α1304*k4+α1305*k5+α1306*k6+α1307*k7+α1308*k8+α1309*k9+α1310*k10+α1311*k11), p, t+c13*dt)
  integrator.destats.nf += 12
  u = uprev + dt*(β1*k1+β6*k6+β7*k7+β8*k8+β9*k9+β10*k10+β11*k11+β12*k12+β13*k13)
  if integrator.opts.adaptive
    utilde = dt*(β1tilde*k1 + β6tilde*k6 + β7tilde*k7 + β8tilde*k8 + β9tilde*k9 + β10tilde*k10 + β11tilde*k11 + β12tilde*k12)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end


function initialize!(integrator, cache::PFRK87Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::PFRK87Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α0201, α0301, α0401, α0501, α0601, α0701, α0302, α0403, α0503, α0504, α0604, α0704, α0605, α0705, α0706, α0908, α1008, α1108, α1208, α1308, α1009, α1109, α1209, α1309, α1110, α1210, α1310, α1211, α1311, β1, β6, β7, β8, β9, β10, β11, β12, β13, β1tilde, β6tilde, β7tilde, β8tilde, β9tilde, β10tilde, β11tilde, β12tilde, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,atmp,k = cache

  alg = unwrap_alg(integrator, false)
  ν = alg.omega*dt
  νsq = ν^2

  c_ν = -0.19781108078634084 - (0.164050909125528499*νsq) + (0.042578310088756321*(νsq^2)) - (0.002300513610963998*(νsq^3)) + (0.000033467244551879287*(νsq^4)) - (7.8661142036921924*(10^(-8))*(νsq^5))
  d_ν = 1-(0.296457092123567400*νsq) + (0.0015793885907465726*(νsq^2)) - (0.00018913011771688527*(νsq^3)) + (0.000017089234650765179*(νsq^4)) - (1.2705211682518626*(10^(-7))*(νsq^5))

  α0807 = c_ν/d_ν
  α0801 = 0.026876256 + 0.0576576*α0807
  α0804 = 0.22464336 - 0.944944*α0807
  α0805 = 0.000369024 - 0.2061696*α0807
  α0806 = 0.21311136 + 0.093456*α0807
  α0901 = 0.07239997637512857 + 0.01913119863380767*α0807
  α0904 = -0.688400520601143 - 0.3135390887207368*α0807
  α0905 = -0.17301267570583073 - 0.06840852844816077*α0807
  α0906 = 0.1440060555560846 + 0.031009360422930017*α0807
  α0907 = 0.9982362892760762 + 0.33180705811215994*α0807
  α1001 = 0.16261514523236525 - 0.12125171966747463*α0807
  α1004 = -2.1255544052061124 + 1.9871809612169453*α0807
  α1005 = -0.216403903283323 + 0.43356675517460624*α0807
  α1006 = -0.060417230254934076 - 0.1965343807796979*α0807
  α1007 = 2.4846281621788395 - 2.102961615944379*α0807
  α1101 = -1.0320124180911034 + 1.061943768952537*α0807
  α1104 = 13.666683232895137 - 17.40407843561103*α0807
  α1105 = 0.25990355211486116 - 3.797253476860588*α0807
  α1106 = -5.759316475814002 + 1.7212824826428488*α0807
  α1107 = -12.822511612651839 + 18.41810566087623*α0807
  α1201 = 0.2478349764611783 - 0.06383934946543009*α0807
  α1204 = -4.593782880309185 + 1.046256005127882*α0807
  α1205 = -0.39566692537411896 + 0.22827403748244698*α0807
  α1206 = -3.0673550479691665 - 0.10347586863902129*α0807
  α1207 = 5.386688702227177 - 1.1072148245058775*α0807
  α1301 = 0.7332242174431163 - 0.5164807626867616*α0807
  α1304 = -10.196728938160977 + 8.464545832921925*α0807
  α1305 = -0.43865244706547707 + 1.846809999910238*α0807
  α1306 = 0.5693856884667226 - 0.8371528845746959*α0807
  α1307 = 10.52865228002416 - 8.957722185570706*α0807

  k1 = cache.fsalfirst
  f(k1, uprev, p, t) 
  @.. tmp = uprev+dt*α0201*k1
  f(k2, tmp, p, t + c2*dt)
  @.. tmp = uprev+dt*(α0301*k1+α0302*k2)
  f(k3, tmp, p, t + c3*dt)
  @.. tmp = uprev+dt*(α0401*k1+α0403*k3)
  f(k4, tmp, p, t + c4*dt)
  @.. tmp = uprev+dt*(α0501*k1+α0503*k3+α0504*k4)
  f(k5, tmp, p, t + c5*dt)
  @.. tmp = uprev+dt*(α0601*k1+α0604*k4+α0605*k5)
  f(k6, tmp, p, t + c6*dt)
  @.. tmp = uprev+dt*(α0701*k1+α0704*k4+α0705*k5+α0706*k6)
  f(k7, tmp, p, t + c7*dt)
  @.. tmp = uprev+dt*(α0801*k1+α0804*k4+α0805*k5+α0806*k6+α0807*k7)
  f(k8, tmp, p, t + c8*dt)
  @.. tmp = uprev+dt*(α0901*k1+α0904*k4+α0905*k5+α0906*k6+α0907*k7+α0908*k8)
  f(k9, tmp, p, t + c9*dt)
  @.. tmp = uprev+dt*(α1001*k1+α1004*k4+α1005*k5+α1006*k6+α1007*k7+α1008*k8+α1009*k9)
  f(k10, tmp, p, t + c10*dt)
  @.. tmp = uprev+dt*(α1101*k1+α1104*k4+α1105*k5+α1106*k6+α1107*k7+α1108*k8+α1109*k9+α1110*k10)
  f(k11, tmp, p, t + c11*dt)
  @.. tmp = uprev+dt*(α1201*k1+α1204*k4+α1205*k5+α1206*k6+α1207*k7+α1208*k8+α1209*k9+α1210*k10+α1211*k11)
  f(k12, tmp, p, t+c12*dt)
  @.. tmp = uprev+dt*(α1301*k1+α1304*k4+α1305*k5+α1306*k6+α1307*k7+α1308*k8+α1309*k9+α1310*k10+α1311*k11)
  f(k13, tmp, p, t+c13*dt)
  @.. u = uprev + dt*(β1*k1+β6*k6+β7*k7+β8*k8+β9*k9+β10*k10+β11*k11+β12*k12+β13*k13)
  integrator.destats.nf += 13
  if integrator.opts.adaptive
    @.. utilde = dt*(β1tilde*k1 + β6tilde*k6 + β7tilde*k7 + β8tilde*k8 + β9tilde*k9 + β10tilde*k10 + β11tilde*k11 + β12tilde*k12)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
  integrator.destats.nf += 1
  return nothing
end
