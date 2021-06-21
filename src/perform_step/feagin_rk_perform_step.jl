function initialize!(integrator, cache::Feagin10ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Feagin10ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = cache
  k1  = integrator.fsalfirst
  a = dt*a0100
  k2  = f(uprev  + a*k1, p, t + c1*dt)
  k3  = f(uprev + dt*(a0200*k1 + a0201*k2), p, t + c2*dt)
  k4  = f(uprev  + dt*(a0300*k1              + a0302*k3), p, t + c3*dt)
  k5  = f(uprev  + dt*(a0400*k1              + a0402*k3 + a0403*k4), p, t + c4*dt)
  k6  = f(uprev  + dt*(a0500*k1                           + a0503*k4 + a0504*k5), p, t + c5*dt)
  k7  = f(uprev  + dt*(a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6), p, t + c6*dt)
  k8  = f(uprev  + dt*(a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7), p, t + c7*dt)
  k9  = f(uprev  + dt*(a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8), p, t + c8*dt)
  k10 = f(uprev  + dt*(a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9), p, t + c9*dt)
  k11 = f(uprev + dt*(a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10), p, t + c10*dt)
  k12 = f(uprev + dt*(a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11), p, t + c11*dt)
  k13 = f(uprev + dt*(a1200*k1                           + a1203*k4 + a1204*k5 + a1205*k6 + a1206*k7 + a1207*k8 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12), p, t + c12*dt)
  k14 = f(uprev + dt*(a1300*k1              + a1302*k3 + a1303*k4              + a1305*k6 + a1306*k7 + a1307*k8 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13), p, t + c13*dt)
  k15 = f(uprev + dt*(a1400*k1 + a1401*k2                           + a1404*k5              + a1406*k7 +                                                                     a1412*k13 + a1413*k14), p, t + c14*dt)
  k16 = f(uprev + dt*(a1500*k1              + a1502*k3                                                                                                                                                     + a1514*k15), p, t + c15*dt)
  k17 = f(uprev + dt*(a1600*k1 + a1601*k2 + a1602*k3              + a1604*k5 + a1605*k6 + a1606*k7 + a1607*k8 + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16), p, t + c16*dt)
  integrator.destats.nf += 16
  u = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b5*k5 + b7*k7 + b9*k9 + b10*k10 + b11*k11 + b12*k12 + b13*k13 + b14*k14 + b15*k15 + b16*k16 + b17*k17)
  if integrator.opts.adaptive
    utilde = @.. dt*(k2 - k16) * adaptiveConst
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  k = f(u, p, t+dt) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Feagin10Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

#=
@muladd function perform_step!(integrator, cache::Feagin10Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,uprev,k = cache
  k1 = cache.fsalfirst
  a =  dt*a0100
  @.. tmp = uprev + a*k1
  f(k2, tmp, p, t + c1*dt)
  @.. tmp = uprev + dt*(a0200*k1 + a0201*k2)
  f(k3, tmp, p, t + c2*dt )
  @.. tmp = uprev + dt*(a0300*k1 + a0302*k3)
  f(k4, tmp, p, t + c3*dt)
  @.. tmp = uprev + dt*(a0400*k1 + a0402*k3 + a0403*k4)
  f(k5, tmp, p, t + c4*dt)
  @.. tmp = uprev + dt*(a0500*k1 + a0503*k4 + a0504*k5)
  f(k6, tmp, p, t + c5*dt)
  @.. tmp = uprev + dt*(a0600*k1 + a0603*k4 + a0604*k5 + a0605*k6)
  f(k7, tmp, p, t + c6*dt)
  @.. tmp = uprev + dt*(a0700*k1 + a0704*k5 + a0705*k6 + a0706*k7)
  f(k8, tmp, p, t + c7*dt)
  @.. tmp = uprev + dt*(a0800*k1 + a0805*k6 + a0806*k7 + a0807*k8)
  f(k9, tmp, p, t + c8*dt)
  @.. tmp = uprev + dt*(a0900*k1 + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9)
  f(k10, tmp, p, t + c9*dt)
  @.. tmp = uprev + dt*(a1000*k1 + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10)
  f(k11, tmp, p, t + c10*dt)
  @.. tmp = uprev + dt*(a1100*k1 + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11)
  f(k12, tmp, p, t + c11*dt)
  @..  tmp = uprev + dt*(a1200*k1 + a1203*k4 + a1204*k5 + a1205*k6 + a1206*k7 + a1207*k8 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12)
  f(k13, tmp, p, t + c12*dt)
  @.. tmp = uprev + dt*(a1300*k1 + a1302*k3 + a1303*k4 + a1305*k6 + a1306*k7 + a1307*k8 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13)
  f(k14, tmp, p, t + c13*dt)
  @.. tmp = uprev + dt*(a1400*k1 + a1401*k2 + a1404*k5 + a1406*k7 + a1412*k13 + a1413*k14)
  f(k15, tmp, p, t + c14*dt)
  @.. tmp = uprev + dt*(a1500*k1 + a1502*k3 + a1514*k15)
  f(k16, tmp, p, t + c15*dt)
  @.. tmp = uprev + dt*(a1600*k1 + a1601*k2 + a1602*k3 + a1604*k5 + a1605*k6 + a1606*k7 + a1607*k8 + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16)
  f(k17, tmp, p, t + c16*dt)
  @.. u = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b5*k5 + b7*k7 + b9*k9 + b10*k10 + b11*k11 + b12*k12 + b13*k13 + b14*k14 + b15*k15 + b16*k16 + b17*k17)
  if integrator.opts.adaptive
    @.. tmp = dt*(k2 - k16) * adaptiveConst
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(integrator.fsallast,u,p,t+dt) # For the interpolation, needs k at the updated point
end
=#

@muladd function perform_step!(integrator, cache::Feagin10Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,uprev,k = cache
  k1 = cache.fsalfirst
  a =  dt*a0100
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + a*k1[i]
  end
  f(k2,tmp,p,t + c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0200*k1[i] + a0201*k2[i])
  end
  f(k3, tmp, p, t + c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0300*k1[i] + a0302*k3[i])
  end
  f(k4, tmp, p, t + c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0400*k1[i] + a0402*k3[i] + a0403*k4[i])
  end
  f(k5, tmp, p, t + c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0500*k1[i] + a0503*k4[i] + a0504*k5[i])
  end
  f(k6, tmp, p, t + c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0600*k1[i] + a0603*k4[i] + a0604*k5[i] + a0605*k6[i])
  end
  f(k7, tmp, p, t + c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0700*k1[i] + a0704*k5[i] + a0705*k6[i] + a0706*k7[i])
  end
  f(k8, tmp, p, t + c7*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0800*k1[i] + a0805*k6[i] + a0806*k7[i] + a0807*k8[i])
  end
  f(k9, tmp, p, t + c8*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0900*k1[i] + a0905*k6[i] + a0906*k7[i] + a0907*k8[i] + a0908*k9[i])
  end
  f(k10, tmp, p, t + c9*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1000*k1[i] + a1005*k6[i] + a1006*k7[i] + a1007*k8[i] + a1008*k9[i] + a1009*k10[i])
  end
  f(k11, tmp, p, t + c10*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1100*k1[i] + a1105*k6[i] + a1106*k7[i] + a1107*k8[i] + a1108*k9[i] + a1109*k10[i] + a1110*k11[i])
  end
  f(k12, tmp, p, t + c11*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1200*k1[i] + a1203*k4[i] + a1204*k5[i] + a1205*k6[i] + a1206*k7[i] + a1207*k8[i] + a1208*k9[i] + a1209*k10[i] + a1210*k11[i] + a1211*k12[i])
  end
  f(k13, tmp, p, t + c12*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1300*k1[i] + a1302*k3[i] + a1303*k4[i] + a1305*k6[i] + a1306*k7[i] + a1307*k8[i] + a1308*k9[i] + a1309*k10[i] + a1310*k11[i] + a1311*k12[i] + a1312*k13[i])
  end
  f(k14, tmp, p, t + c13*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1400*k1[i] + a1401*k2[i] + a1404*k5[i] + a1406*k7[i] + a1412*k13[i] + a1413*k14[i])
  end
  f(k15, tmp, p, t + c14*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1500*k1[i] + a1502*k3[i] + a1514*k15[i])
  end
  f(k16, tmp, p, t + c15*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1600*k1[i] + a1601*k2[i] + a1602*k3[i] + a1604*k5[i] + a1605*k6[i] + a1606*k7[i] + a1607*k8[i] + a1608*k9[i] + a1609*k10[i] + a1610*k11[i] + a1611*k12[i] + a1612*k13[i] + a1613*k14[i] + a1614*k15[i] + a1615*k16[i])
  end
  f(k17, tmp, p, t + c16*dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i] + b7*k7[i] + b9*k9[i] + b10*k10[i] + b11*k11[i] + b12*k12[i] + b13*k13[i] + b14*k14[i] + b15*k15[i] + b16*k16[i] + b17*k17[i])
  end
  integrator.destats.nf += 16
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = dt*(k2[i] - k16[i]) * adaptiveConst
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(integrator.fsallast,u,p,t+dt) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function initialize!(integrator, cache::Feagin12ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Feagin12ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = cache
  k1   = integrator.fsalfirst
  a =  dt*a0100
  k2  = f(uprev  + a*k1, p, t + c1*dt)
  k3  = f(uprev + dt*(a0200*k1 + a0201*k2), p, t + c2*dt)
  k4  = f(uprev  + dt*(a0300*k1              + a0302*k3), p, t + c3*dt)
  k5  = f(uprev  + dt*(a0400*k1              + a0402*k3 + a0403*k4), p, t + c4*dt)
  k6  = f(uprev  + dt*(a0500*k1                           + a0503*k4 + a0504*k5), p, t + c5*dt)
  k7  = f(uprev  + dt*(a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6), p, t + c6*dt)
  k8  = f(uprev  + dt*(a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7), p, t + c7*dt)
  k9  = f(uprev  + dt*(a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8), p, t + c8*dt)
  k10 = f(uprev  + dt*(a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9), p, t + c9*dt)
  k11 = f(uprev + dt*(a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10), p, t + c10*dt)
  k12 = f(uprev + dt*(a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11), p, t + c11*dt)
  k13 = f(uprev + dt*(a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12), p, t + c12*dt)
  k14 = f(uprev + dt*(a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13), p, t + c13*dt)
  k15 = f(uprev + dt*(a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14), p, t + c14*dt)
  k16 = f(uprev + dt*(a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15), p, t + c15*dt)
  k17 = f(uprev + dt*((a1600*k1                                                                                            + a1608*k9 + a1609*k10) + (a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14) + (a1614*k15 + a1615*k16)), p, t + c16*dt)
  k18 = f(uprev + dt*((a1700*k1                                                     + a1705*k6 + a1706*k7) + (a1707*k8 + a1708*k9 + a1709*k10 + a1710*k11) + (a1711*k12 + a1712*k13 + a1713*k14 + a1714*k15) + (a1715*k16 + a1716*k17)), p, t + c17*dt)
  k19 = f(uprev + dt*((a1800*k1                                                     + a1805*k6 + a1806*k7) + (a1807*k8 + a1808*k9 + a1809*k10 + a1810*k11) + (a1811*k12 + a1812*k13 + a1813*k14 + a1814*k15) + (a1815*k16 + a1816*k17 + a1817*k18)), p, t + c18*dt)
  k20 = f(uprev + dt*((a1900*k1                                        + a1904*k5 + a1905*k6) + (a1906*k7              + a1908*k9 + a1909*k10 + a1910*k11) + (a1911*k12 + a1912*k13 + a1913*k14 + a1914*k15) + (a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19)), p, t + c19*dt)
  k21 = f(uprev + dt*((a2000*k1                           + a2003*k4 + a2004*k5) + (a2005*k6              + a2007*k8              + a2009*k10 + a2010*k11)                                                                                     + (a2017*k18 + a2018*k19 + a2019*k20)), p, t + c20*dt)
  k22 = f(uprev + dt*((a2100*k1              + a2102*k3 + a2103*k4)                           + (a2106*k7 + a2107*k8              + a2109*k10 + a2110*k11)                                                                                     + (a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21)), p, t + c21*dt)
  k23 = f(uprev + dt*((a2200*k1 + a2201*k2                           + a2204*k5)              + (a2206*k7                                                                                                                                                                                     + a2220*k21 + a2221*k22)), p, t + c22*dt)
  k24 = f(uprev + dt*(a2300*k1              + a2302*k3                                                                                                                                                                                                                                                                     + a2322*k23), p, t + c23*dt)
  k25 = f(uprev + dt*((a2400*k1 + a2401*k2 + a2402*k3)              + (a2404*k5              + a2406*k7 + a2407*k8 + a2408*k9) + (a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13) + (a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17) + (a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21) + (a2421*k22 + a2422*k23 + a2423*k24)), p, t + c24*dt)
  integrator.destats.nf += 24
  u = uprev + dt*((b1*k1 + b2*k2 + b3*k3 + b5*k5) + (b7*k7 + b8*k8 + b10*k10 + b11*k11) + (b13*k13 + b14*k14 + b15*k15 + b16*k16) + (b17*k17 + b18*k18 + b19*k19 + b20*k20) + (b21*k21 + b22*k22 + b23*k23) + (b24*k24 + b25*k25))
  k = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.fsallast = k
  if integrator.opts.adaptive
    utilde = @.. dt*(k2 - k24) * adaptiveConst
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Feagin12Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

#=
@muladd function perform_step!(integrator, cache::Feagin12Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,uprev,k = cache
  k1 = cache.fsalfirst
  a = dt*a0100
  @.. tmp = uprev + a*k1
  f(k2, tmp, p, t + c1*dt)
  @.. tmp = uprev + dt*(a0200*k1 + a0201*k2)
  f(k3, tmp, p, t + c2*dt )
  @.. tmp = uprev + dt*(a0300*k1 + a0302*k3)
  f(k4, tmp, p, t + c3*dt)
  @.. tmp = uprev + dt*(a0400*k1 + a0402*k3 + a0403*k4)
  f(k5, tmp, p, t + c4*dt)
  @.. tmp = uprev + dt*(a0500*k1 + a0503*k4 + a0504*k5)
  f(k6, tmp, p, t + c5*dt)
  @.. tmp = uprev + dt*(a0600*k1 + a0603*k4 + a0604*k5 + a0605*k6)
  f(k7, tmp, p, t + c6*dt)
  @.. tmp = uprev + dt*(a0700*k1 + a0704*k5 + a0705*k6 + a0706*k7)
  f(k8, tmp, p, t + c7*dt)
  @.. tmp = uprev + dt*(a0800*k1 + a0805*k6 + a0806*k7 + a0807*k8)
  f(k9, tmp, p, t + c8*dt)
  @.. tmp = uprev + dt*(a0900*k1 + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9)
  f(k10, tmp, p, t + c9*dt)
  @.. tmp = uprev + dt*(a1000*k1 + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10)
  f(k11, tmp, p, t + c10*dt)
  @.. tmp = uprev + dt*(a1100*k1 + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11)
  f(k12, tmp, p, t + c11*dt)
  @.. tmp = uprev + dt*(a1200*k1 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12)
  f(k13, tmp, p, t + c12*dt)
  @.. tmp = uprev + dt*(a1300*k1 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13)
  f(k14, tmp, p, t + c13*dt)
  @.. tmp = uprev + dt*(a1400*k1 + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14)
  f(k15, tmp, p, t + c14*dt)
  @.. tmp = uprev + dt*(a1500*k1 + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15)
  f(k16, tmp, p, t + c15*dt)
  @.. tmp = uprev + dt*((a1600*k1 + a1608*k9 + a1609*k10) + (a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14) + (a1614*k15 + a1615*k16))
  f(k17, tmp, p, t + c16*dt)
  @.. tmp = uprev + dt*((a1700*k1 + a1705*k6 + a1706*k7) + (a1707*k8 + a1708*k9 + a1709*k10 + a1710*k11) + (a1711*k12 + a1712*k13 + a1713*k14 + a1714*k15) + (a1715*k16 + a1716*k17))
  f(k18, tmp, p, t + c17*dt)
  @.. tmp = uprev + dt*((a1800*k1 + a1805*k6 + a1806*k7) + (a1807*k8 + a1808*k9 + a1809*k10 + a1810*k11) + (a1811*k12 + a1812*k13 + a1813*k14 + a1814*k15) + (a1815*k16 + a1816*k17 + a1817*k18))
  f(k19, tmp, p, t + c18*dt)
  @.. tmp = uprev + dt*((a1900*k1 + a1904*k5 + a1905*k6) + (a1906*k7 + a1908*k9 + a1909*k10 + a1910*k11) + (a1911*k12 + a1912*k13 + a1913*k14 + a1914*k15) + (a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19))
  f(k20, tmp, p, t + c19*dt)
  @.. tmp = uprev + dt*((a2000*k1 + a2003*k4 + a2004*k5) + (a2005*k6 + a2007*k8 + a2009*k10 + a2010*k11) + (a2017*k18 + a2018*k19 + a2019*k20))
  f(k21, tmp, p, t + c20*dt)
  @.. tmp = uprev + dt*((a2100*k1 + a2102*k3 + a2103*k4) + (a2106*k7 + a2107*k8 + a2109*k10 + a2110*k11) + (a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21))
  f(k22, tmp, p, t + c21*dt)
  @.. tmp = uprev + dt*((a2200*k1 + a2201*k2 + a2204*k5) + (a2206*k7 + a2220*k21 + a2221*k22))
  f(k23, tmp, p, t + c22*dt)
  @.. tmp = uprev + dt*(a2300*k1 + a2302*k3 + a2322*k23)
  f(k24, tmp, p, t + c23*dt)
  @.. tmp = uprev + dt*((a2400*k1 + a2401*k2 + a2402*k3) + (a2404*k5 + a2406*k7 + a2407*k8 + a2408*k9) + (a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13) + (a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17) + (a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21) + (a2421*k22 + a2422*k23 + a2423*k24))
  f(k25, tmp, p, t + c24*dt)
  @.. u = uprev + dt*((b1*k1 + b2*k2 + b3*k3 + b5*k5) + (b7*k7 + b8*k8 + b10*k10 + b11*k11) + (b13*k13 + b14*k14 + b15*k15 + b16*k16) + (b17*k17 + b18*k18 + b19*k19 + b20*k20) + (b21*k21 + b22*k22 + b23*k23) + (b24*k24 + b25*k25))
  if integrator.opts.adaptive
    @.. tmp = dt*(k2 - k24) * adaptiveConst
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
end
=#

@muladd function perform_step!(integrator, cache::Feagin12Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,uprev,k = cache
  k1 = cache.fsalfirst
  a = dt*a0100
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + a*k1[i]
  end
  f(k2, tmp, p, t + c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0200*k1[i] + a0201*k2[i])
  end
  f(k3, tmp, p, t + c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0300*k1[i] + a0302*k3[i])
  end
  f(k4, tmp, p, t + c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0400*k1[i] + a0402*k3[i] + a0403*k4[i])
  end
  f(k5, tmp, p, t + c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0500*k1[i] + a0503*k4[i] + a0504*k5[i])
  end
  f(k6, tmp, p, t + c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0600*k1[i] + a0603*k4[i] + a0604*k5[i] + a0605*k6[i])
  end
  f(k7, tmp, p, t + c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0700*k1[i] + a0704*k5[i] + a0705*k6[i] + a0706*k7[i])
  end
  f(k8, tmp, p, t + c7*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0800*k1[i] + a0805*k6[i] + a0806*k7[i] + a0807*k8[i])
  end
  f(k9, tmp, p, t + c8*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0900*k1[i] + a0905*k6[i] + a0906*k7[i] + a0907*k8[i] + a0908*k9[i])
  end
  f(k10, tmp, p, t + c9*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1000*k1[i] + a1005*k6[i] + a1006*k7[i] + a1007*k8[i] + a1008*k9[i] + a1009*k10[i])
  end
  f(k11, tmp, p, t + c10*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1100*k1[i] + a1105*k6[i] + a1106*k7[i] + a1107*k8[i] + a1108*k9[i] + a1109*k10[i] + a1110*k11[i])
  end
  f(k12, tmp, p, t + c11*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1200*k1[i] + a1208*k9[i] + a1209*k10[i] + a1210*k11[i] + a1211*k12[i])
  end
  f(k13, tmp, p, t + c12*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1300*k1[i] + a1308*k9[i] + a1309*k10[i] + a1310*k11[i] + a1311*k12[i] + a1312*k13[i])
  end
  f(k14, tmp, p, t + c13*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1400*k1[i] + a1408*k9[i] + a1409*k10[i] + a1410*k11[i] + a1411*k12[i] + a1412*k13[i] + a1413*k14[i])
  end
  f(k15, tmp, p, t + c14*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1500*k1[i] + a1508*k9[i] + a1509*k10[i] + a1510*k11[i] + a1511*k12[i] + a1512*k13[i] + a1513*k14[i] + a1514*k15[i])
  end
  f(k16, tmp, p, t + c15*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a1600*k1[i] + a1608*k9[i] + a1609*k10[i]) + (a1610*k11[i] + a1611*k12[i] + a1612*k13[i] + a1613*k14[i]) + (a1614*k15[i] + a1615*k16[i]))
  end
  f(k17, tmp, p, t + c16*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a1700*k1[i] + a1705*k6[i] + a1706*k7[i]) + (a1707*k8[i] + a1708*k9[i] + a1709*k10[i] + a1710*k11[i]) + (a1711*k12[i] + a1712*k13[i] + a1713*k14[i] + a1714*k15[i]) + (a1715*k16[i] + a1716*k17[i]))
  end
  f(k18, tmp, p, t + c17*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a1800*k1[i] + a1805*k6[i] + a1806*k7[i]) + (a1807*k8[i] + a1808*k9[i] + a1809*k10[i] + a1810*k11[i]) + (a1811*k12[i] + a1812*k13[i] + a1813*k14[i] + a1814*k15[i]) + (a1815*k16[i] + a1816*k17[i] + a1817*k18[i]))
  end
  f(k19, tmp, p, t + c18*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a1900*k1[i] + a1904*k5[i] + a1905*k6[i]) + (a1906*k7[i] + a1908*k9[i] + a1909*k10[i] + a1910*k11[i]) + (a1911*k12[i] + a1912*k13[i] + a1913*k14[i] + a1914*k15[i]) + (a1915*k16[i] + a1916*k17[i] + a1917*k18[i] + a1918*k19[i]))
  end
  f(k20, tmp, p, t + c19*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a2000*k1[i] + a2003*k4[i] + a2004*k5[i]) + (a2005*k6[i] + a2007*k8[i] + a2009*k10[i] + a2010*k11[i]) + (a2017*k18[i] + a2018*k19[i] + a2019*k20[i]))
  end
  f(k21, tmp, p, t + c20*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a2100*k1[i] + a2102*k3[i] + a2103*k4[i]) + (a2106*k7[i] + a2107*k8[i] + a2109*k10[i] + a2110*k11[i]) + (a2117*k18[i] + a2118*k19[i] + a2119*k20[i] + a2120*k21[i]))
  end
  f(k22, tmp, p, t + c21*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a2200*k1[i] + a2201*k2[i] + a2204*k5[i]) + (a2206*k7[i] + a2220*k21[i] + a2221*k22[i]))
  end
  f(k23, tmp, p, t + c22*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2300*k1[i] + a2302*k3[i] + a2322*k23[i])
  end
  f(k24, tmp, p, t + c23*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*((a2400*k1[i] + a2401*k2[i] + a2402*k3[i]) + (a2404*k5[i] + a2406*k7[i] + a2407*k8[i] + a2408*k9[i]) + (a2409*k10[i] + a2410*k11[i] + a2411*k12[i] + a2412*k13[i]) + (a2413*k14[i] + a2414*k15[i] + a2415*k16[i] + a2416*k17[i]) + (a2417*k18[i] + a2418*k19[i] + a2419*k20[i] + a2420*k21[i]) + (a2421*k22[i] + a2422*k23[i] + a2423*k24[i]))
  end
  f(k25, tmp, p, t + c24*dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + dt*((b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b13*k13[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b17*k17[i] + b18*k18[i] + b19*k19[i] + b20*k20[i]) + (b21*k21[i] + b22*k22[i] + b23*k23[i]) + (b24*k24[i] + b25*k25[i]))
  end
  integrator.destats.nf += 24
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = dt*(k2[i] - k24[i]) * adaptiveConst
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end

function initialize!(integrator,cache::Feagin14ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Feagin14ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = cache
  k1  = integrator.fsalfirst
  a = dt*a0100
  k2  = f(uprev  + a*k1, p, t + c1*dt)
  k3  = f(uprev + dt*(a0200*k1 + a0201*k2), p, t + c2*dt)
  k4  = f(uprev  + dt*(a0300*k1              + a0302*k3), p, t + c3*dt)
  k5  = f(uprev  + dt*(a0400*k1              + a0402*k3 + a0403*k4), p, t + c4*dt)
  k6  = f(uprev  + dt*(a0500*k1                           + a0503*k4 + a0504*k5), p, t + c5*dt)
  k7  = f(uprev  + dt*(a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6), p, t + c6*dt)
  k8  = f(uprev  + dt*(a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7), p, t + c7*dt)
  k9  = f(uprev  + dt*(a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8), p, t + c8*dt)
  k10 = f(uprev  + dt*(a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9), p, t + c9*dt)
  k11 = f(uprev + dt*(a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10), p, t + c10*dt)
  k12 = f(uprev + dt*(a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11), p, t + c11*dt)
  k13 = f(uprev + dt*(a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12), p, t + c12*dt)
  k14 = f(uprev + dt*(a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13), p, t + c13*dt)
  k15 = f(uprev + dt*(a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14), p, t + c14*dt)
  k16 = f(uprev + dt*(a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15), p, t + c15*dt)
  k17 = f(uprev + dt*(a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16), p, t + c16*dt)
  k18 = f(uprev + dt*(a1700*k1                                                                                                                                                   + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17), p, t + c17*dt)
  k19 = f(uprev + dt*(a1800*k1                                                                                                                                                   + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18), p, t + c18*dt)
  k20 = f(uprev + dt*(a1900*k1                                                                                                                                                   + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19), p, t + c19*dt)
  k21 = f(uprev + dt*(a2000*k1                                                                                                                                                   + a2012*k13 + a2013*k14 + a2014*k15 + a2015*k16 + a2016*k17 + a2017*k18 + a2018*k19 + a2019*k20), p, t + c20*dt)
  k22 = f(uprev + dt*(a2100*k1                                                                                                                                                   + a2112*k13 + a2113*k14 + a2114*k15 + a2115*k16 + a2116*k17 + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21), p, t + c21*dt)
  k23 = f(uprev + dt*(a2200*k1                                                                                                                                                   + a2212*k13 + a2213*k14 + a2214*k15 + a2215*k16 + a2216*k17 + a2217*k18 + a2218*k19 + a2219*k20 + a2220*k21 + a2221*k22), p, t + c22*dt)
  k24 = f(uprev + dt*(a2300*k1                                                                                            + a2308*k9 + a2309*k10 + a2310*k11 + a2311*k12 + a2312*k13 + a2313*k14 + a2314*k15 + a2315*k16 + a2316*k17 + a2317*k18 + a2318*k19 + a2319*k20 + a2320*k21 + a2321*k22 + a2322*k23), p, t + c23*dt)
  k25 = f(uprev + dt*(a2400*k1                                                                                            + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24), p, t + c24*dt)
  k26 = f(uprev + dt*(a2500*k1                                                                                            + a2508*k9 + a2509*k10 + a2510*k11 + a2511*k12 + a2512*k13 + a2513*k14 + a2514*k15 + a2515*k16 + a2516*k17 + a2517*k18 + a2518*k19 + a2519*k20 + a2520*k21 + a2521*k22 + a2522*k23 + a2523*k24 + a2524*k25), p, t + c25*dt)
  k27 = f(uprev + dt*(a2600*k1                                                     + a2605*k6 + a2606*k7 + a2607*k8 + a2608*k9 + a2609*k10 + a2610*k11               + a2612*k13 + a2613*k14 + a2614*k15 + a2615*k16 + a2616*k17 + a2617*k18 + a2618*k19 + a2619*k20 + a2620*k21 + a2621*k22 + a2622*k23 + a2623*k24 + a2624*k25 + a2625*k26), p, t + c26*dt)
  k28 = f(uprev + dt*(a2700*k1                                                     + a2705*k6 + a2706*k7 + a2707*k8 + a2708*k9 + a2709*k10               + a2711*k12 + a2712*k13 + a2713*k14 + a2714*k15 + a2715*k16 + a2716*k17 + a2717*k18 + a2718*k19 + a2719*k20 + a2720*k21 + a2721*k22 + a2722*k23 + a2723*k24 + a2724*k25 + a2725*k26 + a2726*k27), p, t + c27*dt)
  k29 = f(uprev + dt*(a2800*k1                                                     + a2805*k6 + a2806*k7 + a2807*k8 + a2808*k9               + a2810*k11 + a2811*k12               + a2813*k14 + a2814*k15 + a2815*k16                                                                                                   + a2823*k24 + a2824*k25 + a2825*k26 + a2826*k27 + a2827*k28), p, t + c28*dt)
  k30 = f(uprev + dt*(a2900*k1                                        + a2904*k5 + a2905*k6 + a2906*k7                           + a2909*k10 + a2910*k11 + a2911*k12               + a2913*k14 + a2914*k15 + a2915*k16                                                                                                   + a2923*k24 + a2924*k25 + a2925*k26 + a2926*k27 + a2927*k28 + a2928*k29), p, t + c29*dt)
  k31 = f(uprev + dt*(a3000*k1                           + a3003*k4 + a3004*k5 + a3005*k6              + a3007*k8              + a3009*k10 + a3010*k11                             + a3013*k14 + a3014*k15 + a3015*k16                                                                                                   + a3023*k24 + a3024*k25 + a3025*k26               + a3027*k28 + a3028*k29 + a3029*k30), p, t + c30*dt)
  k32 = f(uprev + dt*(a3100*k1              + a3102*k3 + a3103*k4                           + a3106*k7 + a3107*k8              + a3109*k10 + a3110*k11                             + a3113*k14 + a3114*k15 + a3115*k16                                                                                                   + a3123*k24 + a3124*k25 + a3125*k26               + a3127*k28 + a3128*k29 + a3129*k30 + a3130*k31), p, t + c31*dt)
  k33 = f(uprev + dt*(a3200*k1 + a3201*k2                           + a3204*k5              + a3206*k7                                                                                                                                                                                                                                                                                                                                 + a3230*k31 + a3231*k32), p, t + c32*dt)
  k34 = f(uprev + dt*(a3300*k1              + a3302*k3                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k33), p, t + c33*dt)
  k35 = f(uprev + dt*(a3400*k1 + a3401*k2 + a3402*k3              + a3404*k5              + a3406*k7 + a3407*k8              + a3409*k10 + a3410*k11 + a3411*k12 + a3412*k13 + a3413*k14 + a3414*k15 + a3415*k16 + a3416*k17 + a3417*k18 + a3418*k19 + a3419*k20 + a3420*k21 + a3421*k22 + a3422*k23 + a3423*k24 + a3424*k25 + a3425*k26 + a3426*k27 + a3427*k28 + a3428*k29 + a3429*k30 + a3430*k31 + a3431*k32 + a3432*k33 + a3433*k34), p, t + c34*dt)
  integrator.destats.nf += 34
  u = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b5*k5 + b7*k7 + b8*k8 + b10*k10 + b11*k11 + b12*k12 + b14*k14 + b15*k15 + b16*k16 + b18*k18 + b19*k19 + b20*k20 + b21*k21 + b22*k22 + b23*k23 + b24*k24 + b25*k25 + b26*k26 + b27*k27 + b28*k28 + b29*k29 + b30*k30 + b31*k31 + b32*k32 + b33*k33 + b34*k34 + b35*k35)
  if integrator.opts.adaptive
    utilde = @.. dt*(k2 - k34) * adaptiveConst
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  k = f(u, p, t+dt) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Feagin14Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

#=
@muladd function perform_step!(integrator, cache::Feagin14Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,tmp,atmp,uprev,k = cache
  k1 = cache.fsalfirst
  f(k1, uprev, p, t)
  a = dt*a0100
  @.. tmp = uprev + a*k1
  f(k2, tmp, p, t + c1*dt)
  @.. tmp = uprev + dt*(a0200*k1 + a0201*k2)
  f(k3, tmp, p, t + c2*dt )
  @.. tmp = uprev + dt*(a0300*k1 + a0302*k3)
  f(k4, tmp, p, t + c3*dt)
  @.. tmp = uprev + dt*(a0400*k1 + a0402*k3 + a0403*k4)
  f(k5, tmp, p, t + c4*dt)
  @.. tmp = uprev + dt*(a0500*k1 + a0503*k4 + a0504*k5)
  f(k6, tmp, p, t + c5*dt)
  @.. tmp = uprev + dt*(a0600*k1 + a0603*k4 + a0604*k5 + a0605*k6)
  f(k7, tmp, p, t + c6*dt)
  @.. tmp = uprev + dt*(a0700*k1 + a0704*k5 + a0705*k6 + a0706*k7)
  f(k8, tmp, p, t + c7*dt)
  @.. tmp = uprev + dt*(a0800*k1 + a0805*k6 + a0806*k7 + a0807*k8)
  f(k9, tmp, p, t + c8*dt)
  @.. tmp = uprev + dt*(a0900*k1 + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9)
  f(k10, tmp, p, t + c9*dt)
  @.. tmp = uprev + dt*(a1000*k1 + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10)
  f(k11, tmp, p, t + c10*dt)
  @.. tmp = uprev + dt*(a1100*k1 + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11)
  f(k12, tmp, p, t + c11*dt)
  @.. tmp = uprev + dt*(a1200*k1 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12)
  f(k13, tmp, p, t + c12*dt)
  @.. tmp = uprev + dt*(a1300*k1 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13)
  f(k14, tmp, p, t + c13*dt)
  @.. tmp = uprev + dt*(a1400*k1 + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14)
  f(k15, tmp, p, t + c14*dt)
  @.. tmp = uprev + dt*(a1500*k1 + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15)
  f(k16, tmp, p, t + c15*dt)
  @.. tmp = uprev + dt*(a1600*k1 + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16)
  f(k17, tmp, p, t + c16*dt)
  @.. tmp = uprev + dt*(a1700*k1 + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17)
  f(k18, tmp, p, t + c17*dt)
  @.. tmp = uprev + dt*(a1800*k1 + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18)
  f(k19, tmp, p, t + c18*dt)
  @.. tmp = uprev + dt*(a1900*k1 + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19)
  f(k20, tmp, p, t + c19*dt)
  @.. tmp = uprev + dt*(a2000*k1 + a2012*k13 + a2013*k14 + a2014*k15 + a2015*k16 + a2016*k17 + a2017*k18 + a2018*k19 + a2019*k20)
  f(k21, tmp, p, t + c20*dt)
  @.. tmp = uprev + dt*(a2100*k1 + a2112*k13 + a2113*k14 + a2114*k15 + a2115*k16 + a2116*k17 + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21)
  f(k22, tmp, p, t + c21*dt)
  @.. tmp = uprev + dt*(a2200*k1 + a2212*k13 + a2213*k14 + a2214*k15 + a2215*k16 + a2216*k17 + a2217*k18 + a2218*k19 + a2219*k20 + a2220*k21 + a2221*k22)
  f(k23, tmp, p, t + c22*dt)
  @.. tmp = uprev + dt*(a2300*k1 + a2308*k9 + a2309*k10 + a2310*k11 + a2311*k12 + a2312*k13 + a2313*k14 + a2314*k15 + a2315*k16 + a2316*k17 + a2317*k18 + a2318*k19 + a2319*k20 + a2320*k21 + a2321*k22 + a2322*k23)
  f(k24, tmp, p, t + c23*dt)
  @.. tmp = uprev + dt*(a2400*k1 + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24)
  f(k25, tmp, p, t + c24*dt)
  @.. tmp = uprev + dt*(a2500*k1 + a2508*k9 + a2509*k10 + a2510*k11 + a2511*k12 + a2512*k13 + a2513*k14 + a2514*k15 + a2515*k16 + a2516*k17 + a2517*k18 + a2518*k19 + a2519*k20 + a2520*k21 + a2521*k22 + a2522*k23 + a2523*k24 + a2524*k25)
  f(k26, tmp, p, t + c25*dt)
  @.. tmp = uprev + dt*(a2600*k1 + a2605*k6 + a2606*k7 + a2607*k8 + a2608*k9 + a2609*k10 + a2610*k11 + a2612*k13 + a2613*k14 + a2614*k15 + a2615*k16 + a2616*k17 + a2617*k18 + a2618*k19 + a2619*k20 + a2620*k21 + a2621*k22 + a2622*k23 + a2623*k24 + a2624*k25 + a2625*k26)
  f(k27, tmp, p, t + c26*dt)
  @.. tmp = uprev + dt*(a2700*k1 + a2705*k6 + a2706*k7 + a2707*k8 + a2708*k9 + a2709*k10 + a2711*k12 + a2712*k13 + a2713*k14 + a2714*k15 + a2715*k16 + a2716*k17 + a2717*k18 + a2718*k19 + a2719*k20 + a2720*k21 + a2721*k22 + a2722*k23 + a2723*k24 + a2724*k25 + a2725*k26 + a2726*k27)
  f(k28, tmp, p, t + c27*dt)
  @.. tmp = uprev + dt*(a2800*k1 + a2805*k6 + a2806*k7 + a2807*k8 + a2808*k9 + a2810*k11 + a2811*k12 + a2813*k14 + a2814*k15 + a2815*k16 + a2823*k24 + a2824*k25 + a2825*k26 + a2826*k27 + a2827*k28)
  f(k29, tmp, p, t + c28*dt)
  @.. tmp = uprev + dt*(a2900*k1 + a2904*k5 + a2905*k6 + a2906*k7 + a2909*k10 + a2910*k11 + a2911*k12 + a2913*k14 + a2914*k15 + a2915*k16 + a2923*k24 + a2924*k25 + a2925*k26 + a2926*k27 + a2927*k28 + a2928*k29)
  f(k30, tmp, p, t + c29*dt)
  @.. tmp = uprev + dt*(a3000*k1 + a3003*k4 + a3004*k5 + a3005*k6 + a3007*k8 + a3009*k10 + a3010*k11 + a3013*k14 + a3014*k15 + a3015*k16 + a3023*k24 + a3024*k25 + a3025*k26 + a3027*k28 + a3028*k29 + a3029*k30)
  f(k31, tmp, p, t + c30*dt)
  @.. tmp = uprev + dt*(a3100*k1 + a3102*k3 + a3103*k4 + a3106*k7 + a3107*k8 + a3109*k10 + a3110*k11 + a3113*k14 + a3114*k15 + a3115*k16 + a3123*k24 + a3124*k25 + a3125*k26 + a3127*k28 + a3128*k29 + a3129*k30 + a3130*k31)
  f(k32, tmp, p, t + c31*dt)
  @.. tmp = uprev + dt*(a3200*k1 + a3201*k2 + a3204*k5 + a3206*k7 + a3230*k31 + a3231*k32)
  f(k33, tmp, p, t + c32*dt)
  @.. tmp = uprev + dt*(a3300*k1 + a3302*k3 + a3332*k33)
  f(k34, tmp, p, t + c33*dt)
  @.. tmp = uprev + dt*(a3400*k1 + a3401*k2 + a3402*k3 + a3404*k5 + a3406*k7 + a3407*k8 + a3409*k10 + a3410*k11 + a3411*k12 + a3412*k13 + a3413*k14 + a3414*k15 + a3415*k16 + a3416*k17 + a3417*k18 + a3418*k19 + a3419*k20 + a3420*k21 + a3421*k22 + a3422*k23 + a3423*k24 + a3424*k25 + a3425*k26 + a3426*k27 + a3427*k28 + a3428*k29 + a3429*k30 + a3430*k31 + a3431*k32 + a3432*k33 + a3433*k34)
  f(k35, tmp, p, t + c34*dt)
  @.. u = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b5*k5 + b7*k7 + b8*k8 + b10*k10 + b11*k11 + b12*k12 + b14*k14 + b15*k15 + b16*k16 + b18*k18 + b19*k19 + b20*k20 + b21*k21 + b22*k22 + b23*k23 + b24*k24 + b25*k25 + b26*k26 + b27*k27 + b28*k28 + b29*k29 + b30*k30 + b31*k31 + b32*k32 + b33*k33 + b34*k34 + b35*k35)
  if integrator.opts.adaptive
    @.. tmp = dt*(k2 - k34) * adaptiveConst
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(integrator.fsallast,u,p,t+dt) # For the interpolation, needs k at the updated point
end
=#

@muladd function perform_step!(integrator, cache::Feagin14Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = cache.tab
  @unpack k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,tmp,atmp,uprev,k = cache
  k1 = cache.fsalfirst
  a = dt*a0100
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + a*k1[i]
  end
  f(k2, tmp, p, t + c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0200*k1[i] + a0201*k2[i])
  end
  f(k3, tmp, p, t + c2*dt )
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0300*k1[i] + a0302*k3[i])
  end
  f(k4, tmp, p, t + c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0400*k1[i] + a0402*k3[i] + a0403*k4[i])
  end
  f(k5, tmp, p, t + c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0500*k1[i] + a0503*k4[i] + a0504*k5[i])
  end
  f(k6, tmp, p, t + c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0600*k1[i] + a0603*k4[i] + a0604*k5[i] + a0605*k6[i])
  end
  f(k7, tmp, p, t + c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0700*k1[i] + a0704*k5[i] + a0705*k6[i] + a0706*k7[i])
  end
  f(k8, tmp, p, t + c7*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0800*k1[i] + a0805*k6[i] + a0806*k7[i] + a0807*k8[i])
  end
  f(k9, tmp, p, t + c8*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a0900*k1[i] + a0905*k6[i] + a0906*k7[i] + a0907*k8[i] + a0908*k9[i])
  end
  f(k10, tmp, p, t + c9*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1000*k1[i] + a1005*k6[i] + a1006*k7[i] + a1007*k8[i] + a1008*k9[i] + a1009*k10[i])
  end
  f(k11, tmp, p, t + c10*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1100*k1[i] + a1105*k6[i] + a1106*k7[i] + a1107*k8[i] + a1108*k9[i] + a1109*k10[i] + a1110*k11[i])
  end
  f(k12, tmp, p, t + c11*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1200*k1[i] + a1208*k9[i] + a1209*k10[i] + a1210*k11[i] + a1211*k12[i])
  end
  f(k13, tmp, p, t + c12*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1300*k1[i] + a1308*k9[i] + a1309*k10[i] + a1310*k11[i] + a1311*k12[i] + a1312*k13[i])
  end
  f(k14, tmp, p, t + c13*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1400*k1[i] + a1408*k9[i] + a1409*k10[i] + a1410*k11[i] + a1411*k12[i] + a1412*k13[i] + a1413*k14[i])
  end
  f(k15, tmp, p, t + c14*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1500*k1[i] + a1508*k9[i] + a1509*k10[i] + a1510*k11[i] + a1511*k12[i] + a1512*k13[i] + a1513*k14[i] + a1514*k15[i])
  end
  f(k16, tmp, p, t + c15*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1600*k1[i] + a1608*k9[i] + a1609*k10[i] + a1610*k11[i] + a1611*k12[i] + a1612*k13[i] + a1613*k14[i] + a1614*k15[i] + a1615*k16[i])
  end
  f(k17, tmp, p, t + c16*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1700*k1[i] + a1712*k13[i] + a1713*k14[i] + a1714*k15[i] + a1715*k16[i] + a1716*k17[i])
  end
  f(k18, tmp, p, t + c17*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1800*k1[i] + a1812*k13[i] + a1813*k14[i] + a1814*k15[i] + a1815*k16[i] + a1816*k17[i] + a1817*k18[i])
  end
  f(k19, tmp, p, t + c18*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a1900*k1[i] + a1912*k13[i] + a1913*k14[i] + a1914*k15[i] + a1915*k16[i] + a1916*k17[i] + a1917*k18[i] + a1918*k19[i])
  end
  f(k20, tmp, p, t + c19*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2000*k1[i] + a2012*k13[i] + a2013*k14[i] + a2014*k15[i] + a2015*k16[i] + a2016*k17[i] + a2017*k18[i] + a2018*k19[i] + a2019*k20[i])
  end
  f(k21, tmp, p, t + c20*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2100*k1[i] + a2112*k13[i] + a2113*k14[i] + a2114*k15[i] + a2115*k16[i] + a2116*k17[i] + a2117*k18[i] + a2118*k19[i] + a2119*k20[i] + a2120*k21[i])
  end
  f(k22, tmp, p, t + c21*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2200*k1[i] + a2212*k13[i] + a2213*k14[i] + a2214*k15[i] + a2215*k16[i] + a2216*k17[i] + a2217*k18[i] + a2218*k19[i] + a2219*k20[i] + a2220*k21[i] + a2221*k22[i])
  end
  f(k23, tmp, p, t + c22*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2300*k1[i] + a2308*k9[i] + a2309*k10[i] + a2310*k11[i] + a2311*k12[i] + a2312*k13[i] + a2313*k14[i] + a2314*k15[i] + a2315*k16[i] + a2316*k17[i] + a2317*k18[i] + a2318*k19[i] + a2319*k20[i] + a2320*k21[i] + a2321*k22[i] + a2322*k23[i])
  end
  f(k24, tmp, p, t + c23*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2400*k1[i] + a2408*k9[i] + a2409*k10[i] + a2410*k11[i] + a2411*k12[i] + a2412*k13[i] + a2413*k14[i] + a2414*k15[i] + a2415*k16[i] + a2416*k17[i] + a2417*k18[i] + a2418*k19[i] + a2419*k20[i] + a2420*k21[i] + a2421*k22[i] + a2422*k23[i] + a2423*k24[i])
  end
  f(k25, tmp, p, t + c24*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2500*k1[i] + a2508*k9[i] + a2509*k10[i] + a2510*k11[i] + a2511*k12[i] + a2512*k13[i] + a2513*k14[i] + a2514*k15[i] + a2515*k16[i] + a2516*k17[i] + a2517*k18[i] + a2518*k19[i] + a2519*k20[i] + a2520*k21[i] + a2521*k22[i] + a2522*k23[i] + a2523*k24[i] + a2524*k25[i])
  end
  f(k26, tmp, p, t + c25*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2600*k1[i] + a2605*k6[i] + a2606*k7[i] + a2607*k8[i] + a2608*k9[i] + a2609*k10[i] + a2610*k11[i] + a2612*k13[i] + a2613*k14[i] + a2614*k15[i] + a2615*k16[i] + a2616*k17[i] + a2617*k18[i] + a2618*k19[i] + a2619*k20[i] + a2620*k21[i] + a2621*k22[i] + a2622*k23[i] + a2623*k24[i] + a2624*k25[i] + a2625*k26[i])
  end
  f(k27, tmp, p, t + c26*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2700*k1[i] + a2705*k6[i] + a2706*k7[i] + a2707*k8[i] + a2708*k9[i] + a2709*k10[i] + a2711*k12[i] + a2712*k13[i] + a2713*k14[i] + a2714*k15[i] + a2715*k16[i] + a2716*k17[i] + a2717*k18[i] + a2718*k19[i] + a2719*k20[i] + a2720*k21[i] + a2721*k22[i] + a2722*k23[i] + a2723*k24[i] + a2724*k25[i] + a2725*k26[i] + a2726*k27[i])
  end
  f(k28, tmp, p, t + c27*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2800*k1[i] + a2805*k6[i] + a2806*k7[i] + a2807*k8[i] + a2808*k9[i] + a2810*k11[i] + a2811*k12[i] + a2813*k14[i] + a2814*k15[i] + a2815*k16[i] + a2823*k24[i] + a2824*k25[i] + a2825*k26[i] + a2826*k27[i] + a2827*k28[i])
  end
  f(k29, tmp, p, t + c28*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a2900*k1[i] + a2904*k5[i] + a2905*k6[i] + a2906*k7[i] + a2909*k10[i] + a2910*k11[i] + a2911*k12[i] + a2913*k14[i] + a2914*k15[i] + a2915*k16[i] + a2923*k24[i] + a2924*k25[i] + a2925*k26[i] + a2926*k27[i] + a2927*k28[i] + a2928*k29[i])
  end
  f(k30, tmp, p, t + c29*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a3000*k1[i] + a3003*k4[i] + a3004*k5[i] + a3005*k6[i] + a3007*k8[i] + a3009*k10[i] + a3010*k11[i] + a3013*k14[i] + a3014*k15[i] + a3015*k16[i] + a3023*k24[i] + a3024*k25[i] + a3025*k26[i] + a3027*k28[i] + a3028*k29[i] + a3029*k30[i])
  end
  f(k31, tmp, p, t + c30*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a3100*k1[i] + a3102*k3[i] + a3103*k4[i] + a3106*k7[i] + a3107*k8[i] + a3109*k10[i] + a3110*k11[i] + a3113*k14[i] + a3114*k15[i] + a3115*k16[i] + a3123*k24[i] + a3124*k25[i] + a3125*k26[i] + a3127*k28[i] + a3128*k29[i] + a3129*k30[i] + a3130*k31[i])
  end
  f(k32, tmp, p, t + c31*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a3200*k1[i] + a3201*k2[i] + a3204*k5[i] + a3206*k7[i] + a3230*k31[i] + a3231*k32[i])
  end
  f(k33, tmp, p, t + c32*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a3300*k1[i] + a3302*k3[i] + a3332*k33[i])
  end
  f(k34, tmp, p, t + c33*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i] + dt*(a3400*k1[i] + a3401*k2[i] + a3402*k3[i] + a3404*k5[i] + a3406*k7[i] + a3407*k8[i] + a3409*k10[i] + a3410*k11[i] + a3411*k12[i] + a3412*k13[i] + a3413*k14[i] + a3414*k15[i] + a3415*k16[i] + a3416*k17[i] + a3417*k18[i] + a3418*k19[i] + a3419*k20[i] + a3420*k21[i] + a3421*k22[i] + a3422*k23[i] + a3423*k24[i] + a3424*k25[i] + a3425*k26[i] + a3426*k27[i] + a3427*k28[i] + a3428*k29[i] + a3429*k30[i] + a3430*k31[i] + a3431*k32[i] + a3432*k33[i] + a3433*k34[i])
  end
  f(k35, tmp, p, t + c34*dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i] + b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i] + b12*k12[i] + b14*k14[i] + b15*k15[i] + b16*k16[i] + b18*k18[i] + b19*k19[i] + b20*k20[i] + b21*k21[i] + b22*k22[i] + b23*k23[i] + b24*k24[i] + b25*k25[i] + b26*k26[i] + b27*k27[i] + b28*k28[i] + b29*k29[i] + b30*k30[i] + b31*k31[i] + b32*k32[i] + b33*k33[i] + b34*k34[i] + b35*k35[i])
  end
  integrator.destats.nf += 35
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] =  dt*(k2[i] - k34[i]) * adaptiveConst
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  f(integrator.fsallast,u,p,t+dt) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end
