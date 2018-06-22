@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::DP8ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<7 || calcVal
    @unpack c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = cache
    @unpack c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = cache
    @unpack d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = cache
    k1 = f(uprev,p,t)
    k2 = f(uprev+dt*(a0201*k1),p,t+c2*dt)
    k3 = f(uprev+dt*(a0301*k1+a0302*k2),p,t+c3*dt)
    k4 = f(uprev+dt*(a0401*k1       +a0403*k3),p,t+c4*dt)
    k5 = f(uprev+dt*(a0501*k1       +a0503*k3+a0504*k4),p,t+c5*dt)
    k6 = f(uprev+dt*(a0601*k1                +a0604*k4+a0605*k5),p,t+c6*dt)
    k7 = f(uprev+dt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6),p,t+c7*dt)
    k8 = f(uprev+dt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7),p,t+c8*dt)
    k9 = f(uprev+dt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8),p,t+c9*dt)
    k10= f(uprev+dt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9),p,t+c10*dt)
    k11= f(uprev+dt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10),p,t+c11*dt)
    k12= f(uprev+dt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11),p,t+dt)
    kupdate= b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
    utmp = uprev + dt*kupdate
    k13 = f(utmp,p,t+dt)
    k14 = f(uprev+dt*(a1401*k1         +a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13),p,t+c14*dt)
    k15 = f(uprev+dt*(a1501*k1+a1506*k6+a1507*k7+a1508*k8                   +a1511*k11+a1512*k12+a1513*k13+a1514*k14),p,t+c15*dt)
    k16 = f(uprev+dt*(a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9                              +a1613*k13+a1614*k14+a1615*k15),p,t+c16*dt)
    udiff = kupdate
    copyat_or_push!(k,1,udiff)
    bspl = k1 - udiff
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,udiff - k13 - bspl)
    copyat_or_push!(k,4,(d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16))
    copyat_or_push!(k,5,(d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16))
    copyat_or_push!(k,6,(d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16))
    copyat_or_push!(k,7,(d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16))
  end
end

#=
@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::DP8Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<7 || calcVal
    @unpack c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = cache.tab
    @unpack c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = cache.tab
    @unpack d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,kupdate,utilde,tmp = cache
    utmp = utilde
    k = [cache.udiff,cache.bspl,cache.dense_tmp3,cache.dense_tmp4,cache.dense_tmp5,cache.dense_tmp6,cache.dense_tmp7]
    f(k1,uprev,p,t)
    @. tmp = uprev+dt*(a0201*k1)
    f(k2,tmp,p,t+c2*dt)
    @. tmp = uprev+dt*(a0301*k1+a0302*k2)
    f(k3,tmp,p,t+c3*dt)
    @. tmp = uprev+dt*(a0401*k1+a0403*k3)
    f(k4,tmp,p,t+c4*dt)
    @. tmp = uprev+dt*(a0501*k1+a0503*k3+a0504*k4)
    f(k5,tmp,p,t+c5*dt)
    @. tmp = uprev+dt*(a0601*k1+a0604*k4+a0605*k5)
    f(k6,tmp,p,t+c6*dt)
    @. tmp = uprev+dt*(a0701*k1+a0704*k4+a0705*k5+a0706*k6)
    f(k7,tmp,p,t+c7*dt)
    @. tmp = uprev+dt*(a0801*k1+a0804*k4+a0805*k5+a0806*k6+a0807*k7)
    f(k8,tmp,p,t+c8*dt)
    @. tmp = uprev+dt*(a0901*k1+a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8)
    f(k9,tmp,p,t+c9*dt)
    @. tmp = uprev+dt*(a1001*k1+a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9)
    f(k10,tmp,p,t+c10*dt)
    @. tmp = uprev+dt*(a1101*k1+a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
    f(k11,tmp,p,t+c11*dt)
    @. tmp = uprev+dt*(a1201*k1+a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
    f(k12,tmp,p,t+dt)
    @. kupdate = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
    @. utmp = uprev + dt*update
    f(k13,utmp,p,t+dt)
    @. tmp = uprev+dt*(a1401*k1+a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13)
    f(k14,tmp,p,t+c14*dt)
    @. tmp = uprev+dt*(a1501*k1+a1506*k6+a1507*k7+a1508*k8+a1511*k11+a1512*k12+a1513*k13+a1514*k14)
    f(k15,tmp,p,t+c15*dt)
    @. tmp = uprev+dt*(a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9+a1613*k13+a1614*k14+a1615*k15)
    f(k16,tmp,p,t+c16*dt)
    @. udiff= kupdate
    @. bspl = k1 - udiff
    @. k[3] = udiff - k13 - bspl
    @. k[4] = (d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16)
    @. k[5] = (d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16)
    @. k[6] = (d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16)
    @. k[7] = (d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16)
  end
end
=#

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::DP8Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<7 || calcVal
    @unpack c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = cache.tab
    @unpack c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = cache.tab
    @unpack d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,kupdate,utilde,tmp = cache
    utmp = utilde
    k = [cache.udiff,cache.bspl,cache.dense_tmp3,cache.dense_tmp4,cache.dense_tmp5,cache.dense_tmp6,cache.dense_tmp7]
    uidx = eachindex(u)
    f(k1,uprev,p,t)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0201*k1[i])
    end
    f(k2,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0301*k1[i]+a0302*k2[i])
    end
    f(k3,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0401*k1[i]+a0403*k3[i])
    end
    f(k4,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0501*k1[i]+a0503*k3[i]+a0504*k4[i])
    end
    f(k5,tmp,p,t+c5*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0601*k1[i]+a0604*k4[i]+a0605*k5[i])
    end
    f(k6,tmp,p,t+c6*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i])
    end
    f(k7,tmp,p,t+c7*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i])
    end
    f(k8,tmp,p,t+c8*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i])
    end
    f(k9,tmp,p,t+c9*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i])
    end
    f(k10,tmp,p,t+c10*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i])
    end
    f(k11,tmp,p,t+c11*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i])
    end
    f(k12,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds kupdate[i] = b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i]
      @inbounds utmp[i] = uprev[i] + dt*kupdate[i]
    end
    f(k13,utmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1401*k1[i]+a1407*k7[i]+a1408*k8[i]+a1409*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*k13[i])
    end
    f(k14,tmp,p,t+c14*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1501*k1[i]+a1506*k6[i]+a1507*k7[i]+a1508*k8[i]+a1511*k11[i]+a1512*k12[i]+a1513*k13[i]+a1514*k14[i])
    end
    f(k15,tmp,p,t+c15*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a1601*k1[i]+a1606*k6[i]+a1607*k7[i]+a1608*k8[i]+a1609*k9[i]+a1613*k13[i]+a1614*k14[i]+a1615*k15[i])
    end
    f(k16,tmp,p,t+c16*dt)
    @tight_loop_macros for i in uidx
      @inbounds udiff[i]= kupdate[i]
      @inbounds bspl[i] = k1[i] - udiff[i]
      @inbounds k[3][i] = udiff[i] - k13[i] - bspl[i]
      @inbounds k[4][i] = (d401*k1[i]+d406*k6[i]+d407*k7[i]+d408*k8[i]+d409*k9[i]+d410*k10[i]+d411*k11[i]+d412*k12[i]+d413*k13[i]+d414*k14[i]+d415*k15[i]+d416*k16[i])
      @inbounds k[5][i] = (d501*k1[i]+d506*k6[i]+d507*k7[i]+d508*k8[i]+d509*k9[i]+d510*k10[i]+d511*k11[i]+d512*k12[i]+d513*k13[i]+d514*k14[i]+d515*k15[i]+d516*k16[i])
      @inbounds k[6][i] = (d601*k1[i]+d606*k6[i]+d607*k7[i]+d608*k8[i]+d609*k9[i]+d610*k10[i]+d611*k11[i]+d612*k12[i]+d613*k13[i]+d614*k14[i]+d615*k15[i]+d616*k16[i])
      @inbounds k[7][i] = (d701*k1[i]+d706*k6[i]+d707*k7[i]+d708*k8[i]+d709*k9[i]+d710*k10[i]+d711*k11[i]+d712*k12[i]+d713*k13[i]+d714*k14[i]+d715*k15[i]+d716*k16[i])
    end
  end
end
