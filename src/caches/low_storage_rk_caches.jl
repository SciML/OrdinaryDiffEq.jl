
# 2N low storage methods introduced by Williamson
@cache struct LowStorageRK2NCache{uType,rateType,TabType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # tmp acts as second register and fsal both
  tab::TabType
  williamson_condition::Bool
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct LowStorageRK2NConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  A2end::SVector{N,T} # A1 is always zero
  B1::T
  B2end::SVector{N,T}
  c2end::SVector{N,T2} # c1 is always zero
end


function ORK256ConstantCache(T, T2)
  A2 = convert(T, -1.0)
  A3 = convert(T, -1.55798)
  A4 = convert(T, -1.0)
  A5 = convert(T, -0.45031)
  A2end = SVector(A2, A3, A4, A5)

  B1 = convert(T, 0.2)
  B2 = convert(T, 0.83204)
  B3 = convert(T, 0.6)
  B4 = convert(T, 0.35394)
  B5 = convert(T, 0.2)
  B2end = SVector(B2, B3, B4, B5)

  c2 = convert(T2, 0.2)
  c3 = convert(T2, 0.2)
  c4 = convert(T2, 0.8)
  c5 = convert(T2, 0.8)
  c2end = SVector(c2, c3, c4, c5)

  LowStorageRK2NConstantCache{4,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::ORK256,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = ORK256ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::ORK256,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ORK256ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function CarpenterKennedy2N54ConstantCache(T, T2)
  A2 = convert(T, -567301805773//1357537059087)
  A3 = convert(T, -2404267990393//2016746695238)
  A4 = convert(T, -3550918686646//2091501179385)
  A5 = convert(T, -1275806237668//842570457699)
  A2end = SVector(A2, A3, A4, A5)

  B1 = convert(T, 1432997174477//9575080441755)
  B2 = convert(T, 5161836677717//13612068292357)
  B3 = convert(T, 1720146321549//2090206949498)
  B4 = convert(T, 3134564353537//4481467310338)
  B5 = convert(T, 2277821191437//14882151754819)
  B2end = SVector(B2, B3, B4, B5)

  c2 = convert(T2, 1432997174477//9575080441755)
  c3 = convert(T2, 2526269341429//6820363962896)
  c4 = convert(T2, 2006345519317//3224310063776)
  c5 = convert(T2, 2802321613138//2924317926251)
  c2end = SVector(c2, c3, c4, c5)

  LowStorageRK2NConstantCache{4,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = CarpenterKennedy2N54ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CarpenterKennedy2N54ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function SHLDDRK64ConstantCache(T, T2)
  #TODO: Solve the order conditions with more accuracy
  A2 = convert(T, -0.4919575)
  A3 = convert(T, -0.8946264)
  A4 = convert(T, -1.5526678)
  A5 = convert(T, -3.4077973)
  A6 = convert(T, -1.0742640)
  A2end = SVector(A2, A3, A4, A5, A6)

  B1 = convert(T, 0.1453095)
  B2 = convert(T, 0.4653797)
  B3 = convert(T, 0.4675397)
  B4 = convert(T, 0.7795279)
  B5 = convert(T, 0.3574327)
  B6 = convert(T, 0.15)
  B2end = SVector(B2, B3, B4, B5, B6)

  c2 = convert(T2, 0.1453095)
  c3 = convert(T2, 0.3817422)
  c4 = convert(T2, 0.6367813)
  c5 = convert(T2, 0.7560744)
  c6 = convert(T2, 0.9271047)
  c2end = SVector(c2, c3, c4, c5, c6)

  LowStorageRK2NConstantCache{5,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::SHLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = SHLDDRK64ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::SHLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SHLDDRK64ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function DGLDDRK73_CConstantCache(T, T2)
  A2 = convert(T, -0.8083163874983830)
  A3 = convert(T, -1.503407858773331)
  A4 = convert(T, -1.053064525050744)
  A5 = convert(T, -1.463149119280508)
  A6 = convert(T, -0.6592881281087830)
  A7 = convert(T, -1.667891931891068)
  A2end = SVector(A2, A3, A4, A5, A6, A7)

  B1 = convert(T, 0.01197052673097840)
  B2 = convert(T, 0.8886897793820711)
  B3 = convert(T, 0.4578382089261419)
  B4 = convert(T, 0.5790045253338471)
  B5 = convert(T, 0.3160214638138484)
  B6 = convert(T, 0.2483525368264122)
  B7 = convert(T, 0.06771230959408840)
  B2end = SVector(B2, B3, B4, B5, B6, B7)

  c2 = convert(T2, 0.01197052673097840)
  c3 = convert(T2, 0.1823177940361990)
  c4 = convert(T2, 0.5082168062551849)
  c5 = convert(T2, 0.6532031220148590)
  c6 = convert(T2, 0.8534401385678250)
  c7 = convert(T2, 0.9980466084623790)
  c2end = SVector(c2, c3, c4, c5, c6, c7)

  LowStorageRK2NConstantCache{6,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::DGLDDRK73_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = DGLDDRK73_CConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::DGLDDRK73_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  DGLDDRK73_CConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function DGLDDRK84_CConstantCache(T, T2)
  A2  = convert(T, -0.7212962482279240)
  A3  = convert(T, -0.01077336571612980)
  A4  = convert(T, -0.5162584698930970)
  A5  = convert(T, -1.730100286632201)
  A6  = convert(T, -5.200129304403076)
  A7  = convert(T, 0.7837058945416420)
  A8  = convert(T, -0.5445836094332190)
  A2end = SVector(A2, A3, A4, A5, A6, A7, A8)

  B1  = convert(T, 0.2165936736758085)
  B2  = convert(T, 0.1773950826411583)
  B3  = convert(T, 0.01802538611623290)
  B4  = convert(T, 0.08473476372541490)
  B5  = convert(T, 0.8129106974622483)
  B6  = convert(T, 1.903416030422760)
  B7  = convert(T, 0.1314841743399048)
  B8  = convert(T, 0.2082583170674149)
  B2end = SVector(B2, B3, B4, B5, B6, B7, B8)

  c2  = convert(T2, 0.2165936736758085)
  c3  = convert(T2, 0.2660343487538170)
  c4  = convert(T2, 0.2840056122522720)
  c5  = convert(T2, 0.3251266843788570)
  c6  = convert(T2, 0.4555149599187530)
  c7  = convert(T2, 0.7713219317101170)
  c8  = convert(T2, 0.9199028964538660)
  c2end = SVector(c2, c3, c4, c5, c6, c7, c8)

  LowStorageRK2NConstantCache{7,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::DGLDDRK84_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = DGLDDRK84_CConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::DGLDDRK84_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  DGLDDRK84_CConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function DGLDDRK84_FConstantCache(T, T2)
  A2  = convert(T, -0.5534431294501569)
  A3  = convert(T, 0.01065987570203490)
  A4  = convert(T, -0.5515812888932000)
  A5  = convert(T, -1.885790377558741)
  A6  = convert(T, -5.701295742793264)
  A7  = convert(T, 2.113903965664793)
  A8  = convert(T, -0.5339578826675280)
  A2end = SVector(A2, A3, A4, A5, A6, A7, A8)

  B1  = convert(T, 0.08037936882736950)
  B2  = convert(T, 0.5388497458569843)
  B3  = convert(T, 0.01974974409031960)
  B4  = convert(T, 0.09911841297339970)
  B5  = convert(T, 0.7466920411064123)
  B6  = convert(T, 1.679584245618894)
  B7  = convert(T, 0.2433728067008188)
  B8  = convert(T, 0.1422730459001373)
  B2end = SVector(B2, B3, B4, B5, B6, B7, B8)

  c2  = convert(T2, 0.08037936882736950)
  c3  = convert(T2, 0.3210064250338430)
  c4  = convert(T2, 0.3408501826604660)
  c5  = convert(T2, 0.3850364824285470)
  c6  = convert(T2, 0.5040052477534100)
  c7  = convert(T2, 0.6578977561168540)
  c8  = convert(T2, 0.9484087623348481)
  c2end = SVector(c2, c3, c4, c5, c6, c7, c8)

  LowStorageRK2NConstantCache{7,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::DGLDDRK84_F,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = DGLDDRK84_FConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::DGLDDRK84_F,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  DGLDDRK84_FConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function NDBLSRK124ConstantCache(T, T2)
  A2  = convert(T, -0.0923311242368072)
  A3  = convert(T, -0.9441056581158819)
  A4  = convert(T, -4.3271273247576394)
  A5  = convert(T, -2.1557771329026072)
  A6  = convert(T, -0.9770727190189062)
  A7  = convert(T, -0.7581835342571139)
  A8  = convert(T, -1.7977525470825499)
  A9  = convert(T, -2.6915667972700770)
  A10 = convert(T, -4.6466798960268143)
  A11 = convert(T, -0.1539613783825189)
  A12 = convert(T, -0.5943293901830616)
  A2end = SVector(A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12)

  B1  = convert(T, 0.0650008435125904)
  B2  = convert(T, 0.0161459902249842)
  B3  = convert(T, 0.5758627178358159)
  B4  = convert(T, 0.1649758848361671)
  B5  = convert(T, 0.3934619494248182)
  B6  = convert(T, 0.0443509641602719)
  B7  = convert(T, 0.2074504268408778)
  B8  = convert(T, 0.6914247433015102)
  B9  = convert(T, 0.3766646883450449)
  B10 = convert(T, 0.0757190350155483)
  B11 = convert(T, 0.2027862031054088)
  B12 = convert(T, 0.2167029365631842)
  B2end = SVector(B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12)

  c2  = convert(T2, 0.0650008435125904)
  c3  = convert(T2, 0.0796560563081853)
  c4  = convert(T2, 0.1620416710085376)
  c5  = convert(T2, 0.2248877362907778)
  c6  = convert(T2, 0.2952293985641261)
  c7  = convert(T2, 0.3318332506149405)
  c8  = convert(T2, 0.4094724050198658)
  c9  = convert(T2, 0.6356954475753369)
  c10 = convert(T2, 0.6806551557645497)
  c11 = convert(T2, 0.7143773712418350)
  c12 = convert(T2, 0.9032588871651854)
  c2end = SVector(c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12)

  LowStorageRK2NConstantCache{11,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::NDBLSRK124,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = NDBLSRK124ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::NDBLSRK124,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  NDBLSRK124ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function NDBLSRK134ConstantCache(T, T2)
  A2  = convert(T, -0.6160178650170565)
  A3  = convert(T, -0.4449487060774118)
  A4  = convert(T, -1.0952033345276178)
  A5  = convert(T, -1.2256030785959187)
  A6  = convert(T, -0.2740182222332805)
  A7  = convert(T, -0.0411952089052647)
  A8  = convert(T, -0.1797084899153560)
  A9  = convert(T, -1.1771530652064288)
  A10 = convert(T, -0.4078831463120878)
  A11 = convert(T, -0.8295636426191777)
  A12 = convert(T, -4.7895970584252288)
  A13 = convert(T, -0.6606671432964504)
  A2end = SVector(A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13)

  B1  = convert(T, 0.0271990297818803)
  B2  = convert(T, 0.1772488819905108)
  B3  = convert(T, 0.0378528418949694)
  B4  = convert(T, 0.6086431830142991)
  B5  = convert(T, 0.2154313974316100)
  B6  = convert(T, 0.2066152563885843)
  B7  = convert(T, 0.0415864076069797)
  B8  = convert(T, 0.0219891884310925)
  B9  = convert(T, 0.9893081222650993)
  B10 = convert(T, 0.0063199019859826)
  B11 = convert(T, 0.3749640721105318)
  B12 = convert(T, 1.6080235151003195)
  B13 = convert(T, 0.0961209123818189)
  B2end = SVector(B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13)

  c2  = convert(T2, 0.0271990297818803)
  c3  = convert(T2, 0.0952594339119365)
  c4  = convert(T2, 0.1266450286591127)
  c5  = convert(T2, 0.1825883045699772)
  c6  = convert(T2, 0.3737511439063931)
  c7  = convert(T2, 0.5301279418422206)
  c8  = convert(T2, 0.5704177433952291)
  c9  = convert(T2, 0.5885784947099155)
  c10 = convert(T2, 0.6160769826246714)
  c11 = convert(T2, 0.6223252334314046)
  c12 = convert(T2, 0.6897593128753419)
  c13 = convert(T2, 0.9126827615920843)
  c2end = SVector(c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13)

  LowStorageRK2NConstantCache{12,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::NDBLSRK134,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = NDBLSRK134ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::NDBLSRK134,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  NDBLSRK134ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function NDBLSRK144ConstantCache(T, T2)
  A2  = convert(T, -0.7188012108672410)
  A3  = convert(T, -0.7785331173421570)
  A4  = convert(T, -0.0053282796654044)
  A5  = convert(T, -0.8552979934029281)
  A6  = convert(T, -3.9564138245774565)
  A7  = convert(T, -1.5780575380587385)
  A8  = convert(T, -2.0837094552574054)
  A9  = convert(T, -0.7483334182761610)
  A10 = convert(T, -0.7032861106563359)
  A11 = convert(T, 0.0013917096117681)
  A12 = convert(T, -0.0932075369637460)
  A13 = convert(T, -0.9514200470875948)
  A14 = convert(T, -7.1151571693922548)
  A2end = SVector(A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14)

  B1  = convert(T, 0.0367762454319673)
  B2  = convert(T, 0.3136296607553959)
  B3  = convert(T, 0.1531848691869027)
  B4  = convert(T, 0.0030097086818182)
  B5  = convert(T, 0.3326293790646110)
  B6  = convert(T, 0.2440251405350864)
  B7  = convert(T, 0.3718879239592277)
  B8  = convert(T, 0.6204126221582444)
  B9  = convert(T, 0.1524043173028741)
  B10 = convert(T, 0.0760894927419266)
  B11 = convert(T, 0.0077604214040978)
  B12 = convert(T, 0.0024647284755382)
  B13 = convert(T, 0.0780348340049386)
  B14 = convert(T, 5.5059777270269628)
  B2end = SVector(B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14)

  c2  = convert(T2, 0.0367762454319673)
  c3  = convert(T2, 0.1249685262725025)
  c4  = convert(T2, 0.2446177702277698)
  c5  = convert(T2, 0.2476149531070420)
  c6  = convert(T2, 0.2969311120382472)
  c7  = convert(T2, 0.3978149645802642)
  c8  = convert(T2, 0.5270854589440328)
  c9  = convert(T2, 0.6981269994175695)
  c10 = convert(T2, 0.8190890835352128)
  c11 = convert(T2, 0.8527059887098624)
  c12 = convert(T2, 0.8604711817462826)
  c13 = convert(T2, 0.8627060376969976)
  c14 = convert(T2, 0.8734213127600976)
  c2end = SVector(c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14)

  LowStorageRK2NConstantCache{13,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::NDBLSRK144,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = NDBLSRK144ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  tmp = zero(u)
  williamson_condition = alg.williamson_condition
  if calck
    k = zero(rate_prototype)
    williamson_condition = false
  else
    if williamson_condition
      k = tmp
    else
      k = zero(rate_prototype)
    end
  end
  LowStorageRK2NCache(u,uprev,k,tmp,tab,williamson_condition,alg.stage_limiter!,alg.step_limiter!)
end

function alg_cache(alg::NDBLSRK144,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  NDBLSRK144ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end



# 2C low storage methods introduced by Calvo, Franco, Rández (2004)
@cache struct LowStorageRK2CCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct LowStorageRK2CConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  A2end::SVector{N,T} # A1 is always zero
  B1::T
  B2end::SVector{N,T}
  c2end::SVector{N,T2} # c1 is always zero
end



function CFRLDDRK64ConstantCache(T, T2)
  A2 = convert(T, 0.17985400977138)
  A3 = convert(T, 0.14081893152111)
  A4 = convert(T, 0.08255631629428)
  A5 = convert(T, 0.65804425034331)
  A6 = convert(T, 0.31862993413251)
  A2end = SVector(A2, A3, A4, A5, A6)

  B1 = convert(T, 0.10893125722541)
  B2 = convert(T, 0.13201701492152)
  B3 = convert(T, 0.38911623225517)
  B4 = convert(T, -0.59203884581148)
  B5 = convert(T, 0.47385028714844)
  B6 = convert(T, 0.48812405426094)
  B2end = SVector(B2, B3, B4, B5, B6)

  c2 = convert(T2, 0.28878526699679)
  c3 = convert(T2, 0.38176720366804)
  c4 = convert(T2, 0.71262082069639)
  c5 = convert(T2, 0.69606990893393)
  c6 = convert(T2, 0.83050587987157)
  c2end = SVector(c2, c3, c4, c5, c6)

  LowStorageRK2CConstantCache{5,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::CFRLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CFRLDDRK64ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  LowStorageRK2CCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CFRLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CFRLDDRK64ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end


function TSLDDRK74ConstantCache(T, T2)
  A2 = convert(T, 0.241566650129646868)
  A3 = convert(T, 0.0423866513027719953)
  A4 = convert(T, 0.215602732678803776)
  A5 = convert(T, 0.232328007537583987)
  A6 = convert(T, 0.256223412574146438)
  A7 = convert(T, 0.0978694102142697230)
  A2end = SVector(A2, A3, A4, A5, A6, A7)

  B1 = convert(T, 0.0941840925477795334)
  B2 = convert(T, 0.149683694803496998)
  B3 = convert(T, 0.285204742060440058)
  B4 = convert(T, -0.122201846148053668)
  B5 = convert(T, 0.0605151571191401122)
  B6 = convert(T, 0.345986987898399296)
  B7 = convert(T, 0.186627171718797670)
  B2end = SVector(B2, B3, B4, B5, B6, B7)

  c2 = convert(T2, 0.335750742677426401)
  c3 = convert(T2, 0.286254438654048527)
  c4 = convert(T2, 0.744675262090520366)
  c5 = convert(T2, 0.639198690801246909)
  c6 = convert(T2, 0.723609252956949472)
  c7 = convert(T2, 0.91124223849547205)
  c2end = SVector(c2, c3, c4, c5, c6, c7)

  LowStorageRK2CConstantCache{6,T,T2}(A2end, B1, B2end, c2end)
end

function alg_cache(alg::TSLDDRK74,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = TSLDDRK74ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  LowStorageRK2CCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::TSLDDRK74,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  TSLDDRK74ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end



# 3S low storage methods introduced by Ketcheson
@cache struct LowStorageRK3SCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct LowStorageRK3SConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  γ12end::SVector{N,T} # γ11 is always zero
  γ22end::SVector{N,T} # γ21 is always one
  γ32end::SVector{N,T} # γ31 is always zero
  # TODO: γ302 == γ303 == 0 in all emthods implemented below -> possible optimisation?
  δ2end ::SVector{N,T} # δ1  is always one
  β1::T
  β2end::SVector{N,T}
  c2end::SVector{N,T2} # c1 is always zero
end


function ParsaniKetchesonDeconinck3S32ConstantCache(T, T2)
  γ102 = convert(T, -1.2664395576322218e-1)
  γ103 = convert(T, 1.1426980685848858e+0)
  γ12end = SVector(γ102, γ103)

  γ202 = convert(T, 6.5427782599406470e-1)
  γ203 = convert(T, -8.2869287683723744e-2)
  γ22end = SVector(γ202, γ203)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ32end = SVector(γ302, γ303)

  δ02 = convert(T, 7.2196567116037724e-1)
  δ03 = convert(T, 0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03)

  β1  = convert(T, 7.2366074728360086e-1)
  β02 = convert(T, 3.4217876502651023e-1)
  β03 = convert(T, 3.6640216242653251e-1)
  β2end = SVector(β02, β03)

  c02 = convert(T2, 7.2366074728360086e-1)
  c03 = convert(T2, 5.9236433182015646e-1)
  c2end = SVector(c02, c03)

  LowStorageRK3SConstantCache{2,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S32ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S32ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S82ConstantCache(T, T2)
  γ102 = convert(T, 4.2397552118208004e-1)
  γ103 = convert(T, -2.3528852074619033e-1)
  γ104 = convert(T, 7.9598685017877846e-1)
  γ105 = convert(T, -1.3205224623823271e+0)
  γ106 = convert(T, 2.1452956294251941e+0)
  γ107 = convert(T, -9.5532770501880648e-1)
  γ108 = convert(T, 2.5361391125131094e-1)
  γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108)

  γ202 = convert(T, 4.4390665802303775e-1)
  γ203 = convert(T, 7.5333732286056154e-1)
  γ204 = convert(T, 6.5885460813015481e-2)
  γ205 = convert(T, 6.3976199384289623e-1)
  γ206 = convert(T, -7.3823030755143193e-1)
  γ207 = convert(T, 7.0177211879534529e-1)
  γ208 = convert(T, 4.0185379950224559e-1)
  γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ304 = convert(T, 5.8415358412023582e-2)
  γ305 = convert(T, 6.4219008773865116e-1)
  γ306 = convert(T, 6.8770305706885126e-1)
  γ307 = convert(T, 6.3729822311671305e-2)
  γ308 = convert(T, -3.3679429978131387e-1)
  γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308)

  δ02 = convert(T, 2.9762522910396538e-1)
  δ03 = convert(T, 3.4212961014330662e-1)
  δ04 = convert(T, 5.7010739154759105e-1)
  δ05 = convert(T, 4.1350769551529132e-1)
  δ06 = convert(T, -1.4040672669058066e-1)
  δ07 = convert(T, 2.1249567092409008e-1)
  δ08 = convert(T, 0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08)

  β1  = convert(T, 9.9292229393265474e-1)
  β02 = convert(T, 5.2108385130005974e-1)
  β03 = convert(T, 3.8505327083543915e-3)
  β04 = convert(T, 7.9714199213087467e-1)
  β05 = convert(T, -8.1822460276649120e-2)
  β06 = convert(T, 8.4604310411858186e-1)
  β07 = convert(T, -1.0191166090841246e-1)
  β08 = convert(T, 6.3190236038107500e-2)
  β2end = SVector(β02, β03, β04, β05, β06, β07, β08)

  c02 = convert(T2, 9.9292229393265474e-1)
  c03 = convert(T2, 1.0732413280565014e+0)
  c04 = convert(T2, 2.5057060509809409e-1)
  c05 = convert(T2, 1.0496674928979783e+0)
  c06 = convert(T2, -6.7488037049720317e-1)
  c07 = convert(T2, -1.5868411612120166e+0)
  c08 = convert(T2, 2.1138242369563969e+0)
  c2end = SVector(c02, c03, c04, c05, c06, c07, c08)

  LowStorageRK3SConstantCache{7,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S82,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S82ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S82,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S82ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S53ConstantCache(T, T2)
  γ102 = convert(T, 2.5876919610938998e-1)
  γ103 = convert(T, -1.3243708384977859e-1)
  γ104 = convert(T, 5.0556648948362981e-2)
  γ105 = convert(T, 5.6705507883024708e-1)
  γ12end = SVector(γ102, γ103, γ104, γ105)

  γ202 = convert(T, 5.5284013909611196e-1)
  γ203 = convert(T, 6.7318513326032769e-1)
  γ204 = convert(T, 2.8031054965521607e-1)
  γ205 = convert(T, 5.5215115815918758e-1)
  γ22end = SVector(γ202, γ203, γ204, γ205)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ304 = convert(T, 2.7525797946334213e-1)
  γ305 = convert(T, -8.9505445022148511e-1)
  γ32end = SVector(γ302, γ303, γ304, γ305)

  δ02 = convert(T, 3.4076878915216791e-1)
  δ03 = convert(T, 3.4143871647890728e-1)
  δ04 = convert(T, 7.2292984084963252e-1)
  δ05 = convert(T, 0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05)

  β1  = convert(T, 2.3002859824852059e-1)
  β02 = convert(T, 3.0214498165167158e-1)
  β03 = convert(T, 8.0256010238856679e-1)
  β04 = convert(T, 4.3621618871511753e-1)
  β05 = convert(T, 1.1292705979513513e-1)
  β2end = SVector(β02, β03, β04, β05)

  c02 = convert(T2, 2.3002859824852059e-1)
  c03 = convert(T2, 4.0500453764839639e-1)
  c04 = convert(T2, 8.9478204142351003e-1)
  c05 = convert(T2, 7.2351146275625733e-1)
  c2end = SVector(c02, c03, c04, c05)

  LowStorageRK3SConstantCache{4,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S53ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S53ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S173ConstantCache(T, T2)
  γ102 = convert(T, 7.9377023961829174e-1)
  γ103 = convert(T, -8.3475116244241754e-2)
  γ104 = convert(T, -1.6706337980062214e-2)
  γ105 = convert(T, 3.6410691500331427e-1)
  γ106 = convert(T, 6.9178255181542780e-1)
  γ107 = convert(T, 1.4887115004739182e+0)
  γ108 = convert(T, 4.5336125560871188e-1)
  γ109 = convert(T, -1.2705776046458739e-1)
  γ110 = convert(T, 8.3749845457747696e-1)
  γ111 = convert(T, 1.5709218393361746e-1)
  γ112 = convert(T, -5.7768207086288348e-1)
  γ113 = convert(T, -5.7340394122375393e-1)
  γ114 = convert(T, -1.2050734846514470e+0)
  γ115 = convert(T, -2.8100719513641002e+0)
  γ116 = convert(T, 1.6142798657609492e-1)
  γ117 = convert(T, -2.5801264756641613e+0)
  γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108, γ109, γ110, γ111, γ112, γ113, γ114, γ115, γ116, γ117)

  γ202 = convert(T, 3.2857861940811250e-1)
  γ203 = convert(T, 1.1276843361180819e+0)
  γ204 = convert(T, 1.3149447395238016e+0)
  γ205 = convert(T, 5.2062891534209055e-1)
  γ206 = convert(T, 8.8127462325164985e-1)
  γ207 = convert(T, 4.2020606445856712e-1)
  γ208 = convert(T, 7.6532635739246124e-2)
  γ209 = convert(T, 4.4386734924685722e-1)
  γ210 = convert(T, 6.6503093955199682e-2)
  γ211 = convert(T, 1.5850209163184039e+0)
  γ212 = convert(T, 1.1521721573462576e+0)
  γ213 = convert(T, 1.1172750819374575e+0)
  γ214 = convert(T, 7.7630223917584007e-1)
  γ215 = convert(T, 1.0046657060652295e+0)
  γ216 = convert(T, -1.9795868964959054e-1)
  γ217 = convert(T, 1.3350583594705518e+0)
  γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208, γ209, γ210, γ211, γ212, γ213, γ214, γ215, γ216, γ217)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ304 = convert(T, 8.4034574578399479e-1)
  γ305 = convert(T, 8.5047738439705145e-1)
  γ306 = convert(T, 1.4082448501410852e-1)
  γ307 = convert(T, -3.2678802469519369e-1)
  γ308 = convert(T, 5.3716357620635535e-1)
  γ309 = convert(T, 9.0228922115199051e-1)
  γ310 = convert(T, 1.5960226946983552e-1)
  γ311 = convert(T, 1.1038153140686748e+0)
  γ312 = convert(T, 1.0843516423068365e-1)
  γ313 = convert(T, 4.6212710442787724e-1)
  γ314 = convert(T, -3.3448312125108398e-1)
  γ315 = convert(T, 1.1153826567096696e+0)
  γ316 = convert(T, 1.5503248734613539e+0)
  γ317 = convert(T, -1.2200245424704212e+0)
  γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308, γ309, γ310, γ311, γ312, γ313, γ314, γ315, γ316, γ317)

  δ02 = convert(T, -3.7235794357769936e-1)
  δ03 = convert(T, 3.3315440189685536e-1)
  δ04 = convert(T, -8.2667630338402520e-1)
  δ05 = convert(T, -5.4628377681035534e-1)
  δ06 = convert(T, 6.0210777634642887e-1)
  δ07 = convert(T, -5.7528717894031067e-1)
  δ08 = convert(T, 5.0914861529202782e-1)
  δ09 = convert(T, 3.8258114767897194e-1)
  δ10 = convert(T, -4.6279063221185290e-1)
  δ11 = convert(T, -2.0820434288562648e-1)
  δ12 = convert(T, 1.4398056081552713e+0)
  δ13 = convert(T, -2.8056600927348752e-1)
  δ14 = convert(T, 2.2767189929551406e+0)
  δ15 = convert(T, -5.8917530100546356e-1)
  δ16 = convert(T, 9.1328651048418164e-1)
  δ17 = convert(T, 0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08, δ09, δ10, δ11, δ12, δ13, δ14, δ15, δ16, δ17)

  β1  = convert(T, 4.9565403010221741e-2)
  β02 = convert(T, 9.7408718698159397e-2)
  β03 = convert(T, -1.7620737976801870e-1)
  β04 = convert(T, 1.4852069175460250e-1)
  β05 = convert(T, -3.3127657103714951e-2)
  β06 = convert(T, 4.8294609330498492e-2)
  β07 = convert(T, 4.9622612199980112e-2)
  β08 = convert(T, 8.7340766269850378e-1)
  β09 = convert(T, -2.8692804399085370e-1)
  β10 = convert(T, 1.2679897532256112e+0)
  β11 = convert(T, -1.0217436118953449e-2)
  β12 = convert(T, 8.4665570032598350e-2)
  β13 = convert(T, 2.8253854742588246e-2)
  β14 = convert(T, -9.2936733010804407e-2)
  β15 = convert(T, -8.4798124766803512e-2)
  β16 = convert(T, -1.6923145636158564e-2)
  β17 = convert(T, -4.7305106233879957e-2)
  β2end = SVector(β02, β03, β04, β05, β06, β07, β08, β09, β10, β11, β12, β13, β14, β15, β16, β17)

  c02 = convert(T2, 4.9565403010221741e-2)
  c03 = convert(T2, 1.3068799001687578e-1)
  c04 = convert(T2, -1.5883063460310493e-1)
  c05 = convert(T2, 3.5681144740196935e-1)
  c06 = convert(T2, 7.6727123317642698e-2)
  c07 = convert(T2, 1.0812579255374613e-1)
  c08 = convert(T2, 1.8767228084815801e-1)
  c09 = convert(T2, 9.6162976936182631e-1)
  c10 = convert(T2, -2.2760719867560897e-1)
  c11 = convert(T2, 1.1115681606027146e+0)
  c12 = convert(T2, 6.1266845427676520e-1)
  c13 = convert(T2, 1.0729473245077408e+0)
  c14 = convert(T2, 3.7824186468104548e-1)
  c15 = convert(T2, 7.9041891347646720e-1)
  c16 = convert(T2, -1.0406955693161675e+0)
  c17 = convert(T2, -2.4607146824557105e-1)
  c2end = SVector(c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16, c17)

  LowStorageRK3SConstantCache{16,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S173,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S173ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S173,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S173ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S94ConstantCache(T, T2)
  γ102 = convert(T, -4.6556413837561301e+0)
  γ103 = convert(T, -7.7202649689034453e-1)
  γ104 = convert(T, -4.0244202720632174e+0)
  γ105 = convert(T, -2.1296873883702272e-2)
  γ106 = convert(T, -2.4350219407769953e+0)
  γ107 = convert(T,  1.9856336960249132e-2)
  γ108 = convert(T, -2.8107894116913812e-1)
  γ109 = convert(T,  1.6894354373677900e-1)
  γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108, γ109)

  γ202 = convert(T, 2.4992627683300688e+0)
  γ203 = convert(T, 5.8668202764174726e-1)
  γ204 = convert(T, 1.2051419816240785e+0)
  γ205 = convert(T, 3.4747937498564541e-1)
  γ206 = convert(T, 1.3213458736302766e+0)
  γ207 = convert(T, 3.1196363453264964e-1)
  γ208 = convert(T, 4.3514189245414447e-1)
  γ209 = convert(T, 2.3596980658341213e-1)
  γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208, γ209)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ304 = convert(T,  7.6209857891449362e-1)
  γ305 = convert(T, -1.9811817832965520e-1)
  γ306 = convert(T, -6.2289587091629484e-1)
  γ307 = convert(T, -3.7522475499063573e-1)
  γ308 = convert(T, -3.3554373281046146e-1)
  γ309 = convert(T, -4.5609629702116454e-2)
  γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308, γ309)

  δ02 = convert(T,  1.2629238731608268e+0)
  δ03 = convert(T,  7.5749675232391733e-1)
  δ04 = convert(T,  5.1635907196195419e-1)
  δ05 = convert(T, -2.7463346616574083e-2)
  δ06 = convert(T, -4.3826743572318672e-1)
  δ07 = convert(T,  1.2735870231839268e+0)
  δ08 = convert(T, -6.2947382217730230e-1)
  δ09 = convert(T,  0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08, δ09)

  β1  = convert(T,  2.8363432481011769e-1)
  β02 = convert(T,  9.7364980747486463e-1)
  β03 = convert(T,  3.3823592364196498e-1)
  β04 = convert(T, -3.5849518935750763e-1)
  β05 = convert(T, -4.1139587569859462e-3)
  β06 = convert(T,  1.4279689871485013e+0)
  β07 = convert(T,  1.8084680519536503e-2)
  β08 = convert(T,  1.6057708856060501e-1)
  β09 = convert(T,  2.9522267863254809e-1)
  β2end = SVector(β02, β03, β04, β05, β06, β07, β08, β09)

  c02 = convert(T2,  2.8363432481011769e-1)
  c03 = convert(T2,  5.4840742446661772e-1)
  c04 = convert(T2,  3.6872298094969475e-1)
  c05 = convert(T2, -6.8061183026103156e-1)
  c06 = convert(T2,  3.5185265855105619e-1)
  c07 = convert(T2,  1.6659419385562171e+0)
  c08 = convert(T2,  9.7152778807463247e-1)
  c09 = convert(T2,  9.0515694340066954e-1)
  c2end = SVector(c02, c03, c04, c05, c06, c07, c08, c09)

  LowStorageRK3SConstantCache{8,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S94,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S94ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S94,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S94ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S184ConstantCache(T, T2)
  γ102 = convert(T, 1.1750819811951678e+0)
  γ103 = convert(T, 3.0909017892654811e-1)
  γ104 = convert(T, 1.4409117788115862e+0)
  γ105 = convert(T, -4.3563049445694069e-1)
  γ106 = convert(T, 2.0341503014683893e-1)
  γ107 = convert(T, 4.9828356971917692e-1)
  γ108 = convert(T, 3.5307737157745489e+0)
  γ109 = convert(T, -7.9318790975894626e-1)
  γ110 = convert(T, 8.9120513355345166e-1)
  γ111 = convert(T, 5.7091009196320974e-1)
  γ112 = convert(T, 1.6912188575015419e-2)
  γ113 = convert(T, 1.0077912519329719e+0)
  γ114 = convert(T, -6.8532953752099512e-1)
  γ115 = convert(T, 1.0488165551884063e+0)
  γ116 = convert(T, 8.3647761371829943e-1)
  γ117 = convert(T, 1.3087909830445710e+0)
  γ118 = convert(T, 9.0419681700177323e-1)
  γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108, γ109, γ110, γ111, γ112, γ113, γ114, γ115, γ116, γ117, γ118)

  γ202 = convert(T, -1.2891068509748144e-1)
  γ203 = convert(T, 3.5609406666728954e-1)
  γ204 = convert(T, -4.0648075226104241e-1)
  γ205 = convert(T, 6.0714786995207426e-1)
  γ206 = convert(T, 1.0253501186236846e+0)
  γ207 = convert(T, 2.4411240760769423e-1)
  γ208 = convert(T, -1.2813606970134104e+0)
  γ209 = convert(T, 8.1625711892373898e-1)
  γ210 = convert(T, 1.0171269354643386e-1)
  γ211 = convert(T, 1.9379378662711269e-1)
  γ212 = convert(T, 7.4408643544851782e-1)
  γ213 = convert(T, -1.2591764563430008e-1)
  γ214 = convert(T, 1.1996463179654226e+0)
  γ215 = convert(T, 4.5772068865370406e-2)
  γ216 = convert(T, 8.3622292077033844e-1)
  γ217 = convert(T, -1.4179124272450148e+0)
  γ218 = convert(T, 1.3661459065331649e-1)
  γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208, γ209, γ210, γ211, γ212, γ213, γ214, γ215, γ216, γ217, γ218)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ304 = convert(T, 2.5583378537249163e-1)
  γ305 = convert(T, 5.2676794366988289e-1)
  γ306 = convert(T, -2.5648375621792202e-1)
  γ307 = convert(T, 3.1932438003236391e-1)
  γ308 = convert(T, -3.1106815010852862e-1)
  γ309 = convert(T, 4.7631196164025996e-1)
  γ310 = convert(T, -9.8853727938895783e-2)
  γ311 = convert(T, 1.9274726276883622e-1)
  γ312 = convert(T, 3.2389860855971508e-2)
  γ313 = convert(T, 7.5923980038397509e-2)
  γ314 = convert(T, 2.0635456088664017e-1)
  γ315 = convert(T, -8.9741032556032857e-2)
  γ316 = convert(T, 2.6899932505676190e-2)
  γ317 = convert(T, 4.1882069379552307e-2)
  γ318 = convert(T, 6.2016148912381761e-2)
  γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308, γ309, γ310, γ311, γ312, γ313, γ314, γ315, γ316, γ317, γ318)

  δ02 = convert(T, 3.5816500441970289e-1)
  δ03 = convert(T, 5.8208024465093577e-1)
  δ04 = convert(T, -2.2615285894283538e-1)
  δ05 = convert(T, -2.1715466578266213e-1)
  δ06 = convert(T, -4.6990441450888265e-1)
  δ07 = convert(T, -2.7986911594744995e-1)
  δ08 = convert(T, 9.8513926355272197e-1)
  δ09 = convert(T, -1.1899324232814899e-1)
  δ10 = convert(T, 4.2821073124370562e-1)
  δ11 = convert(T, -8.2196355299900403e-1)
  δ12 = convert(T, 5.8113997057675074e-2)
  δ13 = convert(T, -6.1283024325436919e-1)
  δ14 = convert(T, 5.6800136190634054e-1)
  δ15 = convert(T, -3.3874970570335106e-1)
  δ16 = convert(T, -7.3071238125137772e-1)
  δ17 = convert(T, 8.3936016960374532e-2)
  δ18 = convert(T, 0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08, δ09, δ10, δ11, δ12, δ13, δ14, δ15, δ16, δ17, δ18)

  β1  = convert(T, 1.2384169480626298e-1)
  β02 = convert(T, 1.0176262534280349e+0)
  β03 = convert(T, -6.9732026387527429e-2)
  β04 = convert(T, 3.4239356067806476e-1)
  β05 = convert(T, 1.8177707207807942e-2)
  β06 = convert(T, -6.1188746289480445e-3)
  β07 = convert(T, 7.8242308902580354e-2)
  β08 = convert(T, -3.7642864750532951e-1)
  β09 = convert(T, -4.5078383666690258e-2)
  β10 = convert(T, -7.5734228201432585e-1)
  β11 = convert(T, -2.7149222760935121e-1)
  β12 = convert(T, 1.1833684341657344e-3)
  β13 = convert(T, 2.8858319979308041e-2)
  β14 = convert(T, 4.6005267586974657e-1)
  β15 = convert(T, 1.8014887068775631e-2)
  β16 = convert(T, -1.5508175395461857e-2)
  β17 = convert(T, -4.0095737929274988e-1)
  β18 = convert(T, 1.4949678367038011e-1)
  β2end = SVector(β02, β03, β04, β05, β06, β07, β08, β09, β10, β11, β12, β13, β14, β15, β16, β17, β18)

  c02 = convert(T2, 1.2384169480626298e-1)
  c03 = convert(T2, 1.1574324659554065e+0)
  c04 = convert(T2, 5.4372099141546926e-1)
  c05 = convert(T2, 8.8394666834280744e-1)
  c06 = convert(T2, -1.2212042176605774e-1)
  c07 = convert(T2, 4.4125685133082082e-1)
  c08 = convert(T2, 3.8039092095473748e-1)
  c09 = convert(T2, 5.4591107347528367e-2)
  c10 = convert(T2, 4.8731855535356028e-1)
  c11 = convert(T2, -2.3007964303896034e-1)
  c12 = convert(T2, -1.8907656662915873e-1)
  c13 = convert(T2, 8.1059805668623763e-1)
  c14 = convert(T2, 7.7080875997868803e-1)
  c15 = convert(T2, 1.1712158507200179e+0)
  c16 = convert(T2, 1.2755351018003545e+0)
  c17 = convert(T2, 8.0422507946168564e-1)
  c18 = convert(T2, 9.7508680250761848e-1)
  c2end = SVector(c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16, c17, c18)

  LowStorageRK3SConstantCache{17,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S184,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S184ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S184,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S184ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S105ConstantCache(T, T2)
  γ102 = convert(T,  4.0436600785287713e-1)
  γ103 = convert(T, -8.5034274641295027e-1)
  γ104 = convert(T, -6.9508941671218478e+0)
  γ105 = convert(T,  9.2387652252320684e-1)
  γ106 = convert(T, -2.5631780399589106e+0)
  γ107 = convert(T,  2.5457448699988827e-1)
  γ108 = convert(T,  3.1258317336761454e-1)
  γ109 = convert(T, -7.0071148003175443e-1)
  γ110 = convert(T,  4.8396209710057070e-1)
  γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108, γ109, γ110)

  γ202 = convert(T, 6.8714670697294733e-1)
  γ203 = convert(T, 1.0930247604585732e+0)
  γ204 = convert(T, 3.2259753823377983e+0)
  γ205 = convert(T, 1.0411537008416110e+0)
  γ206 = convert(T, 1.2928214888638039e+0)
  γ207 = convert(T, 7.3914627692888835e-1)
  γ208 = convert(T, 1.2391292570651462e-1)
  γ209 = convert(T, 1.8427534793568445e-1)
  γ210 = convert(T, 5.7127889427161162e-2)
  γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208, γ209, γ210)

  γ302 = convert(T,  0.0000000000000000e+0)
  γ303 = convert(T,  0.0000000000000000e+0)
  γ304 = convert(T, -2.3934051593398129e+0)
  γ305 = convert(T, -1.9028544220991284e+0)
  γ306 = convert(T, -2.8200422105835639e+0)
  γ307 = convert(T, -1.8326984641282289e+0)
  γ308 = convert(T, -2.1990945108072310e-1)
  γ309 = convert(T, -4.0824306603783045e-1)
  γ310 = convert(T, -1.3776697911236280e-1)
  γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308, γ309, γ310)

  δ02 = convert(T, -1.3317784091400336e-1)
  δ03 = convert(T,  8.2604227852898304e-1)
  δ04 = convert(T,  1.5137004305165804e+0)
  δ05 = convert(T, -1.3058100631721905e+0)
  δ06 = convert(T,  3.0366787893355149e+0)
  δ07 = convert(T, -1.4494582670831953e+0)
  δ08 = convert(T,  3.8343138733685103e+0)
  δ09 = convert(T,  4.1222939718018692e+0)
  δ10 = convert(T,  0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08, δ09, δ10)

  β1  = convert(T, 2.5978835757039448e-1)
  β02 = convert(T, 1.7770088002098183e-2)
  β03 = convert(T, 2.4816366373161344e-1)
  β04 = convert(T, 7.9417368275785671e-1)
  β05 = convert(T, 3.8853912968701337e-1)
  β06 = convert(T, 1.4550516642704694e-1)
  β07 = convert(T, 1.5875173794655811e-1)
  β08 = convert(T, 1.6506056315937651e-1)
  β09 = convert(T, 2.1180932999328042e-1)
  β10 = convert(T, 1.5593923403495016e-1)
  β2end = SVector(β02, β03, β04, β05, β06, β07, β08, β09, β10)

  c02 = convert(T2, 2.5978835757039448e-1)
  c03 = convert(T2, 9.9045731158085557e-2)
  c04 = convert(T2, 2.1555118823045644e-1)
  c05 = convert(T2, 5.0079500784155040e-1)
  c06 = convert(T2, 5.5922519148547800e-1)
  c07 = convert(T2, 5.4499869734044426e-1)
  c08 = convert(T2, 7.6152246625852738e-1)
  c09 = convert(T2, 8.4270620830633836e-1)
  c10 = convert(T2, 9.1522098071770008e-1)
  c2end = SVector(c02, c03, c04, c05, c06, c07, c08, c09, c10)

  LowStorageRK3SConstantCache{9,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S105,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S105ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S105,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S105ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S205ConstantCache(T, T2)
  γ102 = convert(T, -1.1682479703229380e+0)
  γ103 = convert(T, -2.5112155037089772e+0)
  γ104 = convert(T, -5.5259960154735988e-1)
  γ105 = convert(T, 2.9243033509511740e-3)
  γ106 = convert(T, -4.7948973385386493e+0)
  γ107 = convert(T, -5.3095533497183016e+0)
  γ108 = convert(T, -2.3624194456630736e+0)
  γ109 = convert(T, 2.0068995756589547e-1)
  γ110 = convert(T, -1.4985808661597710e+0)
  γ111 = convert(T, 4.8941228502377687e-1)
  γ112 = convert(T, -1.0387512755259576e-1)
  γ113 = convert(T, -1.3287664273288191e-1)
  γ114 = convert(T, 7.5858678822837511e-1)
  γ115 = convert(T, -4.3321586294096939e+0)
  γ116 = convert(T, 4.8199700138402146e-1)
  γ117 = convert(T, -7.0924756614960671e-3)
  γ118 = convert(T, -8.8422252029506054e-1)
  γ119 = convert(T, -8.9129367099545231e-1)
  γ120 = convert(T, 1.5297157134040762e+0)
  γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108, γ109, γ110, γ111, γ112, γ113, γ114, γ115, γ116, γ117, γ118, γ119, γ120)

  γ202 = convert(T, 8.8952052154583572e-1)
  γ203 = convert(T, 8.8988129100385194e-1)
  γ204 = convert(T, 3.5701564494677057e-1)
  γ205 = convert(T, 2.4232462479216824e-1)
  γ206 = convert(T, 1.2727083024258155e+0)
  γ207 = convert(T, 1.1126977210342681e+0)
  γ208 = convert(T, 5.1360709645409097e-1)
  γ209 = convert(T, 1.1181089682044856e-1)
  γ210 = convert(T, 2.7881272382085232e-1)
  γ211 = convert(T, 4.9032886260666715e-2)
  γ212 = convert(T, 4.1871051065897870e-2)
  γ213 = convert(T, 4.4602463796686219e-2)
  γ214 = convert(T, 1.4897271251154750e-2)
  γ215 = convert(T, 2.6244269699436817e-1)
  γ216 = convert(T, -4.7486056986590294e-3)
  γ217 = convert(T, 2.3219312682036197e-2)
  γ218 = convert(T, 6.2852588972458059e-2)
  γ219 = convert(T, 5.4473719351268962e-2)
  γ220 = convert(T, 2.4345446089014514e-2)
  γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208, γ209, γ210, γ211, γ212, γ213, γ214, γ215, γ216, γ217, γ218, γ219, γ220)

  γ302 = convert(T, 0.0000000000000000e+0)
  γ303 = convert(T, 0.0000000000000000e+0)
  γ304 = convert(T, 1.9595487007932735e-1)
  γ305 = convert(T, -6.9871675039100595e-5)
  γ306 = convert(T, 1.0592231169810050e-1)
  γ307 = convert(T, 1.0730426871909635e+0)
  γ308 = convert(T, 8.9257826744389124e-1)
  γ309 = convert(T, -1.4078912484894415e-1)
  γ310 = convert(T, -2.6869890558434262e-1)
  γ311 = convert(T, -6.5175753568318007e-2)
  γ312 = convert(T, 4.9177812903108553e-1)
  γ313 = convert(T, 4.6017684776493678e-1)
  γ314 = convert(T, -6.4689512947008251e-3)
  γ315 = convert(T, 4.4034728024115377e-1)
  γ316 = convert(T, 6.1086885767527943e-1)
  γ317 = convert(T, 5.0546454457410162e-1)
  γ318 = convert(T, 5.4668509293072887e-1)
  γ319 = convert(T, 7.1414182420995431e-1)
  γ320 = convert(T, -1.0558095282893749e+0)
  γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308, γ309, γ310, γ311, γ312, γ313, γ314, γ315, γ316, γ317, γ318, γ319, γ320)

  δ02 = convert(T, 1.4375468781258596e+0)
  δ03 = convert(T, 1.5081653637261594e+0)
  δ04 = convert(T, -1.4575347066062688e-1)
  δ05 = convert(T, 3.1495761082838158e-1)
  δ06 = convert(T, 3.5505919368536931e-1)
  δ07 = convert(T, 2.3616389374566960e-1)
  δ08 = convert(T, 1.0267488547302055e-1)
  δ09 = convert(T, 3.5991243524519438e+0)
  δ10 = convert(T, 1.5172890003890782e+0)
  δ11 = convert(T, 1.8171662741779953e+0)
  δ12 = convert(T, 2.8762263521436831e+0)
  δ13 = convert(T, 4.6350154228218754e-1)
  δ14 = convert(T, 1.5573122110727220e+0)
  δ15 = convert(T, 2.0001066778080254e+0)
  δ16 = convert(T, 9.1690694855534305e-1)
  δ17 = convert(T, 2.0474618401365854e+0)
  δ18 = convert(T, -3.2336329115436924e-1)
  δ19 = convert(T, 3.2899060754742177e-1)
  δ20 = convert(T, 0.0000000000000000e+0)
  δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08, δ09, δ10, δ11, δ12, δ13, δ14, δ15, δ16, δ17, δ18, δ19, δ20)

  β1  = convert(T, 1.7342385375780556e-1)
  β02 = convert(T, 2.8569004728564801e-1)
  β03 = convert(T, 6.8727044379779589e-1)
  β04 = convert(T, 1.2812121060977319e-1)
  β05 = convert(T, 4.9137180740403122e-4)
  β06 = convert(T, 4.7033584446956857e-2)
  β07 = convert(T, 4.4539998128170821e-1)
  β08 = convert(T, 1.2259824887343720e+0)
  β09 = convert(T, 2.0616463985024421e-2)
  β10 = convert(T, 1.5941162575324802e-1)
  β11 = convert(T, 1.2953803678226099e+0)
  β12 = convert(T, 1.7287352967302603e-3)
  β13 = convert(T, 1.1660483420536467e-1)
  β14 = convert(T, 7.7997036621815521e-2)
  β15 = convert(T, 3.2563250234418012e-1)
  β16 = convert(T, 1.0611520488333197e+0)
  β17 = convert(T, 6.5891625628040993e-4)
  β18 = convert(T, 8.3534647700054046e-2)
  β19 = convert(T, 9.8972579458252483e-2)
  β20 = convert(T, 4.3010116145097040e-2)
  β2end = SVector(β02, β03, β04, β05, β06, β07, β08, β09, β10, β11, β12, β13, β14, β15, β16, β17, β18, β19, β20)

  c02 = convert(T2, 1.7342385375780556e-1)
  c03 = convert(T2, 3.0484982420032158e-1)
  c04 = convert(T2, 5.5271395645729193e-1)
  c05 = convert(T2, 4.7079204549750037e-2)
  c06 = convert(T2, 1.5652540451324129e-1)
  c07 = convert(T2, 1.8602224049074517e-1)
  c08 = convert(T2, 2.8426620035751449e-1)
  c09 = convert(T2, 9.5094727548792268e-1)
  c10 = convert(T2, 6.8046501070096010e-1)
  c11 = convert(T2, 5.9705366562360063e-1)
  c12 = convert(T2, 1.8970821645077285e+0)
  c13 = convert(T2, 2.9742664004529606e-1)
  c14 = convert(T2, 6.0813463700134940e-1)
  c15 = convert(T2, 7.3080004188477765e-1)
  c16 = convert(T2, 9.1656999044951792e-1)
  c17 = convert(T2, 1.4309687554614530e+0)
  c18 = convert(T2, 4.1043824968249148e-1)
  c19 = convert(T2, 8.4898255952298962e-1)
  c20 = convert(T2, 3.3543896258348421e-1)
  c2end = SVector(c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20)

  LowStorageRK3SConstantCache{19,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S205,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = zero(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S205ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S205,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  ParsaniKetchesonDeconinck3S205ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end



# 3S+ low storage methods: 3S methods adding another memory location for the embedded method (non-FSAL version)
# ## References
# - Ranocha, Dalcin, Parsani, Ketcheson (2021)
#   Optimized Runge-Kutta Methods with Automatic Step Size Control for
#   Compressible Computational Fluid Dynamics
#   [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
@cache struct LowStorageRK3SpCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k::rateType
  utilde::uType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

struct LowStorageRK3SpConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  γ12end::SVector{N,T} # γ11 is always zero
  γ22end::SVector{N,T} # γ21 is always one
  γ32end::SVector{N,T} # γ31 is always zero
  # TODO: γ302 == γ303 == 0 in all emthods implemented below -> possible optimization?
  δ2end ::SVector{N,T} # δ1  is always one
  β1::T
  β2end::SVector{N,T}
  c2end::SVector{N,T2} # c1 is always zero
  bhat1::T
  bhat2end::SVector{N,T}
end


function RDPK3Sp35ConstantCache(T, T2)
  γ12end = SVector(
    convert(T, big"2.587669070352079020144955303389306026e-01"),
    convert(T, big"-1.324366873994502973977035353758550057e-01"),
    convert(T, big"5.055601231460399101814291350373559483e-02"),
    convert(T, big"5.670552807902877312521811889846000976e-01"),
  )

  γ22end = SVector(
    convert(T, big"5.528418745102160639901976698795928733e-01"),
    convert(T, big"6.731844400389673824374042790213570079e-01"),
    convert(T, big"2.803103804507635075215805236096803381e-01"),
    convert(T, big"5.521508873507393276457754945308880998e-01"),
  )

  γ32end = SVector(
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"2.752585813446636957256614568573008811e-01"),
    convert(T, big"-8.950548709279785077579454232514633376e-01"),
  )

  δ2end = SVector(
    convert(T, big"3.407687209321455242558804921815861422e-01"),
    convert(T, big"3.414399280584625023244387687873774697e-01"),
    convert(T, big"7.229302732875589702087936723400941329e-01"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
  )

  β1 = convert(T, big"2.300285062878154351930669430512780706e-01")
  β2end = SVector(
    convert(T, big"3.021457892454169700189445968126242994e-01"),
    convert(T, big"8.025601039472704213300183888573974531e-01"),
    convert(T, big"4.362158997637629844305216319994356355e-01"),
    convert(T, big"1.129268494470295369172265188216779157e-01"),
  )

  c2end = SVector(
    convert(T, big"2.300285062878154351930669430512780706e-01"),
    convert(T, big"4.050049049262914975700372321130661410e-01"),
    convert(T, big"8.947823877926760224705450466361360720e-01"),
    convert(T, big"7.235108137218888081489570284485201518e-01"),
  )

  # difference of the usual bhat coefficients and the main b coefficients
  bhat1    = convert(T, big"1.046363371354093758897668305991705199e-01"
                      - big"1.147931563369900682037379182772608287e-01")
  bhat2end = SVector(
    convert(T, big"9.520431574956758809511173383346476348e-02"
             - big"8.933559295232859013880114997436974196e-02"),
    convert(T, big"4.482446645568668405072421350300379357e-01"
             - big"4.355858717379231779899161991033964256e-01"),
    convert(T, big"2.449030295461310135957132640369862245e-01"
             - big"2.473585295257286267503182138232950881e-01"),
    convert(T, big"1.070116530120251819121660365003405564e-01"
             - big"1.129268494470295369172265188216779157e-01"),
  )

  LowStorageRK3SpConstantCache{4,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end)
end

function alg_cache(alg::RDPK3Sp35,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  utilde = zero(u)
  tmp = zero(u)
  if eltype(u) === uEltypeNoUnits
    atmp = utilde # alias the vectors to save memory
  else
    atmp = similar(u,uEltypeNoUnits)
  end
  tab = RDPK3Sp35ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SpCache(u,uprev,fsalfirst,k,utilde,tmp,atmp,tab)
end

function alg_cache(alg::RDPK3Sp35,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  RDPK3Sp35ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function RDPK3Sp49ConstantCache(T, T2)
  γ12end = SVector(
    convert(T, big"-4.655641301259180308677051498071354582e+00"),
    convert(T, big"-7.720264924836063859141482018013692338e-01"),
    convert(T, big"-4.024423213419724605695005429153112050e+00"),
    convert(T, big"-2.129685246739018613087466942802498152e-02"),
    convert(T, big"-2.435022519234470128602335652131234586e+00"),
    convert(T, big"1.985627480986167686791439120784668251e-02"),
    convert(T, big"-2.810790112885283952929218377438668784e-01"),
    convert(T, big"1.689434895835535695524003319503844110e-01"),
  )

  γ22end = SVector(
    convert(T, big"2.499262752607825957145627300817258023e+00"),
    convert(T, big"5.866820365436136799319929406678132638e-01"),
    convert(T, big"1.205141365412670762568835277881144391e+00"),
    convert(T, big"3.474793796700868848597960521248007941e-01"),
    convert(T, big"1.321346140128723105871355808477092220e+00"),
    convert(T, big"3.119636324379370564023292317172847140e-01"),
    convert(T, big"4.351419055894087609560896967082486864e-01"),
    convert(T, big"2.359698299440788299161958168555704234e-01"),
  )

  γ32end = SVector(
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"7.621037111138170045618771082985664430e-01"),
    convert(T, big"-1.981182159087218433914909510116664154e-01"),
    convert(T, big"-6.228960706317566993192689455719570179e-01"),
    convert(T, big"-3.752246993432626328289874575355102038e-01"),
    convert(T, big"-3.355436539000946543242869676125143358e-01"),
    convert(T, big"-4.560963110717484359015342341157302403e-02"),
  )

  δ2end = SVector(
    convert(T, big"1.262923854387806460989545005598562667e+00"),
    convert(T, big"7.574967177560872438940839460448329992e-01"),
    convert(T, big"5.163591158111222863455531895152351544e-01"),
    convert(T, big"-2.746333792042827389548936599648122146e-02"),
    convert(T, big"-4.382674653941770848797864513655752318e-01"),
    convert(T, big"1.273587103668392811985704533534301656e+00"),
    convert(T, big"-6.294740045442794829622796613103492913e-01"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
  )

  β1 = convert(T, big"2.836343531977826022543660465926414772e-01")
  β2end = SVector(
    convert(T, big"9.736497978646965372894268287659773644e-01"),
    convert(T, big"3.382358566377620380505126936670933370e-01"),
    convert(T, big"-3.584937820217850715182820651063453804e-01"),
    convert(T, big"-4.113955814725134294322006403954822487e-03"),
    convert(T, big"1.427968962196019024010757034274849198e+00"),
    convert(T, big"1.808467712038743032991177525728915926e-02"),
    convert(T, big"1.605771316794521018947553625079465692e-01"),
    convert(T, big"2.952226811394310028003810072027839487e-01"),
  )

  c2end = SVector(
    convert(T, big"2.836343531977826022543660465926414772e-01"),
    convert(T, big"5.484073767552486705240014599676811834e-01"),
    convert(T, big"3.687229456675706936558667052479014150e-01"),
    convert(T, big"-6.806119916032093175251948474173648331e-01"),
    convert(T, big"3.518526451892056368706593492732753284e-01"),
    convert(T, big"1.665941920204672094647868254892387293e+00"),
    convert(T, big"9.715276989307335935187466054546761665e-01"),
    convert(T, big"9.051569554420045339601721625247585643e-01"),
  )

  # difference of the usual bhat coefficients and the main b coefficients
  bhat1    = convert(T, big"4.550655927970944948340364817140593012e-02"
                      - big"4.503731969165884304041981629148469971e-02")
  bhat2end = SVector(
    convert(T, big"1.175968310492638562142460384341959193e-01"
             - big"1.859217322011968812563859888433403777e-01"),
    convert(T, big"3.658257330515213200375475084421083608e-02"
             - big"3.329727509207630932171676116314110008e-02"),
    convert(T, big"-5.311555834355629559010061596928357525e-03"
             - big"-4.784222621050198909820741390895649698e-03"),
    convert(T, big"5.178250012713127329531367677410650996e-03"
             - big"4.055848062637567925908043629915811671e-03"),
    convert(T, big"4.954639022118682638697706200022961443e-01"
             - big"4.185027999682794463309031355073933444e-01"),
    convert(T, big"-5.999303132737865921441409466809521699e-03"
             - big"-4.381894507474277848407591859322000026e-03"),
    convert(T, big"9.405093434568315929035250835218733824e-02"
             - big"2.712846097324442608251358061215836749e-02"),
    convert(T, big"2.169318087627035072893925375820310602e-01"
             - big"2.952226811394310028003810072027839487e-01"),
  )

  LowStorageRK3SpConstantCache{8,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end)
end

function alg_cache(alg::RDPK3Sp49,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  utilde = zero(u)
  tmp = zero(u)
  if eltype(u) === uEltypeNoUnits
    atmp = utilde # alias the vectors to save memory
  else
    atmp = similar(u,uEltypeNoUnits)
  end
  tab = RDPK3Sp49ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SpCache(u,uprev,fsalfirst,k,utilde,tmp,atmp,tab)
end

function alg_cache(alg::RDPK3Sp49,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  RDPK3Sp49ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function RDPK3Sp510ConstantCache(T, T2)
  γ12end = SVector(
    convert(T, big"4.043660078504695837542588769963326988e-01"),
    convert(T, big"-8.503427464263185087039788184485627962e-01"),
    convert(T, big"-6.950894167072419998080989313353063399e+00"),
    convert(T, big"9.238765225328278557805080247596562995e-01"),
    convert(T, big"-2.563178039957404359875124580586147888e+00"),
    convert(T, big"2.545744869966347362604059848503340890e-01"),
    convert(T, big"3.125831733863168874151935287174374515e-01"),
    convert(T, big"-7.007114800567584871263283872289072079e-01"),
    convert(T, big"4.839620970980726631935174740648996010e-01"),
  )

  γ22end = SVector(
    convert(T, big"6.871467069752345566001768382316915820e-01"),
    convert(T, big"1.093024760468898686510433898645775908e+00"),
    convert(T, big"3.225975382330161123625348062949430509e+00"),
    convert(T, big"1.041153700841396427100436517666787823e+00"),
    convert(T, big"1.292821488864702752767390075072674807e+00"),
    convert(T, big"7.391462769297006312785029455392854586e-01"),
    convert(T, big"1.239129257039300081860496157739352186e-01"),
    convert(T, big"1.842753479366766790220633908793933781e-01"),
    convert(T, big"5.712788942697077644959290025755003720e-02"),
  )

  γ32end = SVector(
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"-2.393405159342139386425044844626597490e+00"),
    convert(T, big"-1.902854422095986544338294743445530533e+00"),
    convert(T, big"-2.820042210583207174321941694153843259e+00"),
    convert(T, big"-1.832698464130564949123807896975136336e+00"),
    convert(T, big"-2.199094510750697865007677774395365522e-01"),
    convert(T, big"-4.082430660384876496971887725512427800e-01"),
    convert(T, big"-1.377669791121207993339861855818881150e-01"),
  )

  δ2end = SVector(
    convert(T, big"-1.331778409133849616712007380176762548e-01"),
    convert(T, big"8.260422785246030254485064732649153253e-01"),
    convert(T, big"1.513700430513332405798616943654007796e+00"),
    convert(T, big"-1.305810063177048110528482211982726539e+00"),
    convert(T, big"3.036678789342507704281817524408221954e+00"),
    convert(T, big"-1.449458267074592489788800461540171106e+00"),
    convert(T, big"3.834313873320957483471400258279635203e+00"),
    convert(T, big"4.122293971923324492772059928094971199e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
  )

  β1 = convert(T, big"2.597883575710995826783320802193635406e-01")
  β2end = SVector(
    convert(T, big"1.777008800169541694837687556103565007e-02"),
    convert(T, big"2.481636637328140606807905234325691851e-01"),
    convert(T, big"7.941736827560429420202759490815682546e-01"),
    convert(T, big"3.885391296871822541486945325814526190e-01"),
    convert(T, big"1.455051664264339366757555740296587660e-01"),
    convert(T, big"1.587517379462528932413419955691782412e-01"),
    convert(T, big"1.650605631567659573994022720500446501e-01"),
    convert(T, big"2.118093299943235065178000892467421832e-01"),
    convert(T, big"1.559392340339606299335442956580114440e-01"),
  )

  c2end = SVector(
    convert(T, big"2.597883575710995826783320802193635406e-01"),
    convert(T, big"9.904573115730917688557891428202061598e-02"),
    convert(T, big"2.155511882303785204133426661931565216e-01"),
    convert(T, big"5.007950078421880417512789524851012021e-01"),
    convert(T, big"5.592251914858131230054392022144328176e-01"),
    convert(T, big"5.449986973408778242805929551952000165e-01"),
    convert(T, big"7.615224662599497796472095353126697300e-01"),
    convert(T, big"8.427062083059167761623893618875787414e-01"),
    convert(T, big"9.152209807185253394871325258038753352e-01"),
  )

  # difference of the usual bhat coefficients and the main b coefficients
  bhat1    = convert(T, big"5.734588484676193812418453938089759359e-02"
                      - big"-2.280102305596364773323878383881954511e-03")
  bhat2end = SVector(
    convert(T, big"1.971447518039733870541652912891291496e-02"
             - big"1.407393020823230537861040991952849386e-02"),
    convert(T, big"7.215296605683716720707226840456658773e-02"
             - big"2.332691794172822486743039657924919496e-01"),
    convert(T, big"1.739659489807939956977075317768151880e-01"
             - big"4.808266700465181307162297999657715930e-02"),
    convert(T, big"3.703693600445487815015171515640585668e-01"
             - big"4.119003221139622842134291677033040683e-01"),
    convert(T, big"-1.215599039055065009827765147821222534e-01"
             - big"-1.291461071364752805327361051196128312e-01"),
    convert(T, big"1.180372945491121604465067725859678821e-01"
             - big"1.220746011038579789984601943748468541e-01"),
    convert(T, big"4.155688823364870056536983972605056553e-02"
             - big"4.357858803113387764356338334851554715e-02"),
    convert(T, big"1.227886627910379901351569893551486490e-01"
             - big"1.025076875289905073925255867102192694e-01"),
    convert(T, big"1.456284232223684285998448928597043056e-01"
             - big"1.559392340339606299335442956580114440e-01"),
  )

  LowStorageRK3SpConstantCache{9,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end)
end

function alg_cache(alg::RDPK3Sp510,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  utilde = zero(u)
  tmp = zero(u)
  if eltype(u) === uEltypeNoUnits
    atmp = utilde # alias the vectors to save memory
  else
    atmp = similar(u,uEltypeNoUnits)
  end
  tab = RDPK3Sp510ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SpCache(u,uprev,fsalfirst,k,utilde,tmp,atmp,tab)
end

function alg_cache(alg::RDPK3Sp510,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  RDPK3Sp510ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


# 3S+ FSAL low storage methods: 3S methods adding another memory location for the embedded method (FSAL version)
# ## References
# - Ranocha, Dalcin, Parsani, Ketcheson (2021)
#   Optimized Runge-Kutta Methods with Automatic Step Size Control for
#   Compressible Computational Fluid Dynamics
#   [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
@cache struct LowStorageRK3SpFSALCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k::rateType
  utilde::uType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

struct LowStorageRK3SpFSALConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  γ12end::SVector{N,T} # γ11 is always zero
  γ22end::SVector{N,T} # γ21 is always one
  γ32end::SVector{N,T} # γ31 is always zero
  # TODO: γ302 == γ303 == 0 in all emthods implemented below -> possible optimization?
  δ2end ::SVector{N,T} # δ1  is always one
  β1::T
  β2end::SVector{N,T}
  c2end::SVector{N,T2} # c1 is always zero
  bhat1::T
  bhat2end::SVector{N,T}
  bhatfsal::T
end


function RDPK3SpFSAL35ConstantCache(T, T2)
  γ12end = SVector(
    convert(T, big"2.587771979725733308135192812685323706e-01"),
    convert(T, big"-1.324380360140723382965420909764953437e-01"),
    convert(T, big"5.056033948190826045833606441415585735e-02"),
    convert(T, big"5.670532000739313812633197158607642990e-01"),
  )

  γ22end = SVector(
    convert(T, big"5.528354909301389892439698870483746541e-01"),
    convert(T, big"6.731871608203061824849561782794643600e-01"),
    convert(T, big"2.803103963297672407841316576323901761e-01"),
    convert(T, big"5.521525447020610386070346724931300367e-01"),
  )

  γ32end = SVector(
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"2.752563273304676380891217287572780582e-01"),
    convert(T, big"-8.950526174674033822276061734289327568e-01"),
  )

  δ2end = SVector(
    convert(T, big"3.407655879334525365094815965895763636e-01"),
    convert(T, big"3.414382655003386206551709871126405331e-01"),
    convert(T, big"7.229275366787987419692007421895451953e-01"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
  )

  β1 = convert(T, big"2.300298624518076223899418286314123354e-01")
  β2end = SVector(
    convert(T, big"3.021434166948288809034402119555380003e-01"),
    convert(T, big"8.025606185416310937583009085873554681e-01"),
    convert(T, big"4.362158943603440930655148245148766471e-01"),
    convert(T, big"1.129272530455059129782111662594436580e-01"),
  )

  c2end = SVector(
    convert(T, big"2.300298624518076223899418286314123354e-01"),
    convert(T, big"4.050046072094990912268498160116125481e-01"),
    convert(T, big"8.947822893693433545220710894560512805e-01"),
    convert(T, big"7.235136928826589010272834603680114769e-01"),
  )

  # difference of the usual bhat coefficients and the main b coefficients
  bhat1    = convert(T, big"9.484166705035703392326247283838082847e-02"
                      - big"1.147935971023541171733601324486904546e-01")
  bhat2end = SVector(
    convert(T, big"1.726371339430353766966762629176676070e-01"
             - big"8.933442853113315592708384523126474636e-02"),
    convert(T, big"3.998243189084371024483169698618455770e-01"
             - big"4.355871025008616992483722693795608738e-01"),
    convert(T, big"1.718016807580178450618829007973835152e-01"
             - big"2.473576188201451146729725866810402672e-01"),
    convert(T, big"5.881914422155740300718268359027168467e-02"
             - big"1.129272530455059129782111662594436580e-01"),
  )
  bhatfsal = convert(T, big"1.020760551185952388626787099944507877e-01")

  LowStorageRK3SpFSALConstantCache{4,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal)
end

function alg_cache(alg::RDPK3SpFSAL35,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  utilde = zero(u)
  tmp = zero(u)
  if eltype(u) === uEltypeNoUnits
    atmp = utilde # alias the vectors to save memory
  else
    atmp = similar(u,uEltypeNoUnits)
  end
  tab = RDPK3SpFSAL35ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SpFSALCache(u,uprev,fsalfirst,k,utilde,tmp,atmp,tab)
end

function alg_cache(alg::RDPK3SpFSAL35,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  RDPK3SpFSAL35ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function RDPK3SpFSAL49ConstantCache(T, T2)
  γ12end = SVector(
    convert(T, big"-4.655641447335068552684422206224169103e+00"),
    convert(T, big"-7.720265099645871829248487209517314217e-01"),
    convert(T, big"-4.024436690519806086742256154738379161e+00"),
    convert(T, big"-2.129676284018530966221583708648634733e-02"),
    convert(T, big"-2.435022509790109546199372365866450709e+00"),
    convert(T, big"1.985627297131987000579523283542615256e-02"),
    convert(T, big"-2.810791146791038566946663374735713961e-01"),
    convert(T, big"1.689434168754859644351230590422137972e-01"),
  )

  γ22end = SVector(
    convert(T, big"2.499262792574495009336242992898153462e+00"),
    convert(T, big"5.866820377718875577451517985847920081e-01"),
    convert(T, big"1.205146086523094569925592464380295241e+00"),
    convert(T, big"3.474793722186732780030762737753849272e-01"),
    convert(T, big"1.321346060965113109321230804210670518e+00"),
    convert(T, big"3.119636464694193615946633676950358444e-01"),
    convert(T, big"4.351419539684379261368971206040518552e-01"),
    convert(T, big"2.359698130028753572503744518147537768e-01"),
  )

  γ32end = SVector(
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"7.621006678721315291614677352949377871e-01"),
    convert(T, big"-1.981182504339400567765766904309673119e-01"),
    convert(T, big"-6.228959218699007450469629366684127462e-01"),
    convert(T, big"-3.752248380775956442989480369774937099e-01"),
    convert(T, big"-3.355438309135169811915662336248989661e-01"),
    convert(T, big"-4.560955005031121479972862973705108039e-02"),
  )

  δ2end = SVector(
    convert(T, big"1.262923876648114432874834923838556100e+00"),
    convert(T, big"7.574967189685911558308119415539596711e-01"),
    convert(T, big"5.163589453140728104667573195005629833e-01"),
    convert(T, big"-2.746327421802609557034437892013640319e-02"),
    convert(T, big"-4.382673178127944142238606608356542890e-01"),
    convert(T, big"1.273587294602656522645691372699677063e+00"),
    convert(T, big"-6.294740283927400326554066998751383342e-01"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
  )

  β1 = convert(T, big"2.836343005184365275160654678626695428e-01")
  β2end = SVector(
    convert(T, big"9.736500104654741223716056170419660217e-01"),
    convert(T, big"3.382359225242515288768487569778320563e-01"),
    convert(T, big"-3.584943611106183357043212309791897386e-01"),
    convert(T, big"-4.113944068471528211627210454497620358e-03"),
    convert(T, big"1.427968894048586363415504654313371031e+00"),
    convert(T, big"1.808470948394314017665968411915568633e-02"),
    convert(T, big"1.605770645946802213926893453819236685e-01"),
    convert(T, big"2.952227015964591648775833803635147962e-01"),
  )

  c2end = SVector(
    convert(T, big"2.836343005184365275160654678626695428e-01"),
    convert(T, big"5.484076570002894365286665352032296535e-01"),
    convert(T, big"3.687228761669438493478872632332010073e-01"),
    convert(T, big"-6.806126440140844191258463830024463902e-01"),
    convert(T, big"3.518526124230705801739919476290327750e-01"),
    convert(T, big"1.665941994879593315477304663913129942e+00"),
    convert(T, big"9.715279295934715835299192116436237065e-01"),
    convert(T, big"9.051569840159589594903399929316959062e-01"),
  )

  # difference of the usual bhat coefficients and the main b coefficients
  bhat1    = convert(T, big"2.483675912451591196775756814283216443e-02"
                      - big"4.503732627263753698356970706617404465e-02")
  bhat2end = SVector(
    convert(T, big"1.866327774562103796990092260942180726e-01"
             - big"1.859217303699847950262276860012454333e-01"),
    convert(T, big"5.671080795936984495604436622517631183e-02"
             - big"3.329729672569717599759560403851202805e-02"),
    convert(T, big"-3.447695439149287702616943808570747099e-03"
             - big"-4.784204180958975587114459316829942677e-03"),
    convert(T, big"3.602245056516636472203469198006404016e-03"
             - big"4.055835961031310727671557609188874328e-03"),
    convert(T, big"4.545570622145088936800484247980581766e-01"
             - big"4.185027772596074197662616795629003544e-01"),
    convert(T, big"-2.434665289427612407531544765622888855e-04"
             - big"-4.381901968919326084347037216500072323e-03	"),
    convert(T, big"6.642755361103549971517945063138312147e-02"
             - big"2.712843796446089829255188189179448399e-02"),
    convert(T, big"1.613697079523505006226025497715177578e-01"
             - big"2.952227015964591648775833803635147962e-01"),
  )
  bhatfsal = convert(T, big"4.955424859358438183052504342394102722e-02")

  LowStorageRK3SpFSALConstantCache{8,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal)
end

function alg_cache(alg::RDPK3SpFSAL49,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  utilde = zero(u)
  tmp = zero(u)
  if eltype(u) === uEltypeNoUnits
    atmp = utilde # alias the vectors to save memory
  else
    atmp = similar(u,uEltypeNoUnits)
  end
  tab = RDPK3SpFSAL49ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SpFSALCache(u,uprev,fsalfirst,k,utilde,tmp,atmp,tab)
end

function alg_cache(alg::RDPK3SpFSAL49,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  RDPK3SpFSAL49ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function RDPK3SpFSAL510ConstantCache(T, T2)
  γ12end = SVector(
    convert(T, big"4.043660121685749695640462197806189975e-01"),
    convert(T, big"-8.503427289575839690883191973980814832e-01"),
    convert(T, big"-6.950894175262117526410215315179482885e+00"),
    convert(T, big"9.238765192731084931855438934978371889e-01"),
    convert(T, big"-2.563178056509891340215942413817786020e+00"),
    convert(T, big"2.545744879365226143946122067064118430e-01"),
    convert(T, big"3.125831707411998258746812355492206137e-01"),
    convert(T, big"-7.007114414440507927791249989236719346e-01"),
    convert(T, big"4.839621016023833375810172323297465039e-01"),
  )

  γ22end = SVector(
    convert(T, big"6.871467028161416909922221357014564412e-01"),
    convert(T, big"1.093024748914750833700799552463885117e+00"),
    convert(T, big"3.225975379607193001678365742708874597e+00"),
    convert(T, big"1.041153702510101386914019859778740444e+00"),
    convert(T, big"1.292821487912164945157744726076279306e+00"),
    convert(T, big"7.391462755788122847651304143259254381e-01"),
    convert(T, big"1.239129251371800313941948224441873274e-01"),
    convert(T, big"1.842753472370123193132193302369345580e-01"),
    convert(T, big"5.712788998796583446479387686662738843e-02"),
  )

  γ32end = SVector(
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
    convert(T, big"-2.393405133244194727221124311276648940e+00"),
    convert(T, big"-1.902854422421760920850597670305403139e+00"),
    convert(T, big"-2.820042207399977261483046412236557428e+00"),
    convert(T, big"-1.832698465277380999601896111079977378e+00"),
    convert(T, big"-2.199094483084671192328083958346519535e-01"),
    convert(T, big"-4.082430635847870963724591602173546218e-01"),
    convert(T, big"-1.377669797880289713535665985132703979e-01"),
  )

  δ2end = SVector(
    convert(T, big"-1.331778419508803397033287009506932673e-01"),
    convert(T, big"8.260422814750207498262063505871077303e-01"),
    convert(T, big"1.513700425755728332485300719652378197e+00"),
    convert(T, big"-1.305810059935023735972298885749903694e+00"),
    convert(T, big"3.036678802924163246003321318996156380e+00"),
    convert(T, big"-1.449458274398895177922690618003584514e+00"),
    convert(T, big"3.834313899176362315089976408899373409e+00"),
    convert(T, big"4.122293760012985409330881631526514714e+00"),
    convert(T, big"0.000000000000000000000000000000000000e+00"),
  )

  β1 = convert(T, big"2.597883554788674084039539165398464630e-01")
  β2end = SVector(
    convert(T, big"1.777008889438867858759149597539211023e-02"),
    convert(T, big"2.481636629715501931294746189266601496e-01"),
    convert(T, big"7.941736871152005775821844297293296135e-01"),
    convert(T, big"3.885391285642019129575902994397298066e-01"),
    convert(T, big"1.455051657916305055730603387469193768e-01"),
    convert(T, big"1.587517385964749337690916959584348979e-01"),
    convert(T, big"1.650605617880053419242434594242509601e-01"),
    convert(T, big"2.118093284937153836908655490906875007e-01"),
    convert(T, big"1.559392342362059886106995325687547506e-01"),
  )

  c2end = SVector(
    convert(T, big"2.597883554788674084039539165398464630e-01"),
    convert(T, big"9.904573247592460887087003212056568980e-02"),
    convert(T, big"2.155511890524058691860390281856497503e-01"),
    convert(T, big"5.007950088969676776844289399972611534e-01"),
    convert(T, big"5.592251911688643533787800688765883636e-01"),
    convert(T, big"5.449986978853637084972622392134732553e-01"),
    convert(T, big"7.615224694532590139829150720490417596e-01"),
    convert(T, big"8.427062083267360939805493320684741215e-01"),
    convert(T, big"9.152209805057669959657927210873423883e-01"),
  )

  # difference of the usual bhat coefficients and the main b coefficients
  bhat1    = convert(T, big"-2.019255440012066080909442770590267512e-02"
                      - big"-2.280100321836980811830528665041532799e-03")
  bhat2end = SVector(
    convert(T, big"2.737903480959184339932730854141598275e-02"
             - big"1.407393115790186300730580636032878435e-02"),
    convert(T, big"3.028818636145965534365173822296811090e-01"
             - big"2.332691775508456597719992034291118324e-01"),
    convert(T, big"-3.656843880622222190071445247906780540e-02"
             - big"4.808266741353862546318531020856621860e-02"),
    convert(T, big"3.982664774676767729863101188528827405e-01"
             - big"4.119003217706951892385733111000873172e-01"),
    convert(T, big"-5.715959421140685436681459970502471634e-02"
             - big"-1.291461067807736321056740833501596735e-01"),
    convert(T, big"9.849855103848558320961101178888983150e-02"
             - big"1.220746013848710098878384114422516148e-01"),
    convert(T, big"6.654601552456084978615342374581437947e-02"
             - big"4.357858583174420432201228508067333299e-02"),
    convert(T, big"9.073479542748112726465375642050504556e-02"
             - big"1.025076877568080726158907518254273554e-01"),
    convert(T, big"8.432289325330803924891866923939606351e-02"
             - big"1.559392342362059886106995325687547506e-01"),
  )
  bhatfsal = convert(T, big"4.529095628204896774513180907141004447e-02")

  LowStorageRK3SpFSALConstantCache{9,T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal)
end

function alg_cache(alg::RDPK3SpFSAL510,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  utilde = zero(u)
  tmp = zero(u)
  if eltype(u) === uEltypeNoUnits
    atmp = utilde # alias the vectors to save memory
  else
    atmp = similar(u,uEltypeNoUnits)
  end
  tab = RDPK3SpFSAL510ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3SpFSALCache(u,uprev,fsalfirst,k,utilde,tmp,atmp,tab)
end

function alg_cache(alg::RDPK3SpFSAL510,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  RDPK3SpFSAL510ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end



# 2R+ low storage methods introduced by van der Houwen
@cache struct LowStorageRK2RPCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  gprev::uType
  fsalfirst::rateType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

struct LowStorageRK2RPConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  Aᵢ::SVector{N,T}
  Bₗ::T
  B̂ₗ::T
  Bᵢ::SVector{N,T}
  B̂ᵢ::SVector{N,T}
  Cᵢ::SVector{N,T2}
end


function CKLLSRK43_2ConstantCache(T, T2)
  A1 = convert(T,Int128(11847461282814)//Int128(36547543011857))
  A2 = convert(T,Int128(3943225443063)//Int128(7078155732230))
  A3 = convert(T,Int128(-346793006927)//Int128(4029903576067))
  Aᵢ = SVector(A1, A2, A3)

  B1 = convert(T,Int128(1017324711453)//Int128(9774461848756))
  B2 = convert(T,Int128(8237718856693)//Int128(13685301971492))
  B3 = convert(T,Int128(57731312506979)//Int128(19404895981398))
  Bᵢ = SVector(B1, B2, B3)

  B̂1 = convert(T,Int128(15763415370699)//Int128(46270243929542))
  B̂2 = convert(T,Int128(514528521746)//Int128(5659431552419))
  B̂3 = convert(T,Int128(27030193851939)//Int128(9429696342944))
  B̂ᵢ = SVector(B̂1, B̂2, B̂3)

  Bₗ = convert(T,Int128(-101169746363290)//Int128(37734290219643))
  B̂ₗ = convert(T,Int128(-69544964788955)//Int128(30262026368149))

  C1 = convert(T2,Int128(11847461282814)//Int128(36547543011857))                                                  # A1
  C2 = convert(T2,Int128(2079258608735161403527719)//Int128(3144780143828896577027540))                            # A2 + B1
  C3 = convert(T2,Int128(41775191021672206476512620310545281003)//Int128(67383242951014563804622635478530729598))  # A3 + B1 + B2
  Cᵢ = SVector(C1, C2, C3)

  LowStorageRK2RPConstantCache{3,T,T2}(Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK43_2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK43_2ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK2RPCache(u,uprev,k,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK43_2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK43_2ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK54_3CConstantCache(T, T2)
  A1 = convert(T, BigInt(970286171893)//BigInt(4311952581923))
  A2 = convert(T, BigInt(6584761158862)//BigInt(12103376702013))
  A3 = convert(T, BigInt(2251764453980)//BigInt(15575788980749))
  A4 = convert(T, BigInt(26877169314380)//BigInt(34165994151039))
  Aᵢ = SVector(A1, A2, A3, A4)

  B1 = convert(T, BigInt(1153189308089)//BigInt(22510343858157))
  B2 = convert(T, BigInt(1772645290293)//BigInt(4653164025191))
  B3 = convert(T, BigInt(-1672844663538)//BigInt(4480602732383))
  B4 = convert(T, BigInt(2114624349019)//BigInt(3568978502595))
  Bᵢ = SVector(B1, B2, B3, B4)

  B̂1 = convert(T,BigInt(1016888040809)//BigInt(7410784769900))
  B̂2 = convert(T,BigInt(11231460423587)//BigInt(58533540763752))
  B̂3 = convert(T,BigInt(-1563879915014)//BigInt(6823010717585))
  B̂4 = convert(T,BigInt(606302364029)//BigInt(971179775848))
  B̂ᵢ = SVector(B̂1, B̂2, B̂3, B̂4)

  Bₗ = convert(T,BigInt(5198255086312)//BigInt(14908931495163))
  B̂ₗ = convert(T,BigInt(1097981568119)//BigInt(3980877426909))

  C1 = convert(T2,BigInt(970286171893)//BigInt(4311952581923))                                                                                  # A1
  C2 = convert(T2,BigInt(18020302501594987297224499)//BigInt(30272352378568762325374449))                                                       # A2 + B1
  C3 = convert(T2,BigInt(940957347754451928235896289983310398260)//BigInt(1631475460071027605339136597003329167263))                            # A3 + B1 + B2
  C4 = convert(T2,BigInt(8054848232572758807908657851968985615984276476412066)//BigInt(8139155613487734148190408375391604039319069461908135))   # A4 + B1 + B2 + B3
  Cᵢ = SVector(C1, C2, C3, C4)

  LowStorageRK2RPConstantCache{4,T,T2}(Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK54_3C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK54_3CConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK2RPCache(u,uprev,k,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK54_3C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK54_3CConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK95_4SConstantCache(T, T2)
    A1  = convert(T, BigInt(1107026461565)//BigInt(5417078080134))
  A2  = convert(T, BigInt(38141181049399)//BigInt(41724347789894))
  A3  = convert(T, BigInt(493273079041)//BigInt(11940823631197))
  A4  = convert(T, BigInt(1851571280403)//BigInt(6147804934346))
  A5  = convert(T, BigInt(11782306865191)//BigInt(62590030070788))
  A6  = convert(T, BigInt(9452544825720)//BigInt(13648368537481))
  A7  = convert(T, BigInt(4435885630781)//BigInt(26285702406235))
  A8  = convert(T, BigInt(2357909744247)//BigInt(11371140753790))
  Aᵢ = SVector(A1, A2, A3, A4, A5, A6, A7, A8)

  B1  = convert(T, BigInt(2274579626619)//BigInt(23610510767302))
  B2  = convert(T, BigInt(693987741272)//BigInt(12394497460941))
  B3  = convert(T, BigInt(-347131529483)//BigInt(15096185902911))
  B4  = convert(T, BigInt(1144057200723)//BigInt(32081666971178))
  B5  = convert(T, BigInt(1562491064753)//BigInt(11797114684756))
  B6  = convert(T, BigInt(13113619727965)//BigInt(44346030145118))
  B7  = convert(T, BigInt(393957816125)//BigInt(7825732611452))
  B8  = convert(T, BigInt(720647959663)//BigInt(6565743875477))
  Bᵢ = SVector(B1, B2, B3, B4, B5, B6, B7, B8)

  B̂1  = convert(T, BigInt(266888888871)//BigInt(3040372307578))
  B̂2  = convert(T, BigInt(34125631160)//BigInt(2973680843661))
  B̂3  = convert(T, BigInt(-653811289250)//BigInt(9267220972999))
  B̂4  = convert(T, BigInt(323544662297)//BigInt(2461529853637))
  B̂5  = convert(T, BigInt(1105885670474)//BigInt(4964345317203))
  B̂6  = convert(T, BigInt(1408484642121)//BigInt(8758221613943))
  B̂7  = convert(T, BigInt(1454774750537)//BigInt(11112645198328))
  B̂8  = convert(T, BigInt(772137014323)//BigInt(4386814405182))
  B̂ᵢ = SVector(B̂1, B̂2, B̂3, B̂4, B̂5, B̂6, B̂7, B̂8)

  Bₗ = convert(T, BigInt(3559252274877)//BigInt(14424734981077))
  B̂ₗ = convert(T, BigInt(277420604269)//BigInt(1857595682219))

  C1  = convert(T2, BigInt(1107026461565)//BigInt(5417078080134))                                                                                                                                                                                   # A1
  C2  = convert(T2, BigInt(248859529315327119359384971)//BigInt(246283290687986423455311497))                                                                                                                                                       # A2 + B1
  C3  = convert(T2, BigInt(676645811244741430568548054467096184193)//BigInt(3494367591912647069105975861901917224854))                                                                                                                              # A3 + B1 + B2
  C4  = convert(T2, BigInt(974370561662349106845723178377944301517533305964589)//BigInt(2263290880944514209862892217007179742168288737673791))                                                                                                      # A4 + B1 + B2 + B3
  C5  = convert(T2, BigInt(23738915426186839814576142955255044211724736499516359049188590711)//BigInt(67203160149331519751012175988216621571869262839903428488408759604))                                                                           # A5 + B1 + B2 + B3 + B4
  C6  = convert(T2, BigInt(1882683585832901544671586749377753597775777511029847145277760106172106584376955)//BigInt(1901663903553486696887572033100456166564493852721284994300276200102719954709068))                                               # A6 + B1 + B2 + B3 + B4 + B5
  C7  = convert(T2, BigInt(61872982955093233917984290421186995265732234396821660871734841970091372539489172106504162637)//BigInt(81207728164913218881758751120099941603350662788460257311895072645631357391473675997419584220))                     # A7 + B1 + B2 + B3 + B4 + B5 + B6
  C8  = convert(T2, BigInt(197565042693102647130189450792520184956129841555961940530192020871289515369046683661585184411130637357)//BigInt(232196202198018941876505157326935602816917261769279531369710269478309137067357703513986211472070374865)) # A8 + B1 + B2 + B3 + B4 + B5 + B6 + B7
  Cᵢ = SVector(C1, C2, C3, C4, C5, C6, C7, C8)

  LowStorageRK2RPConstantCache{8,T,T2}(Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK95_4S,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK95_4SConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK2RPCache(u,uprev,k,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK95_4S,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK95_4SConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK95_4CConstantCache(T, T2)
  A1  = convert(T, BigInt(2756167973529)//BigInt(16886029417639))
  A2  = convert(T, BigInt(11436141375279)//BigInt(13592993952163))
  A3  = convert(T, BigInt(88551658327)//BigInt(2352971381260))
  A4  = convert(T, BigInt(1882111988787)//BigInt(5590444193957))
  A5  = convert(T, BigInt(846820081679)//BigInt(4754706910573))
  A6  = convert(T, BigInt(4475289710031)//BigInt(6420120086209))
  A7  = convert(T, BigInt(118394748311)//BigInt(9144450320350))
  A8  = convert(T, BigInt(3307377157135)//BigInt(13111544596386))
  Aᵢ = SVector(A1, A2, A3, A4, A5, A6, A7, A8)

  B1  = convert(T, BigInt(1051460336009)//BigInt(14326298067773))
  B2  = convert(T, BigInt(930517604889)//BigInt(7067438519321))
  B3  = convert(T, BigInt(-311910530565)//BigInt(11769786407153))
  B4  = convert(T, BigInt(-410144036239)//BigInt(7045999268647))
  B5  = convert(T, BigInt(16692278975653)//BigInt(83604524739127))
  B6  = convert(T, BigInt(3777666801280)//BigInt(13181243438959))
  B7  = convert(T, BigInt(286682614203)//BigInt(12966190094317))
  B8  = convert(T, BigInt(3296161604512)//BigInt(22629905347183))
  Bᵢ = SVector(B1, B2, B3, B4, B5, B6, B7, B8)

  B̂1  = convert(T, BigInt(3189770262221)//BigInt(35077884776239))
  B̂2  = convert(T, BigInt(780043871774)//BigInt(11919681558467))
  B̂3  = convert(T, BigInt(-483824475979)//BigInt(5387739450692))
  B̂4  = convert(T, BigInt(1306553327038)//BigInt(9528955984871))
  B̂5  = convert(T, BigInt(6521106697498)//BigInt(22565577506855))
  B̂6  = convert(T, BigInt(1400555694605)//BigInt(19784728594468))
  B̂7  = convert(T, BigInt(1183541508418)//BigInt(13436305181271))
  B̂8  = convert(T, BigInt(3036254792728)//BigInt(15493572606329))
  B̂ᵢ = SVector(B̂1, B̂2, B̂3, B̂4, B̂5, B̂6, B̂7, B̂8)

  Bₗ = convert(T, BigInt(2993490409874)//BigInt(13266828321767))
  B̂ₗ = convert(T, BigInt(638483435745)//BigInt(4187244659458))

  C1  = convert(T2, BigInt(2756167973529)//BigInt(16886029417639))                                                                                                                                                                                          # A1
  C2  = convert(T2, BigInt(178130064075748009421121134)//BigInt(194737282992122861693942999))                                                                                                                                                               # A2 + B1
  C3  = convert(T2, BigInt(57818276708998807530478158133449099851)//BigInt(238238895426494403638887583424360627580))                                                                                                                                        # A3 + B1 + B2
  C4  = convert(T2, BigInt(3432454166457135667348375590572529790194124848059104)//BigInt(6662096512485931545803670383440459769502981926779993))                                                                                                             # A4 + B1 + B2 + B3
  C5  = convert(T2, BigInt(11915126765643872062053118401193741919814944004335534493046474237)//BigInt(39923715169802034300462756237193519081954994679332637422466438119))                                                                                   # A5 + B1 + B2 + B3 + B4
  C6  = convert(T2, BigInt(4583883621300589683158355859163890943947800555246686854224916208836514024614442)//BigInt(4506922925096139856045533451931734406235454975594364558624038359246205017801029))                                                       # A6 + B1 + B2 + B3 + B4 + B5
  C7  = convert(T2, BigInt(52423219056629312880725209686636192777075511202228566787042655312097949192300218484424118619)//BigInt(84615702680158836756876794083943762639542619835321175569533203672153042594634924742431352650))                             # A7 + B1 + B2 + B3 + B4 + B5 + B6
  C8  = convert(T2, BigInt(1385843715228499555828057735261132084759031703937678116167963792224108372724503731226480538087331079769069)//BigInt(1573111845759510782008384284066606688388217112071821912231287750254246452350240904652428530379336814559998)) # A8 + B1 + B2 + B3 + B4 + B5 + B6 + B7
  Cᵢ = SVector(C1, C2, C3, C4, C5, C6, C7, C8)

  LowStorageRK2RPConstantCache{8,T,T2}(Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK95_4C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK95_4CConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK2RPCache(u,uprev,k,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK95_4C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK95_4CConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK95_4MConstantCache(T, T2)
  A1  = convert(T, BigInt(5573095071601)//BigInt(11304125995793))
  A2  = convert(T, BigInt(315581365608)//BigInt(4729744040249))
  A3  = convert(T, BigInt(8734064225157)//BigInt(30508564569118))
  A4  = convert(T, BigInt(6457785058448)//BigInt(14982850401353))
  A5  = convert(T, BigInt(5771559441664)//BigInt(18187997215013))
  A6  = convert(T, BigInt(1906712129266)//BigInt(6681214991155))
  A7  = convert(T, BigInt(311585568784)//BigInt(2369973437185))
  A8  = convert(T, BigInt(-4840285693886)//BigInt(7758383361725))
  Aᵢ = SVector(A1, A2, A3, A4, A5, A6, A7, A8)

  B1  = convert(T, BigInt(549666665015)//BigInt(5899839355879))
  B2  = convert(T, BigInt(-548816778320)//BigInt(9402908589133))
  B3  = convert(T, BigInt(1672704946363)//BigInt(13015471661974))
  B4  = convert(T, BigInt(1025420337373)//BigInt(5970204766762))
  B5  = convert(T, BigInt(1524419752016)//BigInt(6755273790179))
  B6  = convert(T, BigInt(-10259399787359)//BigInt(43440802207630))
  B7  = convert(T, BigInt(4242280279850)//BigInt(10722460893763))
  B8  = convert(T, BigInt(1887552771913)//BigInt(6099058196803))
  Bᵢ = SVector(B1, B2, B3, B4, B5, B6, B7, B8)

  B̂1  = convert(T, BigInt(330911065672)//BigInt(9937126492277))
  B̂2  = convert(T, BigInt(-872991930418)//BigInt(11147305689291))
  B̂3  = convert(T, BigInt(2575378033706)//BigInt(14439313202205))
  B̂4  = convert(T, BigInt(3046892121673)//BigInt(11013392356255))
  B̂5  = convert(T, BigInt(1780184658016)//BigInt(8929499316295))
  B̂6  = convert(T, BigInt(10265149063)//BigInt(2098741126425))
  B̂7  = convert(T, BigInt(1643090076625)//BigInt(4891294770654))
  B̂8  = convert(T, BigInt(116106750067)//BigInt(3955800826265))
  B̂ᵢ = SVector(B̂1, B̂2, B̂3, B̂4, B̂5, B̂6, B̂7, B̂8)

  Bₗ = convert(T, BigInt(-453873186647)//BigInt(15285235680030))
  B̂ₗ = convert(T, BigInt(866868642257)//BigInt(42331321870877))

  C1  = convert(T2, BigInt(5573095071601)//BigInt(11304125995793))
  C2  = convert(T2, BigInt(4461661993774357683398167)//BigInt(27904730031895199210773871))
  C3  = convert(T2, BigInt(543425730194107827015264404954831354769)//BigInt(1692482454734045499140692116457071506026))
  C4  = convert(T2, BigInt(6429586327013850295560537918723231687699697140756067)//BigInt(10818243561353065593628044468492745774799533452459554))
  C5  = convert(T2, BigInt(555984804780268998022260997164198311752115182012221553157164786)//BigInt(852213854337283773231630192518719827415190771786411558523853399))
  C6  = convert(T2, BigInt(1789345671284476461332539715762783748132668223013904373945129499237446392572)//BigInt(2114764997945705573761804541148983827155257005191540481884326639410208291635))
  C7  = convert(T2, BigInt(2972211964132922642906704796208250552795647483819924111704054115070043529037601892705217)//BigInt(6517454043294174770082798998332814729652497865130816822916618330047242844192616374937270))
  C8  = convert(T2, BigInt(22038106775746116973750004935225594022265950105933360206617843987546593773108577078867914238620973639)//BigInt(228770596964454885481304478061363897900267080665965044117230250287302271092811814450282133504194141850))
  Cᵢ = SVector(C1, C2, C3, C4, C5, C6, C7, C8)

  LowStorageRK2RPConstantCache{8,T,T2}(Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK95_4M,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK95_4MConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK2RPCache(u,uprev,k,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK95_4M,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK95_4MConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end



# 3R+ low storage methods introduced by van der Houwen
@cache struct LowStorageRK3RPCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  fᵢ₋₂::rateType
  gprev::uType
  fsalfirst::rateType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

struct LowStorageRK3RPConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  Aᵢ₁::SVector{N,T}
  Aᵢ₂::SVector{N,T}
  Bₗ::T
  B̂ₗ::T
  Bᵢ::SVector{N,T}
  B̂ᵢ::SVector{N,T}
  Cᵢ::SVector{N,T2}
end


function CKLLSRK54_3C_3RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(2365592473904)//BigInt(8146167614645))
  A₁2  = convert(T, BigInt(4278267785271)//BigInt(6823155464066))
  A₁3  = convert(T, BigInt(2789585899612)//BigInt(8986505720531))
  A₁4  = convert(T, BigInt(15310836689591)//BigInt(24358012670437))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(-722262345248)//BigInt(10870640012513))
  A₂3  = convert(T, BigInt(1365858020701)//BigInt(8494387045469))
  A₂4  = convert(T, BigInt(3819021186)//BigInt(2763618202291))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4)

  B1   = convert(T, BigInt(846876320697)//BigInt(6523801458457))
  B2   = convert(T, BigInt(3032295699695)//BigInt(12397907741132))
  B3   = convert(T, BigInt(612618101729)//BigInt(6534652265123))
  B4   = convert(T, BigInt(1155491934595)//BigInt(2954287928812))
  Bᵢ   = SVector(B1,B2,B3,B4)

  B̂1   = convert(T, BigInt(1296459667021)//BigInt(9516889378644))
  B̂2   = convert(T, BigInt(2599004989233)//BigInt(11990680747819))
  B̂3   = convert(T, BigInt(1882083615375)//BigInt(8481715831096))
  B̂4   = convert(T, BigInt(1577862909606)//BigInt(5567358792761))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4)

  Bₗ   = convert(T, BigInt(707644755468)//BigInt(5028292464395))
  B̂ₗ   = convert(T, BigInt(328334985361)//BigInt(2316973589007))

  C1   = convert(T2, BigInt(2365592473904)//BigInt(8146167614645))
  C2   = convert(T2, BigInt(41579400703344293287237655)//BigInt(74172066799272566561857858))
  C3   = convert(T2, BigInt(299308060739053880467044545349561265546)//BigInt(497993456493513966629488516767096447823))
  C4   = convert(T2, BigInt(5468330126750791548369684419304733938034170906513585)//BigInt(5444638279732761024893610553331663911104849888809108))
  Cᵢ   = SVector(C1,C2,C3,C4)

  LowStorageRK3RPConstantCache{4,T,T2}(Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK54_3C_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK54_3C_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,fᵢ₋₂,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK54_3C_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK54_3C_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK54_3M_3RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(17396840518954)//BigInt(49788467287365))
  A₁2  = convert(T, BigInt(21253110367599)//BigInt(14558944785238))
  A₁3  = convert(T, BigInt(4293647616769)//BigInt(14519312872408))
  A₁4  = convert(T, BigInt(-8941886866937)//BigInt(7464816931160))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(-12587430488023)//BigInt(11977319897242))
  A₂3  = convert(T, BigInt(6191878339181)//BigInt(13848262311063))
  A₂4  = convert(T, BigInt(19121624165801)//BigInt(12321025968027))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4)

  B1   = convert(T, BigInt(1977388745448)//BigInt(17714523675943))
  B2   = convert(T, BigInt(6528140725453)//BigInt(14879534818174))
  B3   = convert(T, BigInt(4395900531415)//BigInt(55649460397719))
  B4   = convert(T, BigInt(6567440254656)//BigInt(15757960182571))
  Bᵢ   = SVector(B1,B2,B3,B4)

  B̂1   = convert(T, BigInt(390601394181)//BigInt(3503051559916))
  B̂2   = convert(T, BigInt(31150720071161)//BigInt(68604711794052))
  B̂3   = convert(T, BigInt(416927665232)//BigInt(6953044279741))
  B̂4   = convert(T, BigInt(3879867616328)//BigInt(8869216637007))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4)

  Bₗ   = convert(T, BigInt(-436008689643)//BigInt(9453681332953))
  B̂ₗ   = convert(T, BigInt(-163749046041)//BigInt(2599987820560))

  C1  = convert(T2, BigInt(17396840518954)//BigInt(49788467287365))
  C2  = convert(T2, BigInt(2546271293606266795002053)//BigInt(6227754966395669782804057))
  C3  = convert(T2, BigInt(3043453778831534771251734214272440269577)//BigInt(3561810617861654942925591050154818470872))
  C4  = convert(T2, BigInt(10963106193663894855575270257133723083246622141340761)//BigInt(12121458300971454511596914396147459030814063072954120))
  Cᵢ   = SVector(C1,C2,C3,C4)

  LowStorageRK3RPConstantCache{4,T,T2}(Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK54_3M_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK54_3M_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,fᵢ₋₂,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK54_3M_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK54_3M_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK54_3N_3RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(4745337637855)//BigInt(22386579876409))
  A₁2  = convert(T, BigInt(6808157035527)//BigInt(13197844641179))
  A₁3  = convert(T, BigInt(4367509502613)//BigInt(10454198590847))
  A₁4  = convert(T, BigInt(1236962429870)//BigInt(3429868089329))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(546509042554)//BigInt(9152262712923))
  A₂3  = convert(T, BigInt(625707605167)//BigInt(5316659119056))
  A₂4  = convert(T, BigInt(582400652113)//BigInt(7078426004906))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4)

  B1  = convert(T, BigInt(314199625218)//BigInt(7198350928319))
  B2  = convert(T, BigInt(6410344372641)//BigInt(17000082738695))
  B3  = convert(T, BigInt(292278564125)//BigInt(5593752632744))
  B4  = convert(T, BigInt(5010207514426)//BigInt(21876007855139))
  Bᵢ  = SVector(B1,B2,B3,B4)

  B̂1  = convert(T, BigInt(1276689330531)//BigInt(10575835502045))
  B̂2  = convert(T, BigInt(267542835879)//BigInt(1241767155676))
  B̂3  = convert(T, BigInt(1564039648689)//BigInt(9024646069760))
  B̂4  = convert(T, BigInt(3243722451631)//BigInt(13364844673806))
  B̂ᵢ  = SVector(B̂1,B̂2,B̂3,B̂4)

  Bₗ  = convert(T, BigInt(5597675544274)//BigInt(18784428342765))
  B̂ₗ  = convert(T, BigInt(606464709716)//BigInt(2447238536635))

  C1  = convert(T2, BigInt(4745337637855)//BigInt(22386579876409))
  C2  = convert(T2, BigInt(6320253019873211389522417)//BigInt(10980921945492108365568747))
  C3  = convert(T2, BigInt(231699760563456147635097088564862719039)//BigInt(400094496217566390613617613962197753808))
  C4  = convert(T2, BigInt(2565873674791335200443549967376635530873909687156071)//BigInt(2970969302106648098855751120425897741072516011514170))
  Cᵢ  = SVector(C1,C2,C3,C4)

  LowStorageRK3RPConstantCache{4,T,T2}(Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK54_3N_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK54_3N_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,fᵢ₋₂,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK54_3N_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK54_3N_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK85_4C_3RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(141236061735)//BigInt(3636543850841))
  A₁2  = convert(T, BigInt(7367658691349)//BigInt(25881828075080))
  A₁3  = convert(T, BigInt(6185269491390)//BigInt(13597512850793))
  A₁4  = convert(T, BigInt(2669739616339)//BigInt(18583622645114))
  A₁5  = convert(T, BigInt(42158992267337)//BigInt(9664249073111))
  A₁6  = convert(T, BigInt(970532350048)//BigInt(4459675494195))
  A₁7  = convert(T, BigInt(1415616989537)//BigInt(7108576874996))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4,A₁5,A₁6,A₁7)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(-343061178215)//BigInt(2523150225462))
  A₂3  = convert(T, BigInt(-4057757969325)//BigInt(18246604264081))
  A₂4  = convert(T, BigInt(1415180642415)//BigInt(13311741862438))
  A₂5  = convert(T, BigInt(-93461894168145)//BigInt(25333855312294))
  A₂6  = convert(T, BigInt(7285104933991)//BigInt(14106269434317))
  A₂7  = convert(T, BigInt(-4825949463597)//BigInt(16828400578907))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4,A₂5,A₂6,A₂7)

  B1  = convert(T, BigInt(514862045033)//BigInt(4637360145389))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(0)//BigInt(1))
  B4  = convert(T, BigInt(0)//BigInt(1))
  B5  = convert(T, BigInt(2561084526938)//BigInt(7959061818733))
  B6  = convert(T, BigInt(4857652849)//BigInt(7350455163355))
  B7  = convert(T, BigInt(1059943012790)//BigInt(2822036905401))
  Bᵢ  = SVector(B1,B2,B3,B4,B5,B6,B7)

  B̂1  = convert(T, BigInt(1269299456316)//BigInt(16631323494719))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(2153976949307)//BigInt(22364028786708))
  B̂4  = convert(T, BigInt(2303038467735)//BigInt(18680122447354))
  B̂5  = convert(T, BigInt(7354111305649)//BigInt(15643939971922))
  B̂6  = convert(T, BigInt(768474111281)//BigInt(10081205039574))
  B̂7  = convert(T, BigInt(3439095334143)//BigInt(10786306938509))
  B̂ᵢ  = SVector(B̂1,B̂2,B̂3,B̂4,B̂5,B̂6,B̂7)

  Bₗ  = convert(T, BigInt(2987336121747)//BigInt(15645656703944))
  B̂ₗ  = convert(T, BigInt(-3808726110015)//BigInt(23644487528593))

  C1  = convert(T2, BigInt(141236061735)//BigInt(3636543850841))
  C2  = convert(T2, BigInt(4855329627204641469273019)//BigInt(32651870171503411731843480))
  C3  = convert(T2, BigInt(395246570619540395679764439681768625174)//BigInt(1150568172675067443707820382013045349637))
  C4  = convert(T2, BigInt(103533040647279909858308372897770021461)//BigInt(286797987459862321650077169609703051387))
  C5  = convert(T2, BigInt(890342029406775514852349518244920625309)//BigInt(1135377348321966192554675673174478190626))
  C6  = convert(T2, BigInt(82180664649829640456237722943611531408)//BigInt(97244490215364259564723087293866304345))
  C7  = convert(T2, BigInt(1524044277359326675923410465291452002169116939509651)//BigInt(4415279581486844959297591640758696961331751174567964))
  Cᵢ  = SVector(C1,C2,C3,C4,C5,C6,C7)

  LowStorageRK3RPConstantCache{7,T,T2}(Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK85_4C_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK85_4C_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,fᵢ₋₂,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK85_4C_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK85_4C_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK85_4M_3RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(967290102210)//BigInt(6283494269639))
  A₁2  = convert(T, BigInt(852959821520)//BigInt(5603806251467))
  A₁3  = convert(T, BigInt(8043261511347)//BigInt(8583649637008))
  A₁4  = convert(T, BigInt(-115941139189)//BigInt(8015933834062))
  A₁5  = convert(T, BigInt(2151445634296)//BigInt(7749920058933))
  A₁6  = convert(T, BigInt(15619711431787)//BigInt(74684159414562))
  A₁7  = convert(T, BigInt(12444295717883)//BigInt(11188327299274))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4,A₁5,A₁6,A₁7)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(475331134681)//BigInt(7396070923784))
  A₂3  = convert(T, BigInt(-8677837986029)//BigInt(16519245648862))
  A₂4  = convert(T, BigInt(2224500752467)//BigInt(10812521810777))
  A₂5  = convert(T, BigInt(1245361422071)//BigInt(3717287139065))
  A₂6  = convert(T, BigInt(1652079198131)//BigInt(3788458824028))
  A₂7  = convert(T, BigInt(-5225103653628)//BigInt(8584162722535))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4,A₂5,A₂6,A₂7)

  B1  = convert(T, BigInt(83759458317)//BigInt(1018970565139))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(0)//BigInt(1))
  B4  = convert(T, BigInt(0)//BigInt(1))
  B5  = convert(T, BigInt(6968891091250)//BigInt(16855527649349))
  B6  = convert(T, BigInt(783521911849)//BigInt(8570887289572))
  B7  = convert(T, BigInt(3686104854613)//BigInt(11232032898210))
  Bᵢ  = SVector(B1,B2,B3,B4,B5,B6,B7)

  B̂1  = convert(T, BigInt(-2632078767757)//BigInt(9365288548818))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(138832778584802)//BigInt(30360463697573))
  B̂4  = convert(T, BigInt(7424139574315)//BigInt(5603229049946))
  B̂5  = convert(T, BigInt(-32993229351515)//BigInt(6883415042289))
  B̂6  = convert(T, BigInt(-3927384735361)//BigInt(7982454543710))
  B̂7  = convert(T, BigInt(9224293159931)//BigInt(15708162311543))
  B̂ᵢ  = SVector(B̂1,B̂2,B̂3,B̂4,B̂5,B̂6,B̂7)

  Bₗ  = convert(T, BigInt(517396786175)//BigInt(6104475356879))
  B̂ₗ  = convert(T, BigInt(624338737541)//BigInt(7691046757191))

  C1  = convert(T2, BigInt(967290102210)//BigInt(6283494269639))
  C2  = convert(T2, BigInt(8972214919142352493858707)//BigInt(41446148478994088895191128))
  C3  = convert(T2, BigInt(35682660731882055122214991891899678815)//BigInt(72242678055272695781813348615158920272))
  C4  = convert(T2, BigInt(24151963894889409757443700144610337197)//BigInt(88316684951621554188239538678367088186))
  C5  = convert(T2, BigInt(20396803294876689925555603189127802602)//BigInt(29355195069529377650856010387665377655))
  C6  = convert(T2, BigInt(104860372573190455963699691732496938387)//BigInt(144152676952392296448858925279884773652))
  C7  = convert(T2, BigInt(1648260218501227913212294426176971326433416596592133)//BigInt(1649556119556299790473636959153132604082083356090490))
  Cᵢ  = SVector(C1,C2,C3,C4,C5,C6,C7)

  LowStorageRK3RPConstantCache{7,T,T2}(Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK85_4M_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK85_4M_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,fᵢ₋₂,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK85_4M_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK85_4M_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK85_4P_3RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(1298271176151)//BigInt(60748409385661))
  A₁2  = convert(T, BigInt(14078610000243)//BigInt(41877490110127))
  A₁3  = convert(T, BigInt(553998884433)//BigInt(1150223130613))
  A₁4  = convert(T, BigInt(15658478150918)//BigInt(92423611770207))
  A₁5  = convert(T, BigInt(18843935397718)//BigInt(7227975568851))
  A₁6  = convert(T, BigInt(6206560082614)//BigInt(27846110321329))
  A₁7  = convert(T, BigInt(2841125392315)//BigInt(14844217636077))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4,A₁5,A₁6,A₁7)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(-2491873887327)//BigInt(11519757507826))
  A₂3  = convert(T, BigInt(-3833614938189)//BigInt(14183712281236))
  A₂4  = convert(T, BigInt(628609886693)//BigInt(8177399110319))
  A₂5  = convert(T, BigInt(-4943723744483)//BigInt(2558074780976))
  A₂6  = convert(T, BigInt(1024000837540)//BigInt(1998038638351))
  A₂7  = convert(T, BigInt(-2492809296391)//BigInt(9064568868273))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4,A₂5,A₂6,A₂7)

  B1  = convert(T, BigInt(346820227625)//BigInt(3124407780749))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(0)//BigInt(1))
  B4  = convert(T, BigInt(0)//BigInt(1))
  B5  = convert(T, BigInt(814249513470)//BigInt(2521483007009))
  B6  = convert(T, BigInt(195246859987)//BigInt(15831935944600))
  B7  = convert(T, BigInt(3570596951509)//BigInt(9788921605312))
  Bᵢ  = SVector(B1,B2,B3,B4,B5,B6,B7)

  B̂1  = convert(T, BigInt(679447319381)//BigInt(8240332772531))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(798472430005)//BigInt(13882421602211))
  B̂4  = convert(T, BigInt(972791992243)//BigInt(13597677393897))
  B̂5  = convert(T, BigInt(2994516937385)//BigInt(6097853295694))
  B̂6  = convert(T, BigInt(1424705874463)//BigInt(19211220871144))
  B̂7  = convert(T, BigInt(11199564863291)//BigInt(35136367926059))
  B̂ᵢ  = SVector(B̂1,B̂2,B̂3,B̂4,B̂5,B̂6,B̂7)

  Bₗ  = convert(T, BigInt(1886338382073)//BigInt(9981671730680))
  B̂ₗ  = convert(T, BigInt(-1307718103703)//BigInt(13694144003901))

  C1  = convert(T2, BigInt(1298271176151)//BigInt(60748409385661))
  C2  = convert(T2, BigInt(57828749177833338114741189)//BigInt(482418531105044571804353902))
  C3  = convert(T2, BigInt(16431909216114342992530887716659137419)//BigInt(50972944352640941110022041298448213332))
  C4  = convert(T2, BigInt(843711271601954807241466442429582743082)//BigInt(2361379786784371499429045948205315798717))
  C5  = convert(T2, BigInt(45377346645618697840609101263059649515)//BigInt(57769368855607143441437855651622233424))
  C6  = convert(T2, BigInt(147132600561369761792017800077859262701)//BigInt(173834563932749284125206995856250290771))
  C7  = convert(T2, BigInt(123785620236259768586332555932209432529705897037921)//BigInt(353351523019265026737831367789312912172448045683187))
  Cᵢ  = SVector(C1,C2,C3,C4,C5,C6,C7)

  LowStorageRK3RPConstantCache{7,T,T2}(Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK85_4P_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK85_4P_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK3RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,fᵢ₋₂,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK85_4P_3R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK85_4P_3RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


# 4R+ low storage methods introduced by van der Houwen
@cache struct LowStorageRK4RPCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  uᵢ₋₃::uType
  fᵢ₋₂::rateType
  fᵢ₋₃::rateType
  gprev::uType
  fsalfirst::rateType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

struct LowStorageRK4RPConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  Aᵢ₁::SVector{N,T}
  Aᵢ₂::SVector{N,T}
  Aᵢ₃::SVector{N,T}
  Bₗ::T
  B̂ₗ::T
  Bᵢ::SVector{N,T}
  B̂ᵢ::SVector{N,T}
  Cᵢ::SVector{N,T2}
end


function CKLLSRK54_3N_4RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(9435338793489)//BigInt(32856462503258))
  A₁2  = convert(T, BigInt(6195609865473)//BigInt(14441396468602))
  A₁3  = convert(T, BigInt(7502925572378)//BigInt(28098850972003))
  A₁4  = convert(T, BigInt(4527781290407)//BigInt(9280887680514))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(2934593324920)//BigInt(16923654741811))
  A₂3  = convert(T, BigInt(16352725096886)//BigInt(101421723321009))
  A₂4  = convert(T, BigInt(3004243580591)//BigInt(16385320447374))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4)

  A₃1  = convert(T, BigInt(0)//BigInt(1))
  A₃2  = convert(T, BigInt(0)//BigInt(1))
  A₃3  = convert(T, BigInt(390352446067)//BigInt(5989890148791))
  A₃4  = convert(T, BigInt(902830387041)//BigInt(8154716972155))
  Aᵢ₃  = SVector(A₃1,A₃2,A₃3,A₃4)

  B1  = convert(T, BigInt(929310922418)//BigInt(8329727308495))
  B2  = convert(T, BigInt(4343420149496)//BigInt(15735497610667))
  B3  = convert(T, BigInt(885252399220)//BigInt(9490460854667))
  B4  = convert(T, BigInt(3341719902227)//BigInt(13464012733180))
  Bᵢ  = SVector(B1,B2,B3,B4)

  B̂1  = convert(T, BigInt(2929323122013)//BigInt(17725327880387))
  B̂2  = convert(T, BigInt(4379799101587)//BigInt(35838171763617))
  B̂3  = convert(T, BigInt(2267325134734)//BigInt(9725002913543))
  B̂4  = convert(T, BigInt(1519467056643)//BigInt(5852430786130))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4)

  Bₗ   = convert(T, BigInt(2131913067577)//BigInt(7868783702050))
  B̂ₗ   = convert(T, BigInt(3636375423974)//BigInt(16547514622827))

  C1  = convert(T2, BigInt(9435338793489)//BigInt(32856462503258))
  C2  = convert(T2, BigInt(147231987957505837822553443)//BigInt(244401207824228867478118222))
  C3  = convert(T2, BigInt(401086457089554669663078760253749450489)//BigInt(812866282711293513804077001645679258017))
  C4  = convert(T2, BigInt(153823244836258719400905156342054669945035476219421)//BigInt(172160249040778711548900853819650745575758693592285))
  Cᵢ   = SVector(C1,C2,C3,C4)

  LowStorageRK4RPConstantCache{4,T,T2}(Aᵢ₁,Aᵢ₂,Aᵢ₃,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK54_3N_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  uᵢ₋₃ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  fᵢ₋₃ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK54_3N_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK4RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,fᵢ₋₂,fᵢ₋₃,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK54_3N_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK54_3N_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK54_3M_4RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(7142524119)//BigInt(20567653057))
  A₁2  = convert(T, BigInt(20567653057)//BigInt(89550000000))
  A₁3  = convert(T, BigInt(7407775)//BigInt(2008982))
  A₁4  = convert(T, BigInt(-4577300)//BigInt(867302297))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(15198616943)//BigInt(89550000000))
  A₂3  = convert(T, BigInt(-226244183627)//BigInt(80359280000))
  A₂4  = convert(T, BigInt(33311687500)//BigInt(8703531091))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4)

  A₃1  = convert(T, BigInt(0)//BigInt(1))
  A₃2  = convert(T, BigInt(0)//BigInt(1))
  A₃3  = convert(T, BigInt(9890667227)//BigInt(80359280000))
  A₃4  = convert(T, BigInt(-20567653057)//BigInt(6979191486))
  Aᵢ₃  = SVector(A₃1,A₃2,A₃3,A₃4)

  B1  = convert(T, BigInt(297809)//BigInt(2384418))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(156250000)//BigInt(270591503))
  B4  = convert(T, BigInt(5030000)//BigInt(888933))
  Bᵢ  = SVector(B1,B2,B3,B4)

  B̂1  = convert(T, BigInt(121286694859)//BigInt(931793198518))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(9680751416357)//BigInt(17201392077364))
  B̂4  = convert(T, BigInt(6633076090000)//BigInt(1042143269349))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4)

  Bₗ   = convert(T, BigInt(-2927)//BigInt(546))
  B̂ₗ   = convert(T, BigInt(-127961558623)//BigInt(21123456354))

  C1  = convert(T2, BigInt(7142524119)//BigInt(20567653057))
  C2  = convert(T2, BigInt(1997)//BigInt(5000))
  C3  = convert(T2, BigInt(199)//BigInt(200))
  C4  = convert(T2, BigInt(1)//BigInt(1))
  Cᵢ   = SVector(C1,C2,C3,C4)

  LowStorageRK4RPConstantCache{4,T,T2}(Aᵢ₁,Aᵢ₂,Aᵢ₃,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK54_3M_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  uᵢ₋₃ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  fᵢ₋₃ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK54_3M_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK4RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,fᵢ₋₂,fᵢ₋₃,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK54_3M_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK54_3M_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK65_4M_4RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(1811061732419)//BigInt(6538712036350))
  A₁2  = convert(T, BigInt(936386506953)//BigInt(6510757757683))
  A₁3  = convert(T, BigInt(8253430823511)//BigInt(9903985211908))
  A₁4  = convert(T, BigInt(4157325866175)//BigInt(11306150349782))
  A₁5  = convert(T, BigInt(3299942024581)//BigInt(13404534943033))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4,A₁5)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(968127049827)//BigInt(6993254963231))
  A₂3  = convert(T, BigInt(-4242729801665)//BigInt(12001587034923))
  A₂4  = convert(T, BigInt(1960956671631)//BigInt(3017447659538))
  A₂5  = convert(T, BigInt(2088737530132)//BigInt(14638867961951))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4,A₂5)

  A₃1  = convert(T, BigInt(0)//BigInt(1))
  A₃2  = convert(T, BigInt(0)//BigInt(1))
  A₃3  = convert(T, BigInt(332803037697)//BigInt(7529436905221))
  A₃4  = convert(T, BigInt(-19590089343957)//BigInt(51581831082203))
  A₃5  = convert(T, BigInt(3811366828049)//BigInt(10653298326636))
  Aᵢ₃  = SVector(A₃1,A₃2,A₃3,A₃4,A₃5)

  B1  = convert(T, BigInt(1437717300581)//BigInt(14622899446031))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(3070006287879)//BigInt(9321175678070))
  B4  = convert(T, BigInt(2276970273632)//BigInt(7940670647385))
  B5  = convert(T, BigInt(-1056149936631)//BigInt(7427907425983))
  Bᵢ  = SVector(B1,B2,B3,B4,B5)

  B̂1  = convert(T, BigInt(399352205828)//BigInt(2843676810815))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(460449895996)//BigInt(4301836608005))
  B̂4  = convert(T, BigInt(15965746118666)//BigInt(21690343195681))
  B̂5  = convert(T, BigInt(-19281717001664)//BigInt(29911607353389))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4,B̂5)

  Bₗ   = convert(T, BigInt(2571845656138)//BigInt(6012342010435))
  B̂ₗ   = convert(T, BigInt(5058427127221)//BigInt(7651806618075))

  C1  = convert(T2, BigInt(1811061732419)//BigInt(6538712036350))
  C2  = convert(T2, BigInt(12851630287335503073915984)//BigInt(45531389003311376172753773))
  C3  = convert(T2, BigInt(468994575306978457607500930904657513641)//BigInt(894975528626103930282351283769588361564))
  C4  = convert(T2, BigInt(4735520442856752193881763097298943558246492547269018)//BigInt(6433166018040288425494806218280078848936316641536447))
  C5  = convert(T2, BigInt(25828983228256103590265182981008154883102570637999497)//BigInt(30568689961801519095090666149791133914967119469889228))
  Cᵢ   = SVector(C1,C2,C3,C4,C5)

  LowStorageRK4RPConstantCache{5,T,T2}(Aᵢ₁,Aᵢ₂,Aᵢ₃,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK65_4M_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  uᵢ₋₃ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  fᵢ₋₃ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK65_4M_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK4RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,fᵢ₋₂,fᵢ₋₃,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK65_4M_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK65_4M_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


function CKLLSRK85_4FM_4RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(319960152914)//BigInt(39034091721739))
  A₁2  = convert(T, BigInt(16440040368765)//BigInt(7252463661539))
  A₁3  = convert(T, BigInt(1381950791880)//BigInt(6599155371617))
  A₁4  = convert(T, BigInt(18466735994895)//BigInt(7394178462407))
  A₁5  = convert(T, BigInt(2786140924985)//BigInt(14262827431161))
  A₁6  = convert(T, BigInt(28327099865656)//BigInt(21470840267743))
  A₁7  = convert(T, BigInt(0)//BigInt(1))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4,A₁5,A₁6,A₁7)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(-16195115415565)//BigInt(7808461210678))
  A₂3  = convert(T, BigInt(-1316066362688)//BigInt(10261382634081))
  A₂4  = convert(T, BigInt(-23893000145797)//BigInt(9614512377075))
  A₂5  = convert(T, BigInt(6556893593075)//BigInt(12530787773541))
  A₂6  = convert(T, BigInt(-5015572218207)//BigInt(5719938983072))
  A₂7  = convert(T, BigInt(0)//BigInt(1))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4,A₂5,A₂6,A₂7)

  A₃1  = convert(T, BigInt(0)//BigInt(1))
  A₃2  = convert(T, BigInt(0)//BigInt(1))
  A₃3  = convert(T, BigInt(334167490531)//BigInt(1677017272502))
  A₃4  = convert(T, BigInt(4579492417936)//BigInt(7930641522963))
  A₃5  = convert(T, BigInt(-2255846922213)//BigInt(30066310003000))
  A₃6  = convert(T, BigInt(3212719728776)//BigInt(7037340048693))
  A₃7  = convert(T, BigInt(0)//BigInt(1))
  Aᵢ₃  = SVector(A₃1,A₃2,A₃3,A₃4,A₃5,A₃6,A₃7)

  B1  = convert(T, BigInt(1147876221211)//BigInt(13910763665259))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(182134362610)//BigInt(9852075053293))
  B4  = convert(T, BigInt(3396705055007)//BigInt(8495597747463))
  B5  = convert(T, BigInt(363006049056)//BigInt(22366003978609))
  B6  = convert(T, BigInt(6078825123673)//BigInt(15200143133108))
  B7  = convert(T, BigInt(583593328277)//BigInt(7028929464160))
  Bᵢ  = SVector(B1,B2,B3,B4,B5,B6,B7)

  B̂1  = convert(T, BigInt(2023383632057)//BigInt(26525303340911))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(480990062147)//BigInt(12694528747923))
  B̂4  = convert(T, BigInt(14502014597821)//BigInt(36979005529861))
  B̂5  = convert(T, BigInt(-3883966523914)//BigInt(63014133260123))
  B̂6  = convert(T, BigInt(1643296191892)//BigInt(3432451463915))
  B̂7  = convert(T, BigInt(2576984903812)//BigInt(11692468803935))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4,B̂5,B̂6,B̂7)

  Bₗ   = convert(T, BigInt(0)//BigInt(1))
  B̂ₗ   = convert(T, BigInt(-2393889703871)//BigInt(16641202878460))

  C1  = convert(T2, BigInt(319960152914)//BigInt(39034091721739))
  C2  = convert(T2, BigInt(10916931475666701983218135)//BigInt(56630581182979020764713442))
  C3  = convert(T2, BigInt(31845189551971545944223680050155078355)//BigInt(113561670251926090809438891701398790454))
  C4  = convert(T2, BigInt(585892393366635581491792016142825500310911249371223)//BigInt(871432942801472160798333604371480303171919616321325))
  C5  = convert(T2, BigInt(6030664727234996630401450278844701818157369618311237)//BigInt(8305630304762506786823923305099106403075216590053000))
  C6  = convert(T2, BigInt(190737487565451971541550207118478711767748834018874068552898297)//BigInt(190737487565451971541550204260359567420033302718711745345318816))
  C7  = convert(T2, BigInt(194373043039840208108258122050794558876)//BigInt(388106905684556737922360607016380520227))
  Cᵢ   = SVector(C1,C2,C3,C4,C5,C6,C7)

  LowStorageRK4RPConstantCache{7,T,T2}(Aᵢ₁,Aᵢ₂,Aᵢ₃,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK85_4FM_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  uᵢ₋₃ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  fᵢ₋₃ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK85_4FM_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK4RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,fᵢ₋₂,fᵢ₋₃,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK85_4FM_4R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK85_4FM_4RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end


# 5R+ low storage methods introduced by van der Houwen
@cache struct LowStorageRK5RPCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  uᵢ₋₃::uType
  uᵢ₋₄::uType
  fᵢ₋₂::rateType
  fᵢ₋₃::rateType
  fᵢ₋₄::rateType
  gprev::uType
  fsalfirst::rateType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

struct LowStorageRK5RPConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  Aᵢ₁::SVector{N,T}
  Aᵢ₂::SVector{N,T}
  Aᵢ₃::SVector{N,T}
  Aᵢ₄::SVector{N,T}
  Bₗ::T
  B̂ₗ::T
  Bᵢ::SVector{N,T}
  B̂ᵢ::SVector{N,T}
  Cᵢ::SVector{N,T2}
end


function CKLLSRK75_4M_5RConstantCache(T, T2)
  A₁1  = convert(T, BigInt(984894634849)//BigInt(6216792334776))
  A₁2  = convert(T, BigInt(984894634849)//BigInt(5526037630912))
  A₁3  = convert(T, BigInt(13256335809797)//BigInt(10977774807827))
  A₁4  = convert(T, BigInt(5386479425293)//BigInt(11045691190948))
  A₁5  = convert(T, BigInt(-1717767168952)//BigInt(11602237717369))
  A₁6  = convert(T, BigInt(-10054679524430)//BigInt(10306851287569))
  Aᵢ₁  = SVector(A₁1,A₁2,A₁3,A₁4,A₁5,A₁6)

  A₂1  = convert(T, BigInt(0)//BigInt(1))
  A₂2  = convert(T, BigInt(890852251480)//BigInt(14995156510369))
  A₂3  = convert(T, BigInt(-18544705752398)//BigInt(18426539884027))
  A₂4  = convert(T, BigInt(1115398761892)//BigInt(28058504699217))
  A₂5  = convert(T, BigInt(5538441135605)//BigInt(13014942352969))
  A₂6  = convert(T, BigInt(23855853001162)//BigInt(20968156556405))
  Aᵢ₂  = SVector(A₂1,A₂2,A₂3,A₂4,A₂5,A₂6)

  A₃1  = convert(T, BigInt(0)//BigInt(1))
  A₃2  = convert(T, BigInt(0)//BigInt(1))
  A₃3  = convert(T, BigInt(1722683259617)//BigInt(5669183367476))
  A₃4  = convert(T, BigInt(342961171087)//BigInt(6505721096888))
  A₃5  = convert(T, BigInt(-14472869285404)//BigInt(19736045536601))
  A₃6  = convert(T, BigInt(-8169744035288)//BigInt(5424738459363))
  Aᵢ₃  = SVector(A₃1,A₃2,A₃3,A₃4,A₃5,A₃6)

  A₄1  = convert(T, BigInt(0)//BigInt(1))
  A₄2  = convert(T, BigInt(0)//BigInt(1))
  A₄3  = convert(T, BigInt(0)//BigInt(1))
  A₄4  = convert(T, BigInt(762111618422)//BigInt(5198184381557))
  A₄5  = convert(T, BigInt(2896263505307)//BigInt(6364015805096))
  A₄6  = convert(T, BigInt(60049403517654)//BigInt(26787923986853))
  Aᵢ₄  = SVector(A₄1,A₄2,A₄3,A₄4,A₄5,A₄6)

  B1  = convert(T, BigInt(1008141064049)//BigInt(9867084721348))
  B2  = convert(T, BigInt(0)//BigInt(1))
  B3  = convert(T, BigInt(8222186491841)//BigInt(18352662300888))
  B4  = convert(T, BigInt(514621697208)//BigInt(8712119383831))
  B5  = convert(T, BigInt(1808964136873)//BigInt(4546032443428))
  B6  = convert(T, BigInt(-362754645297)//BigInt(3989911846061))
  Bᵢ   = SVector(B1,B2,B3,B4,B5,B6)

  B̂1  = convert(T, BigInt(1633918545125)//BigInt(12016465907206))
  B̂2  = convert(T, BigInt(0)//BigInt(1))
  B̂3  = convert(T, BigInt(5614864639673)//BigInt(10804025076427))
  B̂4  = convert(T, BigInt(229286380958)//BigInt(6920724258831))
  B̂5  = convert(T, BigInt(5960415897193)//BigInt(14726168927560))
  B̂6  = convert(T, BigInt(-4042532386559)//BigInt(22820216867423))
  B̂ᵢ   = SVector(B̂1,B̂2,B̂3,B̂4,B̂5,B̂6)

  Bₗ   = convert(T, BigInt(599706619333)//BigInt(7161178965783))
  B̂ₗ   = convert(T, BigInt(930770261899)//BigInt(11134660916874))

  C1  = convert(T2, BigInt(984894634849)//BigInt(6216792334776))
  C2  = convert(T2, BigInt(19691532261044641782999041)//BigInt(82863799157714161922926528))
  C3  = convert(T2, BigInt(579140763944732527715749105230082493541)//BigInt(1146776047854201324825397010814855303604))
  C4  = convert(T2, BigInt(1904235205010770769196995566618512437342488019008993)//BigInt(2620260981179174237577004881164696841381017975634264))
  C5  = convert(T2, BigInt(4745866356039511505795256436748010529615723318082554645080208661)//BigInt(46784744516176933667763632070461960177241008032286254911869725672))
  C6  = convert(T2, BigInt(309879595293732553069368807532997606922999693101104106883289601491)//BigInt(309879595293732553069368804305686805880909932549908997963514738540))
  Cᵢ   = SVector(C1,C2,C3,C4,C5,C6)

  LowStorageRK5RPConstantCache{6,T,T2}(Aᵢ₁,Aᵢ₂,Aᵢ₃,Aᵢ₄,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ)
end

function alg_cache(alg::CKLLSRK75_4M_5R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  tmp  = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  k    = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  uᵢ₋₃ = zero(u)
  uᵢ₋₄ = zero(u)
  fᵢ₋₂ = zero(rate_prototype)
  fᵢ₋₃ = zero(rate_prototype)
  fᵢ₋₄ = zero(rate_prototype)
  gprev    = zero(u)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CKLLSRK75_4M_5RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  LowStorageRK5RPCache(u,uprev,k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,uᵢ₋₄,fᵢ₋₂,fᵢ₋₃,fᵢ₋₄,gprev,fsalfirst,tmp,atmp,tab)
end

function alg_cache(alg::CKLLSRK75_4M_5R,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  CKLLSRK75_4M_5RConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end
