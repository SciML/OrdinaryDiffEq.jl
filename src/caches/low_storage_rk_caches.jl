
# 2N low storage methods introduced by Williamson
@cache struct LowStorageRK2NCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct LowStorageRK2NConstantCache{N,T,T2} <: OrdinaryDiffEqConstantCache
  A2end::SVector{N,T} # A1 is always zero
  B1::T
  B2end::SVector{N,T}
  c2end::SVector{N,T2} # c1 is always zero
end


function ORK256ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ORK256,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ORK256ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ORK256,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ORK256ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function CarpenterKennedy2N54ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function HSLDDRK64ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::HSLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = HSLDDRK64ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::HSLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  HSLDDRK64ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function DGLDDRK73_CConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::DGLDDRK73_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = DGLDDRK73_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::DGLDDRK73_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  DGLDDRK73_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function DGLDDRK84_CConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::DGLDDRK84_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = DGLDDRK84_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::DGLDDRK84_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  DGLDDRK84_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function DGLDDRK84_FConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::DGLDDRK84_F,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = DGLDDRK84_FConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::DGLDDRK84_F,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  DGLDDRK84_FConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function NDBLSRK124ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::NDBLSRK124,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = NDBLSRK124ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::NDBLSRK124,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  NDBLSRK124ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function NDBLSRK134ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::NDBLSRK134,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = NDBLSRK134ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::NDBLSRK134,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  NDBLSRK134ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function NDBLSRK144ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::NDBLSRK144,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = NDBLSRK144ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2NCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::NDBLSRK144,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  NDBLSRK144ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
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



function CFRLDDRK64ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::CFRLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = CFRLDDRK64ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2CCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CFRLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  CFRLDDRK64ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


function TSLDDRK74ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::TSLDDRK74,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = TSLDDRK74ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  LowStorageRK2CCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::TSLDDRK74,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  TSLDDRK74ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
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


function ParsaniKetchesonDeconinck3S32ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S32ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S32ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S82ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S82,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S82ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S82,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S82ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S53ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S53ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S53ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S173ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S173,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S173ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S173,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S173ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S94ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S94,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S94ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S94,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S94ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S184ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S184,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S184ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S184,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S184ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S105ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S105,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S105ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S105,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S105ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end


function ParsaniKetchesonDeconinck3S205ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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

function alg_cache(alg::ParsaniKetchesonDeconinck3S205,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = ParsaniKetchesonDeconinck3S205ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  LowStorageRK3SCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S205,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S205ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end
