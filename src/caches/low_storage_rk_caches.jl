
@cache struct CarpenterKennedy2N54Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct CarpenterKennedy2N54ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  A2end::SVector{4,T} # A1 is always zero
  B1::T
  B2end::SVector{4,T}
  c2end::SVector{4,T2} # c1 is always zero

  function CarpenterKennedy2N54ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    A2 = T(-567301805773/1357537059087)
    A3 = T(-2404267990393/2016746695238)
    A4 = T(-3550918686646/2091501179385)
    A5 = T(-1275806237668/842570457699)
    A2end = SVector{4,T}(A2, A3, A4, A5)

    B1 = T(1432997174477/9575080441755)
    B2 = T(5161836677717/13612068292357)
    B3 = T(1720146321549/2090206949498)
    B4 = T(3134564353537/4481467310338)
    B5 = T(2277821191437/14882151754819)
    B2end = SVector{4,T}(B2, B3, B4, B5)

    c2 = T2(1432997174477/9575080441755)
    c3 = T2(2526269341429/6820363962896)
    c4 = T2(2006345519317/3224310063776)
    c5 = T2(2802321613138/2924317926251)
    c2end = SVector{4,T2}(c2, c3, c4, c5)

    new{T,T2}(A2end, B1, B2end, c2end)
  end
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  CarpenterKennedy2N54Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end



@cache struct ParsaniKetchesonDeconinck3S184Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct ParsaniKetchesonDeconinck3S184ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  γ12end::SVector{17,T} # γ11 is always zero
  γ22end::SVector{17,T} # γ21 is always one
  γ32end::SVector{17,T} # γ31 is always zero
  δ2end ::SVector{17,T} # δ1  is always one
  β1::T
  β2end::SVector{17,T}
  c2end::SVector{17,T2} # c1 is always zero

  function ParsaniKetchesonDeconinck3S184ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    γ102 = T(1.1750819811951678e+0)
    γ103 = T(3.0909017892654811e-1)
    γ104 = T(1.4409117788115862e+0)
    γ105 = T(-4.3563049445694069e-1)
    γ106 = T(2.0341503014683893e-1)
    γ107 = T(4.9828356971917692e-1)
    γ108 = T(3.5307737157745489e+0)
    γ109 = T(-7.9318790975894626e-1)
    γ110 = T(8.9120513355345166e-1)
    γ111 = T(5.7091009196320974e-1)
    γ112 = T(1.6912188575015419e-2)
    γ113 = T(1.0077912519329719e+0)
    γ114 = T(-6.8532953752099512e-1)
    γ115 = T(1.0488165551884063e+0)
    γ116 = T(8.3647761371829943e-1)
    γ117 = T(1.3087909830445710e+0)
    γ118 = T(9.0419681700177323e-1)
    γ12end = SVector(γ102, γ103, γ104, γ105, γ106, γ107, γ108, γ109, γ110, γ111, γ112, γ113, γ114, γ115, γ116, γ117, γ118)

    γ202 = T(-1.2891068509748144e-1)
    γ203 = T(3.5609406666728954e-1)
    γ204 = T(-4.0648075226104241e-1)
    γ205 = T(6.0714786995207426e-1)
    γ206 = T(1.0253501186236846e+0)
    γ207 = T(2.4411240760769423e-1)
    γ208 = T(-1.2813606970134104e+0)
    γ209 = T(8.1625711892373898e-1)
    γ210 = T(1.0171269354643386e-1)
    γ211 = T(1.9379378662711269e-1)
    γ212 = T(7.4408643544851782e-1)
    γ213 = T(-1.2591764563430008e-1)
    γ214 = T(1.1996463179654226e+0)
    γ215 = T(4.5772068865370406e-2)
    γ216 = T(8.3622292077033844e-1)
    γ217 = T(-1.4179124272450148e+0)
    γ218 = T(1.3661459065331649e-1)
    γ22end = SVector(γ202, γ203, γ204, γ205, γ206, γ207, γ208, γ209, γ210, γ211, γ212, γ213, γ214, γ215, γ216, γ217, γ218)

    γ302 = T(0.0000000000000000e+0)
    γ303 = T(0.0000000000000000e+0)
    γ304 = T(2.5583378537249163e-1)
    γ305 = T(5.2676794366988289e-1)
    γ306 = T(-2.5648375621792202e-1)
    γ307 = T(3.1932438003236391e-1)
    γ308 = T(-3.1106815010852862e-1)
    γ309 = T(4.7631196164025996e-1)
    γ310 = T(-9.8853727938895783e-2)
    γ311 = T(1.9274726276883622e-1)
    γ312 = T(3.2389860855971508e-2)
    γ313 = T(7.5923980038397509e-2)
    γ314 = T(2.0635456088664017e-1)
    γ315 = T(-8.9741032556032857e-2)
    γ316 = T(2.6899932505676190e-2)
    γ317 = T(4.1882069379552307e-2)
    γ318 = T(6.2016148912381761e-2)
    γ32end = SVector(γ302, γ303, γ304, γ305, γ306, γ307, γ308, γ309, γ310, γ311, γ312, γ313, γ314, γ315, γ316, γ317, γ318)

    δ02 = T(3.5816500441970289e-1)
    δ03 = T(5.8208024465093577e-1)
    δ04 = T(-2.2615285894283538e-1)
    δ05 = T(-2.1715466578266213e-1)
    δ06 = T(-4.6990441450888265e-1)
    δ07 = T(-2.7986911594744995e-1)
    δ08 = T(9.8513926355272197e-1)
    δ09 = T(-1.1899324232814899e-1)
    δ10 = T(4.2821073124370562e-1)
    δ11 = T(-8.2196355299900403e-1)
    δ12 = T(5.8113997057675074e-2)
    δ13 = T(-6.1283024325436919e-1)
    δ14 = T(5.6800136190634054e-1)
    δ15 = T(-3.3874970570335106e-1)
    δ16 = T(-7.3071238125137772e-1)
    δ17 = T(8.3936016960374532e-2)
    δ18 = T(0.0000000000000000e+0)
    δ2end = SVector(δ02, δ03, δ04, δ05, δ06, δ07, δ08, δ09, δ10, δ11, δ12, δ13, δ14, δ15, δ16, δ17, δ18)

    β1 = T(1.2384169480626298e-1)
    β02 = T(1.0176262534280349e+0)
    β03 = T(-6.9732026387527429e-2)
    β04 = T(3.4239356067806476e-1)
    β05 = T(1.8177707207807942e-2)
    β06 = T(-6.1188746289480445e-3)
    β07 = T(7.8242308902580354e-2)
    β08 = T(-3.7642864750532951e-1)
    β09 = T(-4.5078383666690258e-2)
    β10 = T(-7.5734228201432585e-1)
    β11 = T(-2.7149222760935121e-1)
    β12 = T(1.1833684341657344e-3)
    β13 = T(2.8858319979308041e-2)
    β14 = T(4.6005267586974657e-1)
    β15 = T(1.8014887068775631e-2)
    β16 = T(-1.5508175395461857e-2)
    β17 = T(-4.0095737929274988e-1)
    β18 = T(1.4949678367038011e-1)
    β2end = SVector(β02, β03, β04, β05, β06, β07, β08, β09, β10, β11, β12, β13, β14, β15, β16, β17, β18)

    c02 = T2(1.2384169480626298e-1)
    c03 = T2(1.1574324659554065e+0)
    c04 = T2(5.4372099141546926e-1)
    c05 = T2(8.8394666834280744e-1)
    c06 = T2(-1.2212042176605774e-1)
    c07 = T2(4.4125685133082082e-1)
    c08 = T2(3.8039092095473748e-1)
    c09 = T2(5.4591107347528367e-2)
    c10 = T2(4.8731855535356028e-1)
    c11 = T2(-2.3007964303896034e-1)
    c12 = T2(-1.8907656662915873e-1)
    c13 = T2(8.1059805668623763e-1)
    c14 = T2(7.7080875997868803e-1)
    c15 = T2(1.1712158507200179e+0)
    c16 = T2(1.2755351018003545e+0)
    c17 = T2(8.0422507946168564e-1)
    c18 = T2(9.7508680250761848e-1)
    c2end = SVector(c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16, c17, c18)

    new{T,T2}(γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end)
  end
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S184,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = ParsaniKetchesonDeconinck3S184ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  ParsaniKetchesonDeconinck3S184Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::ParsaniKetchesonDeconinck3S184,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ParsaniKetchesonDeconinck3S184ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end
