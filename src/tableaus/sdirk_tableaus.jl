struct TRBDF2Tableau{T,T2}
  γ::T2
  d::T
  ω::T
  btilde1::T
  btilde2::T
  btilde3::T
  α1::T
  α2::T
end

#=
Tableau:

c = [0; γ; 1]
A = [0  0  0
     d  d  0
     ω  ω  d]
b = [ω; ω; d]
bhat = [(1-ω)/3; (3ω+1)/3; d/3]

where γ=2-√2, d=γ/2, and ω=(√2)/4. Hence

btilde = bhat-b = [(1-4ω)/3; 1/3; -2d/3] = [(1-√2)/3; 1/3; (√2-2)/3]

Let uᵧ and u be approximations to u(t+γ*dt) and u(t+dt). Then initial guesses of
approximations zᵧ and z to dt*f(t+γ*dt,uᵧ) and dt*f(t+dt,u) according to Shampine:

zᵧ = zprev

z = (1.5+√2)*zprev + (2.5+2√2)*zᵧ - (6+4.5√2)*(uᵧ-uprev) = (by def. of uᵧ)
  = (1.5+√2)*zprev + (2.5+2√2)*zᵧ - (6+4.5√2)*(d*zᵧ + d*zprev) =
  = (-√2)/2*zprev + (1+(√2)/2)*zᵧ
=#
function TRBDF2Tableau(T,T2)
  γ = convert(T2,2-sqrt(2))
  d = convert(T,1-sqrt(2)/2)
  ω = convert(T,sqrt(2)/4)
  btilde1 = convert(T,(1-sqrt(2))/3)
  btilde2 = convert(T,1//3)
  btilde3 = convert(T,(sqrt(2)-2)/3)
  α1 = convert(T,-sqrt(2)/2)
  α2 = convert(T,1+sqrt(2)/2)
  TRBDF2Tableau(γ,d,ω,btilde1,btilde2,btilde3,α1,α2)
end

struct ESDIRK4Tableau{T,T2}
    γ::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    c3::T2
    α31::T
    α32::T
    α41::T
    α42::T
end

#=
Derivative of Hermite Polynomial
k[1] + Θ*(-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(3*dt*k[1] + 3*dt*k[2] + 6*y₀ - 6*y₁) + 6*y₁)/dt

Extrapolation for ESDIRK interior step 3
dt = c2 since interval is [c1,c2] and c1 = 0
θ =  c3/c2 the extrapolation point
z = dt*k

z₁ + Θ*(-4dt*z₁ - 2dt*z₂ - 6y₀ + Θ*(3dt*z₁ + 3z₂ + 6y₀ - 6y₁ ) + 6y₁)/dt

Test Expression on TRBDF2:
c2 = 2 - sqrt(2)
c3 = 1
θ = c3/c2; dt = c2

Coefficient on z₁:
(1 + (-4θ + 3θ^2))*z₁
1 + (-4θ + 3θ^2) - (1.5 + sqrt(2)) # 1.5 + sqrt(2) given by Shampine

Coefficient on z₂:
(-2θ + 3θ^2)*z₂
(-2θ + 3θ^2) - (2.5 + 2sqrt(2)) # 2.5 + 2sqrt(2) given by Shampine

Coefficient on y₁-y₀:
θ*(θ*6(y₀-y₁)+6(y₁-y₀))/dt
θ*(-6θ(y₁-y₀)+6(y₁-y₀))/dt
(y₁-y₀)6θ*(1-θ)/dt

(6θ*(1-θ)/dt)*(y₁-y₀)

6θ*(1-θ)/dt - (- (6 + 4.5sqrt(2)))  # - (6 + 4.5sqrt(2)) given by Shampine

# Write only in terms of z primatives
y₀ = uprev
y₁ = uprev + γ*z₁ + γ*z₂
y₁-y₀ = γ*z₁ + γ*z₂

# Full Expression
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/dt)*γ)*z₁ + ((-2θ + 3θ^2) + (6θ*(1-θ)/dt)*γ)*z₂
=#

#=
# Kvaerno3
# Predict z4 from yhat

yhat = uprev + a31*z1 + a32*z2 + γ*z3
z₄ = yhat - uprev = a31*z1 + a32*z2 + γ*z3

# Note Hermite is too small of an interval for this one!!!
=#

function Kvaerno3Tableau(T,T2)
  γ   = convert(T,0.4358665215)
  a31 = convert(T,0.490563388419108)
  a32 = convert(T,0.073570090080892)
  a41 = convert(T,0.308809969973036)
  a42 = convert(T,1.490563388254106)
  a43 = -convert(T,1.235239879727145)
  # bhat1 = convert(T,0.490563388419108)
  # bhat2 = convert(T,0.073570090080892)
  # bhat3 = convert(T,0.4358665215)
  # bhat4 = convert(T,0.0)
  btilde1 = convert(T,0.181753418446072) # bhat1-a41
  btilde2 = convert(T,-1.416993298173214) # bhat2-a42
  btilde3 = convert(T,1.671106401227145) # bhat3-a43
  btilde4 = -γ # bhat4-γ
  c3 = convert(T2,1)
  c2 = 2γ
  θ = c3/c2
  α31 = ((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
  α32 = ((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
  α41 = convert(T,0.0)
  α42 = convert(T,0.0)
  ESDIRK4Tableau(γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32,α41,α42)
end


struct KenCarp3Tableau{T,T2}
    γ::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    c3::T2
    α31::T
    α32::T
    α41::T
    α42::T
    ea21::T
    ea31::T
    ea32::T
    ea41::T
    ea42::T
    ea43::T
    eb1::T
    eb2::T
    eb3::T
    eb4::T
    ebtilde1::T
    ebtilde2::T
    ebtilde3::T
    ebtilde4 ::T
end

#=
# KenCarp3
# Predict z4 from Hermite z2 and z1
# Not z3 because c3 < c2 !

θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
θ = c4/c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
=#
function KenCarp3Tableau(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
  γ  = convert(T,0.435866521508459)
  a31 = convert(T,0.2576482460664272)
  a32 = -convert(T,0.09351476757488625)
  a41 = convert(T,0.18764102434672383)
  a42 = -convert(T,0.595297473576955)
  a43 = convert(T,0.9717899277217721)
  # bhat1 = convert(T,2756255671327//12835298489170)
  # bhat2 = -convert(T,10771552573575//22201958757719)
  # bhat3 = convert(T,9247589265047//10645013368117)
  # bhat4 = convert(T,2193209047091//5459859503100)
  btilde1 = convert(T,0.027099261876665316) # bhat1-a41
  btilde2 = convert(T,0.11013520969201586) # bhat2-a42
  btilde3 = convert(T,-0.10306492520138458) # bhat3-a43
  btilde4 = convert(T,-0.0341695463672966) # bhat4-γ
  c3 = convert(T2,0.6)
  c2 = 2γ
  θ = c3/c2
  α31 = ((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
  α32 = ((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
  θ = 1/c2
  α41 = ((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
  α42 = ((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

  ea21 = convert(T,0.871733043016918)
  ea31 = convert(T,0.5275890119763004)
  ea32 = convert(T,0.0724109880236996)
  ea41 = convert(T,0.3990960076760701)
  ea42 = -convert(T,0.4375576546135194)
  ea43 = convert(T,1.0384616469374492)
  eb1 = convert(T,0.18764102434672383)
  eb2 = convert(T,-0.595297473576955)
  eb3 = convert(T,0.9717899277217721)
  eb4 = convert(T,0.435866521508459)
  ebtilde1 = convert(T,0.027099261876665316)
  ebtilde2 = convert(T,0.11013520969201586)
  ebtilde3 = -convert(T,0.10306492520138458)
  ebtilde4 = -convert(T,0.0341695463672966)

  KenCarp3Tableau(γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,
                  α32,α41,α42,ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4,
                  ebtilde1,ebtilde2,ebtilde3,ebtilde4)
end



function KenCarp3Tableau(T,T2)
  γ  = convert(T,1767732205903//4055673282236)
  a31 = convert(T,2746238789719//10658868560708)
  a32 = -convert(T,640167445237//6845629431997)
  a41 = convert(T,1471266399579//7840856788654)
  a42 = -convert(T,4482444167858//7529755066697)
  a43 = convert(T,11266239266428//11593286722821)
  # bhat1 = convert(T,2756255671327//12835298489170)
  # bhat2 = -convert(T,10771552573575//22201958757719)
  # bhat3 = convert(T,9247589265047//10645013368117)
  # bhat4 = convert(T,2193209047091//5459859503100)
  btilde1 = convert(T,BigInt(681815649026867975666107)//BigInt(25159934323302256049469295)) # bhat1-a41
  btilde2 = convert(T,BigInt(18411887981491912264464127)//BigInt(167175311446532472108584143)) # bhat2-a42
  btilde3 = convert(T,BigInt(-12719313754959329011138489)//BigInt(123410692144842870217698057)) # bhat3-a43
  btilde4 = convert(T,BigInt(-47289384293135913063989)//BigInt(1383962894467812063558225)) # bhat4-γ
  c3 = convert(T2,3//5)
  c2 = 2γ
  θ = c3/c2
  α31 = ((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
  α32 = ((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
  θ = 1/c2
  α41 = ((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
  α42 = ((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

  # Explicit Tableau
  ea21 = convert(T,1767732205903//2027836641118)
  ea31 = convert(T,5535828885825//10492691773637)
  ea32 = convert(T,788022342437//10882634858940)
  ea41 = convert(T,6485989280629//16251701735622)
  ea42 = -convert(T,4246266847089//9704473918619)
  ea43 = convert(T,10755448449292//10357097424841)
  eb1 = convert(T,1471266399579//7840856788654)
  eb2 = convert(T,-4482444167858//7529755066697 )
  eb3 = convert(T,11266239266428//11593286722821)
  eb4 = convert(T,1767732205903//4055673282236)
  ebtilde1 = convert(T,BigInt(681815649026867975666107)//BigInt(25159934323302256049469295))
  ebtilde2 = convert(T,BigInt(18411887981491912264464127)//BigInt(167175311446532472108584143))
  ebtilde3 = -convert(T,BigInt(12719313754959329011138489)//BigInt(123410692144842870217698057))
  ebtilde4 = -convert(T,BigInt(47289384293135913063989)//BigInt(1383962894467812063558225))
  KenCarp3Tableau(γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,
                  α32,α41,α42,ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4,
                  ebtilde1,ebtilde2,ebtilde3,ebtilde4)
end

# Flip them all!

# ebtilde1 = big(1471266399579)//7840856788654 - big(2756255671327)//12835298489170
# ebtilde2 = -big(4482444167858)//7529755066697 + big(10771552573575)//22201958757719
# ebtilde3 = big(11266239266428)//11593286722821 - big(9247589265047)//10645013368117
# ebtilde4 = big(1767732205903)//4055673282236 - big(2193209047091)//5459859503100

struct Cash4Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  b1hat1::T
  b2hat1::T
  b3hat1::T
  b4hat1::T
  b1hat2::T
  b2hat2::T
  b3hat2::T
  b4hat2::T
  c2::T2
  c3::T2
  c4::T2
end

#=
Extrapolation for Cash interior step 3
dt = c1-c2 since interval is [c2,c1] and c1 = 0
c2 < c1, so z₂ is the left
θ =  (c3-c1)/dt the extrapolation point
z = dt*k

z₂ + Θ*(-4dt*z₂ - 2dt*z₁ - 6y₀ + Θ*(3dt*z₂ + 3z₁ + 6y₀ - 6y₁ ) + 6y₁)/dt

Coefficient on z₁:
(-2θ + 3θ^2)

Coefficient on z₂:
(1 + (-4θ + 3θ^2))

Coefficient on y₁-y₀:
(6θ*(1-θ)/dt)

# Write only in terms of z primatives
y₁ = uprev + a21*z₁ + γ*z₂
y₀ = uprev + γ*z₁
y₁-y₀ = (a21-γ)*z₁ + γ*z₂

θ = 1.1
Full z₁ coefficient: (-2θ + 3θ^2) + (6θ*(1-θ)/dt)*(a21-γ)
Full z₂ coefficient: (1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/dt)*γ

f(θ)= (-2θ + 3θ^2) + (6θ*(1-θ)/dt)*(a21-γ)
g(θ) = (1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/dt)*γ
t = linspace(0,1.5,100)
y = f.(t)
z = g.(t)
plot(t,y)
plot!(t,z)

The extrapolation is really bad that far
Hairer's extrapolation is no better.
Using constant extrapolations
=#

function Cash4Tableau(T,T2)
  γ = convert(T,0.435866521508)
  a21 = convert(T,-1.13586652150)
  a31 = convert(T,1.08543330679)
  a32 = -convert(T,0.721299828287)
  a41 = convert(T,0.416349501547)
  a42 = convert(T,0.190984004184)
  a43 = convert(T,-0.118643265417)
  a51 = convert(T,0.896869652944)
  a52 = convert(T,0.0182725272734)
  a53 = -convert(T,0.0845900310706)
  a54 = -convert(T,0.266418670647)
  b1hat1 = convert(T,1.05646216107052)
  b2hat1 = -convert(T,0.0564621610705236)
  b3hat1 = convert(T,0)
  b4hat1 = convert(T,0)
  b1hat2 = convert(T,0.776691932910)
  b2hat2 = convert(T,0.0297472791484)
  b3hat2 = -convert(T,0.0267440239074)
  b4hat2 = convert(T,0.220304811849)
  c2 = -convert(T2,0.7)
  c3 = convert(T2,0.8)
  c4 = convert(T2,0.924556761814)
  Cash4Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,
               b1hat1,b2hat1,b3hat1,b4hat1,b1hat2,b2hat2,b3hat2,b4hat2,c2,c3,c4)
end

struct SFSDIRK4Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  c2::T2
  c3::T2
  c4::T2
end

function SFSDIRK4Tableau(T,T2)
  γ = convert(T,0.097961082941)
  a21 = convert(T,0.262318069183)
  a31 = convert(T,0.230169419019)
  a32 = convert(T,0.294466719347)
  a41 = convert(T,0.210562684389)
  a42 = convert(T,0.269382888280)
  a43 = convert(T,0.307008634881)
  a51 = convert(T,0.222119403264)
  a52 = convert(T,0.282060762166)
  a53 = convert(T,0.236881213175)
  a54 = convert(T,0.258938621395)
  c2 = convert(T2,0.360279152124)
  c3 = convert(T2,0.622597221298)
  c4 = convert(T2,0.884915290491)
  SFSDIRK4Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4)
end

struct SFSDIRK5Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  a61::T
  a62::T
  a63::T
  a64::T
  a65::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
end

function SFSDIRK5Tableau(T,T2)
  γ = convert(T,0.078752939968)
  a21 = convert(T,0.222465723027)
  a31 = convert(T,0.203192361700)
  a32 = convert(T,0.230847263068)
  a41 = convert(T,0.188022704389)
  a42 = convert(T,0.191735630027)
  a43 = convert(T,0.209922288451)
  a51 = convert(T,0.188025114093)
  a52 = convert(T,0.191739898281)
  a53 = convert(T,0.209907601860)
  a54 = convert(T,0.252726086329)
  a61 = convert(T,0.192143833571)
  a62 = convert(T,0.200935182974)
  a63 = convert(T,0.205799262036)
  a64 = convert(T,0.200553844640)
  a65 = convert(T,0.200567876778)
  c2 = convert(T2,0.301218662995)
  c3 = convert(T2,0.512792564736)
  c4 = convert(T2,0.668433562835)
  c5 = convert(T2,0.921151640531)
  SFSDIRK5Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,c2,c3,c4,c5)
end

struct SFSDIRK6Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  a61::T
  a62::T
  a63::T
  a64::T
  a65::T
  a71::T
  a72::T
  a73::T
  a74::T
  a75::T
  a76::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
end

function SFSDIRK6Tableau(T,T2)
  γ = convert(T,0.067410767219)
  a21 = convert(T,0.194216850802)
  a31 = convert(T,0.194216850802)
  a32 = convert(T,0.199861501713)
  a41 = convert(T,0.162188551749)
  a42 = convert(T,0.166902343330)
  a43 = convert(T,0.145120313717)
  a51 = convert(T,0.165176818500)
  a52 = convert(T,0.169977460026)
  a53 = convert(T,0.150227711763)
  a54 = convert(T,0.181214258555)
  a61 = convert(T,0.165176818500)
  a62 = convert(T,0.169977460026)
  a63 = convert(T,0.150227711763)
  a64 = convert(T,0.181214258555)
  a65 = convert(T,0.199861501713)
  a71 = convert(T,0.168954170460)
  a72 = convert(T,0.173864595628)
  a73 = convert(T,0.156683775305)
  a74 = convert(T,0.157643002581)
  a75 = convert(T,0.173864725004)
  a76 = convert(T,0.168989731022)
  c2 = convert(T2,0.261627618021)
  c3 = convert(T2,0.461489119734)
  c4 = convert(T2,0.541621976015)
  c5 = convert(T2,0.734007016063)
  c6 = convert(T2,0.933868517776)
  SFSDIRK6Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,c2,c3,c4,c5,c6)
end

struct SFSDIRK7Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  a61::T
  a62::T
  a63::T
  a64::T
  a65::T
  a71::T
  a72::T
  a73::T
  a74::T
  a75::T
  a76::T
  a81::T
  a82::T
  a83::T
  a84::T
  a85::T
  a86::T
  a87::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
end

function SFSDIRK7Tableau(T,T2)
  γ = convert(T,0.056879041592)
  a21 = convert(T,0.172205581756)
  a31 = convert(T,0.135485903539)
  a32 = convert(T,0.135485903539)
  a41 = convert(T,0.133962606568)
  a42 = convert(T,0.133962606568)
  a43 = convert(T,0.170269437596)
  a51 = convert(T,0.133962606568)
  a52 = convert(T,0.133962606568)
  a53 = convert(T,0.170269437596)
  a54 = convert(T,0.172205581756)
  a61 = convert(T,0.138004377067)
  a62 = convert(T,0.133084723451)
  a63 = convert(T,0.152274237527)
  a64 = convert(T,0.154005757170)
  a65 = convert(T,0.154005757170)
  a71 = convert(T,0.139433665640)
  a72 = convert(T,0.134719607258)
  a73 = convert(T,0.145910607076)
  a74 = convert(T,0.147569765489)
  a75 = convert(T,0.147569765489)
  a76 = convert(T,0.165009008641)
  a81 = convert(T,0.138370770799)
  a82 = convert(T,0.134572540279)
  a83 = convert(T,0.150642940425)
  a84 = convert(T,0.152355910489)
  a85 = convert(T,0.152355910489)
  a86 = convert(T,0.132951737506)
  a87 = convert(T,0.138750190012)
  c2 = convert(T2,0.229084623348)
  c3 = convert(T2,0.32785084867)
  c4 = convert(T2,0.495073692324)
  c5 = convert(T2,0.66727927408)
  c6 = convert(T2,0.788253893977)
  c7 = convert(T2,0.937091461185)
  SFSDIRK7Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a82,a83,a84,a85,a86,a87,c2,c3,c4,c5,c6,c7)
end

struct SFSDIRK8Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  a61::T
  a62::T
  a63::T
  a64::T
  a65::T
  a71::T
  a72::T
  a73::T
  a74::T
  a75::T
  a76::T
  a81::T
  a82::T
  a83::T
  a84::T
  a85::T
  a86::T
  a87::T
  a91::T
  a92::T
  a93::T
  a94::T
  a95::T
  a96::T
  a97::T
  a98::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2
end

function SFSDIRK8Tableau(T,T2)
  γ = convert(T,0.050353353407)
  a21 = convert(T,0.147724666662)
  a31 = convert(T,0.114455029802)
  a32 = convert(T,0.114455029802)
  a41 = convert(T,0.114147680771)
  a42 = convert(T,0.114147680771)
  a43 = convert(T,0.147327977820)
  a51 = convert(T,0.114163314686)
  a52 = convert(T,0.114163314686)
  a53 = convert(T,0.147259379853)
  a54 = convert(T,0.147655883990)
  a61 = convert(T,0.114163314686)
  a62 = convert(T,0.114163314686)
  a63 = convert(T,0.147259379853)
  a64 = convert(T,0.147655883990)
  a65 = convert(T,0.147724666662)
  a71 = convert(T,0.118472990244)
  a72 = convert(T,0.118472990244)
  a73 = convert(T,0.128349529304)
  a74 = convert(T,0.128695117609)
  a75 = convert(T,0.128755067770)
  a76 = convert(T,0.128755067770)
  a81 = convert(T,0.118472990244)
  a82 = convert(T,0.118472990244)
  a83 = convert(T,0.128349529304)
  a84 = convert(T,0.128695117609)
  a85 = convert(T,0.128755067770)
  a86 = convert(T,0.128755067770)
  a87 = convert(T,0.147724666662)
  a91 = convert(T,0.117592883046)
  a92 = convert(T,0.117592883046)
  a93 = convert(T,0.132211234288)
  a94 = convert(T,0.132567220450)
  a95 = convert(T,0.132628974356)
  a96 = convert(T,0.132293123539)
  a97 = convert(T,0.117556840638)
  a98 = convert(T,0.117556840638)
  c2 = convert(T2,0.198078020069)
  c3 = convert(T2,0.279263413011)
  c4 = convert(T2,0.425976692769)
  c5 = convert(T2,0.573595246622)
  c6 = convert(T2,0.721319913284)
  c7 = convert(T2,0.801854116348)
  c8 = convert(T2,0.94957878301)
  SFSDIRK8Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a82,a83,a84,a85,a86,a87,a91,a92,a93,a94,a95,a96,a97,a98,c2,c3,c4,c5,c6,c7,c8)
end

struct Hairer4Tableau{T,T2}
  γ::T
  a21::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  bhat1::T
  bhat2::T
  bhat3::T
  bhat4::T
  btilde1::T
  btilde2::T
  btilde3::T
  btilde4::T
  btilde5::T
  c2::T2
  c3::T2
  c4::T2
  r11::T
  r12::T
  r13::T
  r14::T
  r21::T
  r22::T
  r23::T
  r24::T
  r31::T
  r32::T
  r33::T
  r34::T
  r41::T
  r42::T
  r43::T
  r44::T
  r51::T
  r52::T
  r53::T
  r54::T
  α21::T
  α31::T
  α32::T
  α41::T
  α43::T
end

function Hairer4Tableau(T,T2)
  γ = convert(T,1//4)
  c2 = convert(T,3//4)
  c3 = convert(T,11//20)
  c4 = convert(T,1//2)

  #=
  α21 = convert(T,2)
  α31 = convert(T,42//25)
  α32 = -convert(T,4//25)
  α41 = convert(T,89//68)
  α42 = -convert(T,25//136)
  α43 = convert(T,15//136)
  α51 = -convert(T,37//12)
  α52 = -convert(T,103//24)
  α53 = convert(T,275//8)
  α54 = -convert(T,85//3)

  alpha = -inv(A)*γ

  α = [-1 0 0 0 0
       α21 -1 0 0 0
       α31 α32 -1 0 0
       α41 α42 α43 -1 0
       α51 α52 α53 α54 -1]
  A α = -γ
  A = -γ*inv(α)

  Now zⱼ = fⱼ + ∑α_ⱼᵢ zᵢ
  =#

  a11 = convert(T,1//4)
  a21 = convert(T,1//2)
  a22 = convert(T,1//4)
  a31 = convert(T,17//50)
  a32 = convert(T,-1//25)
  a33 = convert(T,1//4)
  a41 = convert(T,371//1360)
  a42 = convert(T,-137//2720)
  a43 = convert(T,15//544)
  a44 = convert(T,1//4)
  a51 = convert(T,25//24)
  a52 = convert(T,-49//48)
  a53 = convert(T,125//16)
  a54 = convert(T,-85//12)

  #=
  e1 = -convert(T,23//6)
  e2 = -convert(T,17//12)
  e3 = convert(T,125//4)
  e4 = -convert(T,85//3)
  E = [e1 e2 e3 e4 0]

  bhat = [59//48,-17//96,225//32,-85//12,0]

  α = [-1 0 0 0 0
       α21 -1 0 0 0
       α31 α32 -1 0 0
       α41 α42 α43 -1 0
       e1 e2 e3 e4 -1]

  A = [γ 0 0 0 0
       a21 γ 0 0 0
       a31 a32 γ 0 0
       a41 a42 a43 γ 0
       a51 a52 a53 a54 γ]

  E = bhat'*inv(A)
  bhat = E*A
  =#

  bhat1 = convert(T,59//48)
  bhat2 = convert(T,-17//96)
  bhat3 = convert(T,225//32)
  bhat4 = convert(T,-85//12)

  btilde1 = convert(T,3//16) # bhat1-a51
  btilde2 = convert(T,27//32) # bhat2-a52
  btilde3 = convert(T,-25//32) # bhat3-a53
  btilde4 = convert(T,0) # bhat4-a54
  btilde5 = -γ

  #=
  d11 = convert(T,61//27)
  d12 = convert(T,-185//54)
  d13 = convert(T,2525//18)
  d14 = convert(T,-3740//27)
  d15 = convert(T,-44//9)
  d21 = convert(T,2315//81)
  d22 = convert(T,1049//162)
  d23 = convert(T,-27725//54)
  d24 = convert(T,40460//81)
  d25 = convert(T,557//27)
  d31 = convert(T,-6178//81)
  d32 = convert(T,-1607//81)
  d33 = convert(T,20075//27)
  d34 = convert(T,-56440//81)
  d35 = convert(T,-718//27)
  d41 = convert(T,3680//81)
  d42 = convert(T,1360//81)
  d43 = convert(T,-10000//27)
  d44 = convert(T,27200//81)
  d45 = convert(T,320//27)

  D = [d11 d12 d13 d14 d15
       d21 d22 d23 d24 d25
       d31 d32 d33 d34 d35
       d41 d42 d43 d44 d45]
  R = (D*A)'
  =#

  r11 = convert(T,11//3)
  r12 = convert(T,-463//72)
  r13 = convert(T,217//36)
  r14 = convert(T,-20//9)
  r21 = convert(T,11//2)
  r22 = convert(T,-385//16)
  r23 = convert(T,661//24)
  r24 = convert(T,-10//1)
  r31 = convert(T,-125//18)
  r32 = convert(T,20125//432)
  r33 = convert(T,-8875//216)
  r34 = convert(T,250//27)
  r41 = convert(T,0)
  r42 = convert(T,-85//4)
  r43 = convert(T,85//6)
  r44 = convert(T,0//1)
  r51 = convert(T,-11//9)
  r52 = convert(T,557//108)
  r53 = convert(T,-359//54)
  r54 = convert(T,80//27)

  # c2/γ
  α21 = convert(T,3)
  #=
  # Prediction alphas from Hairer
  # Predict z3 from z1 and z2
  A = [c1 c2
  γ*c1 a21*c1+γ*c2]
  b = [c3,a31*c1+a32*c2+γ*c3]
  A\b
  =#
  α31 = convert(T,88//100)
  α32 = convert(T,44//100)
  #=
  # Predict z4 from z1 and z3
  A = [c1   c3
       γ*c1 a31*c1+a32*c2+γ*c3]
  b = [c4,a41*c1+a42*c2+a43*c3+γ*c4]
  A\b
  =#
  α41 = convert(T,3//17)
  α43 = convert(T,155//187)

 Hairer4Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,
                bhat1,bhat2,bhat3,bhat4,btilde1,btilde2,btilde3,
                btilde4,btilde5,c2,c3,c4,r11,r12,r13,r14,
                r21,r22,r23,r24,r31,r32,r33,r34,r41,r42,r43,r44,r51,
                r52,r53,r54,α21,α31,α32,α41,α43)
end



function Hairer42Tableau(T,T2)
  γ  = convert(T,4//15)
  c2 = convert(T2,23//30)
  c3 = convert(T2,17//30)
  c4 = convert(T2,2881//28965)+γ
  #=
  α21 = convert(T,15//8)
  α31 = convert(T,1577061//922880)
  α32 = -convert(T,23427//115360)
  α41 = convert(T,647163682356923881//2414496535205978880)
  α42 = -convert(T,593512117011179//3245291041943520)
  α43 = convert(T,559907973726451//1886325418129671)
  α51 = convert(T,724545451//796538880)
  α52 = -convert(T,830832077//267298560)
  α53 = convert(T,30957577//2509272)
  α54 = -convert(T,69863904375173//6212571137048)

  α = [-1 0 0 0 0
       α21 -1 0 0 0
       α31 α32 -1 0 0
       α41 α42 α43 -1 0
       α51 α52 α53 α54 -1]
  A = -γ*inv(α)
  =#
  a21 = convert(T,1//2)
  a31 = convert(T,51069//144200)
  a32 = convert(T,-7809//144200)
  a41 = convert(T,12047244770625658//141474406359725325)
  a42 = convert(T,-3057890203562191//47158135453241775)
  a43 = convert(T,2239631894905804//28294881271945065)
  a51 = convert(T,181513//86430)
  a52 = convert(T,-89074//116015)
  a53 = convert(T,83636//34851)
  a54 = convert(T,-69863904375173//23297141763930)

  #=

  A = [γ 0 0 0 0
       a21 γ 0 0 0
       a31 a32 γ 0 0
       a41 a42 a43 γ 0
       a51 a52 a53 a54 γ]

  A = convert(Array{Rational{BigInt}},A)

  e1 = convert(T,7752107607//11393456128)
  e2 = -convert(T,17881415427//11470078208)
  e3 = convert(T,2433277665//179459416)
  e4 = -convert(T,96203066666797//6212571137048)
  E = [e1 e2 e3 e4 0]

  bhat = E*A
  =#

  bhat1 = convert(T,33665407//11668050)
  bhat2 = convert(T,-2284766//15662025)
  bhat3 = convert(T,11244716//4704885)
  bhat4 = convert(T,-96203066666797//23297141763930)

  btilde1 = convert(T,4580576//5834025) # bhat1-a51
  btilde2 = convert(T,9740224//15662025) # bhat2-a52
  btilde3 = convert(T,-46144//4704885) # bhat3-a53
  btilde4 = convert(T,-13169581145812//11648570881965) # bhat4-a54
  btilde5 = -γ

  #=
  d11 = convert(T,24.74416644927758)
  d12 = -convert(T,4.325375951824688)
  d13 = convert(T,41.39683763286316)
  d14 = convert(T,-61.04144619901784)
  d15 = convert(T,-3.391332232917013)
  d21 = convert(T,-51.98245719616925)
  d22 = convert(T,10.52501981094525)
  d23 = convert(T,-154.2067922191855)
  d24 = convert(T,214.3082125319825)
  d25 = convert(T,14.71166018088679)
  d31 = convert(T,33.14347947522142)
  d32 = convert(T,-19.72986789558523)
  d33 = convert(T,230.4878502285804)
  d34 = convert(T,-287.6629744338197)
  d35 = convert(T,-18.99932366302254)
  d41 = convert(T,-5.905188728329743)
  d42 = convert(T,13.53022403646467)
  d43 = convert(T,-117.6778956422581)
  d44 = convert(T,134.3962081008550)
  d45 = convert(T,8.678995715052762)
  =#

  r11 = convert(T,6.776439256624082)
  r12 = convert(T,-14.066831911883533)
  r13 = convert(T,16.204808856162565)
  r14 = convert(T,-6.8143005003361985)
  r21 = convert(T,3.166691382949011)
  r22 = convert(T,-14.034196189427504)
  r23 = convert(T,15.497198116229603)
  r24 = convert(T,-5.3974733381957005)
  r31 = convert(T,-1.9310469085972866)
  r32 = convert(T,11.146663701107887)
  r33 = convert(T,-6.9009212321038405)
  r34 = convert(T,0.085120800673252)
  r41 = convert(T,-6.107728468864632)
  r42 = convert(T,13.031255018633459)
  r43 = convert(T,-19.734599430149146)
  r44 = convert(T,9.812254180511282)
  r51 = convert(T,-0.9043552621112034)
  r52 = convert(T,3.9231093815698106)
  r53 = convert(T,-5.066486310139344)
  r54 = convert(T,2.3143988573474035)

  # c2/γ
  α21 = convert(T,23//8)
  α31 = convert(T,0.9838473040915402)
  α32 = convert(T,0.3969226768377252)
  α41 = convert(T,0.6563374010466914)
  α43 = convert(T,0.3372498196189311)

  Hairer4Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,
                 bhat1,bhat2,bhat3,bhat4,btilde1,btilde2,btilde3,btilde4,btilde5,
                 c2,c3,c4,r11,r12,r13,r14,
                 r21,r22,r23,r24,r31,r32,r33,r34,r41,r42,r43,r44,r51,
                 r52,r53,r54,α21,α31,α32,α41,α43)
end

struct Kvaerno4Tableau{T,T2}
  γ::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  btilde1::T
  btilde2::T
  btilde3::T
  btilde4::T
  btilde5::T
  c3::T2
  c4::T2
  α21::T
  α31::T
  α32::T
  α41::T
  α42::T
end

#=
# Kvaerno4
# Predict z3 from Hermite z2 and z1

c2 = 2γ
θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z4 from Hermite z2 and z1

θ = c4/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
=#

function Kvaerno4Tableau(T,T2)
  γ = convert(T,0.4358665215)
  a31 = convert(T,0.140737774731968)
  a32 = convert(T,-0.108365551378832)
  a41 = convert(T,0.102399400616089)
  a42 = convert(T,-0.376878452267324)
  a43 = convert(T,0.838612530151233)
  a51 = convert(T,0.157024897860995)
  a52 = convert(T,0.117330441357768)
  a53 = convert(T,0.61667803039168)
  a54 = convert(T,-0.326899891110444)
  btilde1 = convert(T,-0.054625497244906) # a41 - a51
  btilde2 = convert(T,-0.494208893625092) # a42 - a52
  btilde3 = convert(T,0.221934499759553) # a43 - a53
  btilde4 = convert(T,0.762766412610444) # γ - a54
  btilde5 = -γ
  c3 = convert(T2,0.468238744853136)
  c4 = convert(T2,1)
  α21 = convert(T,2) # c2/γ
  α31 = convert(T,0.462864521870446)
  α32 = convert(T,0.537135478129554)
  α41 = convert(T,-0.14714018016178376)
  α42 = convert(T,1.1471401801617838)
  Kvaerno4Tableau(γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,
                  btilde1,btilde2,btilde3,btilde4,btilde5,
                  c3,c4,α21,α31,α32,α41,α42)
end

struct KenCarp4Tableau{T,T2}
  γ::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  a61::T
  a63::T
  a64::T
  a65::T
  btilde1::T
  btilde3::T
  btilde4::T
  btilde5::T
  btilde6::T
  c3::T2
  c4::T2
  c5::T2
  α21::T
  α31::T
  α32::T
  α41::T
  α42::T
  α51::T
  α52::T
  α53::T
  α54::T
  α61::T
  α62::T
  α63::T
  α64::T
  α65::T
  ea21::T
  ea31::T
  ea32::T
  ea41::T
  ea42::T
  ea43::T
  ea51::T
  ea52::T
  ea53::T
  ea54::T
  ea61::T
  ea62::T
  ea63::T
  ea64::T
  ea65::T
  eb1::T
  eb3::T
  eb4::T
  eb5::T
  eb6::T
  ebtilde1::T
  ebtilde3::T
  ebtilde4::T
  ebtilde5::T
  ebtilde6::T
end

struct CFNLIRK3Tableau{T,T2}
    γ::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    c2::T2
    c3::T2
    ea21::T
    ea31::T
    ea32::T
    ea41::T
    ea42::T
    ea43::T
    eb1::T
    eb2::T
    eb3::T
    eb4::T
end

function CFNLIRK3Tableau(T,T2)
  #Implicit Tableau
  γ  = convert(T,0.43586652150846)
  a31 = convert(T,0.0)
  a32 = convert(T,(1-γ)/2)
  a41 = convert(T,0.0)
  a42 = convert(T,1.20849664917601276)
  a43 = convert(T,-0.64436317068447276)
  c3 = (1+γ)/2
  c2 = γ

  # Explicit Tableau
  ea21 = convert(T,γ)
  ea31 = convert(T,(1.7+γ)/2)
  ea32 = convert(T,-0.35)
  ea41 = convert(T,0.0)
  ea42 = convert(T,1.9891757246798590)
  ea43 = convert(T,-0.9891757246798590)
  eb1 = convert(T,0.0)
  eb2 = convert(T,1.20849664917601276)
  eb3 = convert(T,-0.64436317068447276)
  eb4 = convert(T,γ)
  CFNLIRK3Tableau(γ,a31,a32,a41,a42,a43,c2,c3,ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4)
end

#=
# KenCarp4
# Predict z3 from Hermite z2 and z1

c2 = 2γ
θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z4 from Hermite z2 and z1
θ = c4/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z5 from Hermite z2 and z1
θ = c5/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z5 from z1 and z4

θ = c5/c4
dt = c4

y₀ = uprev
y₁ = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄
y₁-y₀ = a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

(1 + (-4θ + 3θ^2) + a31*(6θ*(1-θ)/dt))*z₁ +
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))*z₅
+ (6θ*(1-θ)/dt)*(a52*z₂ + a53*z₃ + a54*z₄)

(1 + (-4θ + 3θ^2) + a41*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a42
(6θ*(1-θ)/dt)*a43
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))

# Predict last stage from z1 and z5

θ = 1/c5
dt = c5

y₀ = uprev
y₁ = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅
y₁-y₀ = a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅

(1 + (-4θ + 3θ^2) + a31*(6θ*(1-θ)/dt))*z₁ +
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))*z₅
+ (6θ*(1-θ)/dt)*(a52*z₂ + a53*z₃ + a54*z₄)

(1 + (-4θ + 3θ^2) + a51*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a52
(6θ*(1-θ)/dt)*a53
(6θ*(1-θ)/dt)*a54
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))

=#

struct CFNLIRK4Tableau{T,T2}
    γ::T
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    a51::T
    a52::T
    a53::T
    a54::T
    a61::T
    a62::T
    a63::T
    a64::T
    a65::T
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    ea21::T
    ea31::T
    ea32::T
    ea41::T
    ea42::T
    ea43::T
    ea51::T
    ea52::T
    ea53::T
    ea54::T
    ea61::T
    ea62::T
    ea63::T
    ea64::T
    ea65::T
    eb1::T
    eb2::T
    eb3::T
    eb4::T
    eb5::T
    eb6::T
end

function CFNLIRK4Tableau(T,T2)
  #Implicit Tableau
  γ  = convert(T,1/4)
  a21  = convert(T,0)
  a31  = convert(T,0)
  a32  = convert(T,1/2)
  a41  = convert(T,0)
  a42  = convert(T,17/50)
  a43  = convert(T,-1/25)
  a51  = convert(T,0)
  a52  = convert(T,371/1360)
  a53  = convert(T,-137/2720)
  a54  = convert(T,15/544)
  a61  = convert(T,0.0)
  a62  = convert(T,25/24)
  a63  = convert(T,-49/48)
  a64  = convert(T,125/16)
  a65  = convert(T,-85/12)
  c2  = convert(T2,1/4)
  c3  = convert(T2,3/4)
  c4  = convert(T2,11/20)
  c5  = convert(T2,1/2)
  ea21  = convert(T,1/4)
  ea31  = convert(T,-1/4)
  ea32  = convert(T,1)
  ea41  = convert(T,-13/100)
  ea42  = convert(T,43/75)
  ea43  = convert(T,8/75)
  ea51  = convert(T,-6/85)
  ea52  = convert(T,42/85)
  ea53  = convert(T,179/1360)
  ea54  = convert(T,-15/272)
  ea61  = convert(T,0.0)
  ea62  = convert(T,79/24)
  ea63  = convert(T,-5/8)
  ea64  = convert(T,25/2)
  ea65  = convert(T,-85/6)
  eb1  = convert(T,0.0)
  eb2  = convert(T,25/24)
  eb3  = convert(T,-49/48)
  eb4  = convert(T,125/16)
  eb5  = convert(T,-85/12)
  eb6  = convert(T,1/4)
  CFNLIRK4Tableau(γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,c2,c3,c4,c5,ea21,ea31,ea32,ea41,ea42,ea43,ea51,ea52,ea53,ea54,ea61,ea62,ea63,ea64,ea65,eb1,eb2,eb3,eb4,eb5,eb6)
end

function KenCarp4Tableau(T,T2)
  γ = convert(T,1//4)
  a31 = convert(T,8611//62500)
  a32 = -convert(T,1743//31250)
  a41 = convert(T,5012029//34652500)
  a42 = -convert(T,654441//2922500)
  a43 = convert(T,174375//388108)
  a51 = convert(T,15267082809//155376265600)
  a52 = -convert(T,71443401//120774400)
  a53 = convert(T,730878875//902184768)
  a54 = convert(T,2285395//8070912)
  a61 = convert(T,82889//524892)
  a63 = convert(T,15625//83664)
  a64 = convert(T,69875//102672)
  a65 = -convert(T,2260//8211)
  # bhat1 = convert(T,4586570599//29645900160)
  # bhat3 = convert(T,178811875//945068544)
  # bhat4 = convert(T,814220225//1159782912)
  # bhat5 = -convert(T,3700637//11593932)
  # bhat6 = convert(T,61727//225920)
  btilde1 = convert(T,-31666707//9881966720) # bhat1-a61
  btilde3 = convert(T,256875//105007616) # bhat3-a63
  btilde4 = convert(T,2768025//128864768) # bhat4-a64
  btilde5 = convert(T,-169839//3864644) # bhat5-a65
  btilde6 = convert(T,5247//225920) # bhat6-γ
  c3 = convert(T2,83//250)
  c4 = convert(T2,31//50)
  c5 = convert(T2,17//20)
  α21 = convert(T,2) # c2/γ
  α31 = convert(T,42//125)
  α32 = convert(T,83//125)
  α41 = convert(T,-6//25)
  α42 = convert(T,31//25)
  α51 = convert(T,914470432//2064665255)
  α52 = convert(T,798813//724780)
  α53 = convert(T,-824765625//372971788)
  α54 = convert(T,49640//29791)
  α61 = convert(T,288521442795//954204491116)
  α62 = convert(T,2224881//2566456)
  α63 = convert(T,-1074821875//905317354)
  α64 = convert(T,-3360875//8098936)
  α65 = convert(T,7040//4913)

  ea21 = convert(T,1//2)
  ea31 = convert(T,13861//62500)
  ea32 = convert(T,6889//62500)
  ea41 = -convert(T,116923316275//2393684061468)
  ea42 = -convert(T,2731218467317//15368042101831)
  ea43 = convert(T,9408046702089//11113171139209)
  ea51 = -convert(T,451086348788//2902428689909)
  ea52 = -convert(T,2682348792572//7519795681897)
  ea53 = convert(T,12662868775082//11960479115383)
  ea54 = convert(T,3355817975965//11060851509271)
  ea61 = convert(T,647845179188//3216320057751)
  ea62 = convert(T,73281519250//8382639484533)
  ea63 = convert(T,552539513391//3454668386233)
  ea64 = convert(T,3354512671639//8306763924573)
  ea65 = convert(T,4040//17871)

  eb1 = convert(T,82889//524892)
  eb3 = convert(T,15625//83664)
  eb4 = convert(T,69875//102672)
  eb5 = -convert(T,2260//8211)
  eb6 = convert(T,1//4)

  ebtilde1 = -convert(T,31666707//9881966720)
  ebtilde3 =  convert(T,256875//105007616)
  ebtilde4 =  convert(T,2768025//128864768)
  ebtilde5 = -convert(T,169839//3864644)
  ebtilde6 =  convert(T,5247//225920)

  KenCarp4Tableau(γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,
                  btilde1,btilde3,btilde4,btilde5,btilde6,
                  c3,c4,c5,
                  α21,α31,α32,α41,α42,α51,α52,α53,α54,α61,α62,α63,α64,α65,
                  ea21,ea31,ea32,ea41,ea42,ea43,ea51,ea52,ea53,ea54,ea61,ea62,
                  ea63,ea64,ea65,eb1,eb3,eb4,eb5,eb6,ebtilde1,ebtilde3,ebtilde4,
                  ebtilde5,ebtilde6)
end

# Flip them all!

# ebtilde1 = 82889//524892 - 4586570599//29645900160
# ebtilde3 = 15625//83664 - 178811875//945068544
# ebtilde4 = 69875//102672 - 814220225//1159782912
# ebtilde5 = -2260//8211 + 3700637//11593932
# ebtilde6 = 1//4 - 61727//225920

struct Kvaerno5Tableau{T,T2}
  γ::T
  a31::T
  a32::T
  a41::T
  a42::T
  a43::T
  a51::T
  a52::T
  a53::T
  a54::T
  a61::T
  a63::T
  a64::T
  a65::T
  a71::T
  a73::T
  a74::T
  a75::T
  a76::T
  btilde1::T
  btilde3::T
  btilde4::T
  btilde5::T
  btilde6::T
  btilde7::T
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  α31::T
  α32::T
  α41::T
  α42::T
  α43::T
  α51::T
  α52::T
  α53::T
  α61::T
  α62::T
  α63::T
end

#=
# Kvaerno5
# Predict z3 from Hermite z2 and z1

c2 = 2γ
θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict others from z1 and z3 since it covers [0,1.23]

dt = c3 since interval is [c1,c3] and c1 = 0
θ =  c4/c3, c5/c3, c6/c3, c7/c3
z = dt*k

z₁ + Θ*(-4dt*z₁ - 2dt*z₃ - 6y₀ + Θ*(3dt*z₁ + 3z₃ + 6y₀ - 6y₁ ) + 6y₁)/dt

(1 + (-4θ + 3θ^2))*z₁ + (-2θ + 3θ^2)*z₃ + (6θ*(1-θ)/dt)*(y₁-y₀)

y₀ = uprev
y₁ = uprev + a31*z₁ + a32*z₂ + γ*z₃
y₁-y₀ = a31*z₁ + a32*z₂ + γ*z₃

(1 + (-4θ + 3θ^2) + a31*(6θ*(1-θ)/dt))*z₁ +
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))*z₃ + (6θ*(1-θ)/dt)*a32*z₂

dt = c3
θ = c4/c3
(1 + (-4θ + 3θ^2) + a31*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a32
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))

θ = c5/c3
(1 + (-4θ + 3θ^2) + a31*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a32
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))

θ = c6/c3
(1 + (-4θ + 3θ^2) + a31*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a32
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))
=#

function Kvaerno5Tableau(T,T2)
  γ = convert(T,0.26)
  a31 = convert(T,0.13)
  a32 = convert(T,0.84033320996790809)
  a41 = convert(T,0.22371961478320505)
  a42 = convert(T,0.47675532319799699)
  a43 = -convert(T,0.06470895363112615)
  a51 = convert(T,0.16648564323248321)
  a52 = convert(T,0.10450018841591720)
  a53 = convert(T,0.03631482272098715)
  a54 = -convert(T,0.13090704451073998)
  a61 = convert(T,0.13855640231268224)
  a63 = -convert(T,0.04245337201752043)
  a64 = convert(T,0.02446657898003141)
  a65 = convert(T,0.61943039072480676)
  a71 = convert(T,0.13659751177640291)
  a73 = -convert(T,0.05496908796538376)
  a74 = -convert(T,0.04118626728321046)
  a75 = convert(T,0.62993304899016403)
  a76 = convert(T,0.06962479448202728)
  btilde1 = convert(T,0.00195889053627933) # a61-a71
  btilde3 = convert(T,0.01251571594786333) # a63-a73
  btilde4 = convert(T,0.06565284626324187) # a64-a74
  btilde5 = -convert(T,0.01050265826535727) # a65-a75
  btilde6 = convert(T,0.19037520551797272) # γ-a76
  btilde7 = -γ
  α21 = convert(T,2) # c2/γ
  α31 = convert(T,-1.366025403784441)
  α32 = convert(T,2.3660254037844357)
  α41 = convert(T,-0.19650552613122207)
  α42 = convert(T,0.8113579546496623)
  α43 = convert(T,0.38514757148155954)
  α51 = convert(T,0.10375304369958693)
  α52 = convert(T,0.937994698066431)
  α53 = convert(T,-0.04174774176601781)
  α61 = convert(T,-0.17281112873898072)
  α62 = convert(T,0.6235784481025847)
  α63 = convert(T,0.5492326806363959)
  c3 = convert(T,1.230333209967908)
  c4 = convert(T,0.895765984350076)
  c5 = convert(T,0.436393609858648)
  c6 = convert(T,1)
  Kvaerno5Tableau(γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,
                  a61,a63,a64,a65,a71,a73,a74,a75,a76,
                  btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,
                  c3,c4,c5,c6,α31,α32,α41,α42,α43,α51,α52,α53,
                  α61,α62,α63)
end

struct KenCarp5Tableau{T,T2}
  γ::T
  a31::T
  a32::T
  a41::T
  a43::T
  a51::T
  a53::T
  a54::T
  a61::T
  a63::T
  a64::T
  a65::T
  a71::T
  a73::T
  a74::T
  a75::T
  a76::T
  a81::T
  a84::T
  a85::T
  a86::T
  a87::T
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  α31::T
  α32::T
  α41::T
  α42::T
  α51::T
  α52::T
  α61::T
  α62::T
  α71::T
  α72::T
  α73::T
  α74::T
  α75::T
  α81::T
  α82::T
  α83::T
  α84::T
  α85::T
  btilde1::T
  btilde4::T
  btilde5::T
  btilde6::T
  btilde7::T
  btilde8::T
  ea21::T
  ea31::T
  ea32::T
  ea41::T
  ea43::T
  ea51::T
  ea53::T
  ea54::T
  ea61::T
  ea63::T
  ea64::T
  ea65::T
  ea71::T
  ea73::T
  ea74::T
  ea75::T
  ea76::T
  ea81::T
  ea83::T
  ea84::T
  ea85::T
  ea86::T
  ea87::T
  eb1::T
  eb4::T
  eb5::T
  eb6::T
  eb7::T
  eb8::T
  ebtilde1::T
  ebtilde4::T
  ebtilde5::T
  ebtilde6::T
  ebtilde7::T
  ebtilde8::T
end

#=
# KenCarp5
# Predict z3 from Hermite z2 and z1

c2 = 2γ
θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z4 from z2 and z1
θ = c4/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z5 from z2 and z1
θ = c5/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z6 from z2 and z1
θ = c6/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)

# Predict z7 from z5 and z1
θ = c7/c5
dt = c5

(1 + (-4θ + 3θ^2))*z₁ + (-2θ + 3θ^2)*z₃ + (6θ*(1-θ)/dt)*(y₁-y₀)
y₁-y₀ = a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄ + γ*z₅

(1 + (-4θ + 3θ^2) + a51*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a52
(6θ*(1-θ)/dt)*a53
(6θ*(1-θ)/dt)*a54
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))

# Predict z8 from z5 and z1
θ = 1/c5
dt = c5

(1 + (-4θ + 3θ^2) + a51*(6θ*(1-θ)/dt))
(6θ*(1-θ)/dt)*a52
(6θ*(1-θ)/dt)*a53
(6θ*(1-θ)/dt)*a54
(-2θ + 3θ^2 + γ*(6θ*(1-θ)/dt))

=#

function KenCarp5Tableau(T,T2)
  γ = convert(T,41//200)
  a31 = convert(T,41//400)
  a32 = -convert(T,567603406766//11931857230679)
  a41 = convert(T,683785636431//9252920307686)
  a43 = -convert(T,110385047103//1367015193373)
  a51 = convert(T,3016520224154//10081342136671)
  a53 = convert(T,30586259806659//12414158314087)
  a54 = -convert(T,22760509404356//11113319521817)
  a61 = convert(T,218866479029//1489978393911)
  a63 = convert(T,638256894668//5436446318841)
  a64 = -convert(T,1179710474555//5321154724896)
  a65 = -convert(T,60928119172//8023461067671)
  a71 = convert(T,1020004230633//5715676835656)
  a73 = convert(T,25762820946817//25263940353407)
  a74 = -convert(T,2161375909145//9755907335909)
  a75 = -convert(T,211217309593//5846859502534)
  a76 = -convert(T,4269925059573//7827059040749)
  a81 = -convert(T,872700587467//9133579230613)
  a84 = convert(T,22348218063261//9555858737531)
  a85 = -convert(T,1143369518992//8141816002931)
  a86 = -convert(T,39379526789629//19018526304540)
  a87 = convert(T,32727382324388//42900044865799)
  # bhat1 = -convert(T,975461918565//9796059967033)
  # bhat4 = convert(T,78070527104295//32432590147079)
  # bhat5 = -convert(T,548382580838//3424219808633)
  # bhat6 = -convert(T,33438840321285//15594753105479)
  # bhat7 = convert(T,3629800801594//4656183773603)
  # bhat8 = convert(T,4035322873751//18575991585200)
  btilde1 = -convert(T,360431431567533808054934//89473089856732078284381229) # bhat1-a81
  btilde4 = convert(T,21220331609936495351431026//309921249937726682547321949) # bhat4-a84
  btilde5 = -convert(T,42283193605833819490634//2144566741190883522084871) # bhat5-a85
  btilde6 = -convert(T,21843466548811234473856609//296589222149359214696574660) # bhat6-a86
  btilde7 = convert(T,3333910710978735057753642//199750492790973993533703797) # bhat7-a87
  btilde8 = convert(T,45448919757//3715198317040) # bhat8-γ
  c3 = convert(T2,2935347310677//11292855782101)
  c4 = convert(T2,1426016391358//7196633302097)
  c5 = convert(T2,92//100)
  c6 = convert(T2,24//100)
  c7 = convert(T2,3//5)
  α31 = convert(T,169472355998441//463007087066141)
  α32 = convert(T,293534731067700//463007087066141)
  α41 = convert(T,152460326250177//295061965385977)
  α42 = convert(T,142601639135800//295061965385977)
  α51 = convert(T,-51//41)
  α52 = convert(T,92//41)
  α61 = convert(T,17//41)
  α62 = convert(T,24//41)
  α71 = convert(T,13488091065527792//122659689776876057)
  α72 = convert(T,-3214953045//3673655312)
  α73 = convert(T,550552676519862000//151043064207496529)
  α74 = convert(T,-409689169278408000//135215758621947439)
  α75 = convert(T,3345//12167)
  α81 = convert(T,1490668709762032//122659689776876057)
  α82 = convert(T,5358255075//14694621248)
  α83 = convert(T,-229396948549942500//151043064207496529)
  α84 = convert(T,170703820532670000//135215758621947439)
  α85 = convert(T,30275//24334)

  ea21 = convert(T,41//100)
  ea31 = convert(T,367902744464//2072280473677)
  ea32 = convert(T,677623207551//8224143866563)
  ea41 = convert(T,1268023523408//10340822734521)
  ea43 = convert(T,1029933939417//13636558850479)
  ea51 = convert(T,14463281900351//6315353703477)
  ea53 = convert(T,66114435211212//5879490589093)
  ea54 = -convert(T,54053170152839//4284798021562)
  ea61 = convert(T,14090043504691//34967701212078)
  ea63 = convert(T,15191511035443//11219624916014)
  ea64 = -convert(T,18461159152457//12425892160975)
  ea65 = -convert(T,281667163811//9011619295870)
  ea71 = convert(T,19230459214898//13134317526959)
  ea73 = convert(T,21275331358303//2942455364971)
  ea74 = -convert(T,38145345988419//4862620318723)
  ea75 = -convert(T,1//8)
  ea76 = -convert(T,1//8)
  ea81 = -convert(T,19977161125411//11928030595625)
  ea83 = -convert(T,40795976796054//6384907823539)
  ea84 = convert(T,177454434618887//12078138498510)
  ea85 = convert(T,782672205425//8267701900261)
  ea86 = -convert(T,69563011059811//9646580694205)
  ea87 = convert(T,7356628210526//4942186776405)

  eb1 = -convert(T,872700587467//9133579230613)
  eb4 = convert(T,22348218063261//9555858737531)
  eb5 = -convert(T,1143369518992//8141816002931)
  eb6 = -convert(T,39379526789629//19018526304540)
  eb7 = convert(T,32727382324388//42900044865799)
  eb8 = convert(T,41//200)

  ebtilde1 = -convert(T,360431431567533808054934//89473089856732078284381229)
  ebtilde4 = convert(T,21220331609936495351431026//309921249937726682547321949)
  ebtilde5 = -convert(T,42283193605833819490634//2144566741190883522084871)
  ebtilde6 = -convert(T,21843466548811234473856609//296589222149359214696574660)
  ebtilde7 = convert(T,3333910710978735057753642//199750492790973993533703797)
  ebtilde8 = convert(T,45448919757//3715198317040)

  KenCarp5Tableau(γ,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,
                  a71,a73,a74,a75,a76,a81,a84,a85,a86,a87,
                  c3,c4,c5,c6,c7,α31,α32,α41,α42,α51,α52,
                  α61,α62,α71,α72,α73,α74,α75,α81,α82,α83,α84,α85,
                  btilde1,btilde4,btilde5,btilde6,btilde7,btilde8,
                  ea21,ea31,ea32,ea41,ea43,ea51,ea53,ea54,ea61,ea63,
                  ea64,ea65,ea71,ea73,ea74,ea75,ea76,ea81,ea83,ea84,
                  ea85,ea86,ea87,eb1,eb4,eb5,eb6,eb7,eb8,ebtilde1,
                  ebtilde4,ebtilde5,ebtilde6,ebtilde7,ebtilde8)
end

# Flip them all!

# ebtilde1 = -Int128(872700587467)//9133579230613 + Int128(975461918565)//9796059967033
# ebtilde4 = Int128(22348218063261)//9555858737531 - Int128(78070527104295)//32432590147079
# ebtilde5 = -Int128(1143369518992)//8141816002931 + Int128(548382580838)//3424219808633
# ebtilde6 = -Int128(39379526789629)//19018526304540 + Int128(33438840321285)//15594753105479
# ebtilde7 = Int128(32727382324388)//42900044865799 - Int128(3629800801594)//4656183773603
# ebtilde8 = Int128(41)//200 - Int128(4035322873751)//18575991585200

struct ESDIRK54I8L2SATableau{T,T2}
  γ::T
  a31::T; a32::T
  a41::T; a42::T; a43::T
  a51::T; a52::T; a53::T; a54::T
  a61::T; a62::T; a63::T; a64::T; a65::T
  a71::T; a72::T; a73::T; a74::T; a75::T; a76::T
  a81::T; a82::T; a83::T; a84::T; a85::T; a86::T; a87::T
  c3::T2; c4::T2; c5::T2; c6::T2; c7::T2
  btilde1::T; btilde2::T; btilde3::T; btilde4::T; btilde5::T; btilde6::T; btilde7::T; btilde8::T
end

function ESDIRK54I8L2SATableau(T, T2)
  γ  = convert(T, 1//4)
  a31 = convert(T, 1748874742213//5795261096931)
  a32 = convert(T, 1748874742213//5795261096931)
  a41 = convert(T, 2426486750897//12677310711630)
  a42 = convert(T, 2426486750897//12677310711630)
  a43 = convert(T, -783385356511//7619901499812)
  a51 = convert(T, 1616209367427//5722977998639)
  a52 = convert(T, 1616209367427//5722977998639)
  a53 = convert(T, -211896077633//5134769641545)
  a54 = convert(T, 464248917192//17550087120101)
  a61 = convert(T, 1860464898611//7805430689312)
  a62 = convert(T, 1825204367749//7149715425471)
  a63 = convert(T, -1289376786583//6598860380111)
  a64 = convert(T, 55566826943//2961051076052)
  a65 = convert(T, 1548994872005//13709222415197)
  a71 = convert(T, 1783640092711//14417713428467)
  a72 = convert(T, -5781183663275//18946039887294)
  a73 = convert(T, -57847255876685//10564937217081)
  a74 = convert(T, 29339178902168//9787613280015)
  a75 = convert(T, 122011506936853//12523522131766)
  a76 = convert(T, -60418758964762//9539790648093)
  a81 = convert(T, 3148564786223//23549948766475)
  a82 = convert(T, -4152366519273//20368318839251)
  a83 = convert(T, -143958253112335//33767350176582)
  a84 = convert(T, 16929685656751//6821330976083)
  a85 = convert(T, 37330861322165//4907624269821)
  a86 = convert(T, -103974720808012//20856851060343)
  a87 = convert(T, -93596557767//4675692258479)
  c3  = convert(T2, (2+sqrt(convert(T, 2)))/4)
  c4  = convert(T2, 53//100)
  c5  = convert(T2, 4//5)
  c6  = convert(T2, 17//25)
  c7  = convert(T2, 1)
  btilde1 = convert(T, -206948709334490044469698//10480001459192469358387375)
  btilde2 = convert(T, -38800234036520698435148405//193732560740235652447264213)
  btilde3 = convert(T, -42312118141829119927717945//17651296211698462951718982)
  btilde4 = convert(T, 77601425937783402908082927//76091823354023603930374324)
  btilde5 = convert(T, 7549135156215231570800855//1787280796764804347348433)
  btilde6 = convert(T, -401514321964993460839314379//150599863859115530598736650)
  btilde7 = convert(T, 17761325247710183915293664//33262552787523086832167825)
  btilde8 = convert(T, -25249389576073//51072051291964)
  ESDIRK54I8L2SATableau(γ,
                        a31, a32,
                        a41, a42, a43,
                        a51, a52, a53, a54,
                        a61, a62, a63, a64, a65,
                        a71, a72, a73, a74, a75, a76,
                        a81, a82, a83, a84, a85, a86, a87,
                                  c3,  c4,  c5,  c6,  c7,
                        btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, btilde8)
end

struct SDIRK22Tableau{T}
  a::T
  α::T
  β::T
end

function SDIRK22Tableau(T)
  a = convert(T, 1-1/sqrt(2))
  α = convert(T, -sqrt(2))
  β = convert(T, 1+sqrt(2))
  SDIRK22Tableau(a, α, β)
end
