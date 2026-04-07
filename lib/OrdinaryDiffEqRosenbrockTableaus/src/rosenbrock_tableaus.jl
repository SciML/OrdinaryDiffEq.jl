const RODAS4A = [
    0 0 0 0 0 0
    1.544 0 0 0 0 0
    0.9466785280815826 0.2557011698983284 0 0 0 0
    3.314825187068521 2.896124015972201 0.9986419139977817 0 0 0
    1.221224509226641 6.019134481288629 12.53708332932087 -0.687886036105895 0 0
    1.221224509226641 6.019134481288629 12.53708332932087 -0.687886036105895 1 0
]
const RODAS4C = [
    0 0 0 0 0
    -5.6688 0 0 0 0
    -2.430093356833875 -0.2063599157091915 0 0 0
    -0.1073529058151375 -9.594562251023355 -20.47028614809616 0 0
    7.496443313967647 -10.24680431464352 -33.99990352819905 11.7089089320616 0
    8.083246795921522 -7.981132988064893 -31.52159432874371 16.31930543123136 -6.058818238834054
]
const RODAS4c = [0, 0.386, 0.21, 0.63, 1, 1]
const RODAS4d = [0.25, -0.1043, 0.1035, -0.0362, 0, 0]
const RODAS4H = [
    10.12623508344586 -7.487995877610167 -34.80091861555747 -7.992771707568823 1.025137723295662 0
    -0.6762803392801253 6.087714651680015 16.43084320892478 24.76722511418386 -6.594389125716872 0
]
"""
    Rodas4Tableau(T, T2)

The tableau for the 4th order L-stable Rosenbrock method Rodas4. 
It is a 6-stage method with a built-in error estimate.
Reference: Hairer, E., Nørsett, S. P., & Wanner, G. (1996). Solving Ordinary Differential Equations II.
"""
function Rodas4Tableau(T, T2)
    gamma = 0.25
    b = T[RODAS4A[6, 1], RODAS4A[6, 2], RODAS4A[6, 3], RODAS4A[6, 4], RODAS4A[6, 5], one(T)]
    btilde = T[zero(T), zero(T), zero(T), zero(T), zero(T), one(T)]
    return RodasTableau{T, T2, Vector{T}}(RODAS4A, RODAS4C, gamma, RODAS4c, RODAS4d, RODAS4H, b, btilde)
end

const RODAS42A = [
    0 0 0 0 0 0
    1.4028884 0 0 0 0 0
    0.6581212688557198 -1.320936088384301 0 0 0 0
    7.131197445744498 16.02964143958207 -5.561572550509766 0 0 0
    22.73885722420363 67.38147284535289 -31.2187749303856 0.7285641833203814 0 0
    22.73885722420363 67.38147284535289 -31.2187749303856 0.7285641833203814 1 0
]
const RODAS42C = [
    0 0 0 0 0
    -5.1043536 0 0 0 0
    -2.899967805418783 4.040399359702244 0 0 0
    -32.64449927841361 -99.35311008728094 49.99119122405989 0 0
    -76.46023087151691 -278.5942120829058 153.9294840910643 10.97101866258358 0
    -76.29701586804983 -294.2795630511232 162.0029695867566 23.6516690309527 -7.652977706771382
]
const RODAS42c = [0, 0.3507221, 0.2557041, 0.681779, 1, 1]
const RODAS42d = [0.25, -0.0690221, -0.0009672, -0.087979, 0, 0]
const RODAS42H = [
    -38.71940424117216 -135.8025833007622 64.51068857505875 -4.192663174613162 -2.53193205033506 0
    -14.99268484949843 -76.30242396627033 58.65928432851416 16.61359034616402 -0.6758691794084156 0
]
"""
    Rodas42Tableau(T, T2)

A 4th order L-stable Rosenbrock method with 6 stages, often used as an alternative to Rodas4.
Reference: Hairer, E., & Wanner, G. (1996). Solving Ordinary Differential Equations II:
    Stiff and Differential-Algebraic Problems. Springer-Verlag, 2nd Edition.
"""
function Rodas42Tableau(T, T2)
    gamma = 0.25
    b = T[RODAS42A[6, 1], RODAS42A[6, 2], RODAS42A[6, 3], RODAS42A[6, 4], RODAS42A[6, 5], one(T)]
    btilde = T[zero(T), zero(T), zero(T), zero(T), zero(T), one(T)]
    return RodasTableau{T, T2, Vector{T}}(RODAS42A, RODAS42C, gamma, RODAS42c, RODAS42d, RODAS42H, b, btilde)
end

const RODAS4PA = [
    0 0 0 0 0 0
    3 0 0 0 0 0
    1.831036793486759 0.4955183967433795 0 0 0 0
    2.304376582692669 -0.05249275245743001 -1.176798761832782 0 0 0
    -7.170454962423024 -4.741636671481785 -16.31002631330971 -1.062004044111401 0 0
    -7.170454962423024 -4.741636671481785 -16.31002631330971 -1.062004044111401 1 0
]
const RODAS4PC = [
    0 0 0 0 0
    -12 0 0 0 0
    -8.791795173947035 -2.207865586973518 0 0 0
    10.81793056857153 6.780270611428266 19.5348594464241 0 0
    34.19095006749676 15.49671153725963 54.7476087596413 14.16005392148534 0
    34.62605830930532 15.30084976114473 56.99955578662667 18.40807009793095 -5.714285714285717
]
const RODAS4Pc = [0, 0.75, 0.21, 0.63, 1, 1]
const RODAS4Pd = [0.25, -0.5, -0.023504, -0.0362, 0, 0]
const RODAS4PH = [
    25.09876703708589 11.62013104361867 28.49148307714626 -5.664021568594133 0 0
    1.638054557396973 -0.7373619806678748 8.47791821923899 15.9925314877952 -1.882352941176471 0
]
"""
    Rodas4PTableau(T, T2)

A 4th order L-stable Rosenbrock method with 6 stages, emphasizing stability for parabolic problems.
Reference: Steinebach, G. (1995). Order-reduction of ROW-methods for DAEs and
    method of lines applications. Preprint-Nr. 1741, FB Mathematik, TH Darmstadt.
"""
function Rodas4PTableau(T, T2)
    gamma = 0.25
    b = T[RODAS4PA[6, 1], RODAS4PA[6, 2], RODAS4PA[6, 3], RODAS4PA[6, 4], RODAS4PA[6, 5], one(T)]
    btilde = T[zero(T), zero(T), zero(T), zero(T), zero(T), one(T)]
    return RodasTableau{T, T2, Vector{T}}(RODAS4PA, RODAS4PC, gamma, RODAS4Pc, RODAS4Pd, RODAS4PH, b, btilde)
end

const RODAS4P2A = [
    0 0 0 0 0 0
    3 0 0 0 0 0
    0.906377755268814 -0.189707390391685 0 0 0 0
    3.758617027739064 1.161741776019525 -0.849258085312803 0 0 0
    7.089566927282776 4.573591406461604 -8.423496976860259 -0.959280113459775 0 0
    7.089566927282776 4.573591406461604 -8.423496976860259 -0.959280113459775 1 0
]
const RODAS4P2C = [
    0 0 0 0 0
    -12 0 0 0 0
    -6.354581592719008 0.338972550544623 0 0 0
    -8.575016317114033 -7.606483992117508 12.22499765012482 0 0
    -5.888975457523102 -8.157396617841821 24.805546872612922 12.790401512796979 0
    -4.408651676063871 -6.692003137674639 24.625568527593117 16.627521966636085 -5.714285714285718
]
const RODAS4P2c = [0, 0.75, 0.321448134013046, 0.519745732277726, 1, 1]
const RODAS4P2d = [0.25, -0.5, -0.189532918363016, 0.085612108792769, 0, 0]
const RODAS4P2H = [
    -5.323528268423303 -10.042123754867493 17.175254928256965 -5.079931171878093 -0.016185991706112 0
    6.984505741529879 6.914061169603662 -0.849178943070653 18.104410789349338 -3.516963011559032 0
]
"""
    Rodas4P2Tableau(T, T2)

An improved version of the Rodas4P 4th order L-stable Rosenbrock method.
Reference: Steinebach, G. (2020). Improvement of Rosenbrock-Wanner method RODASP.
    In: Progress in Differential-Algebraic Equations II, pp. 165–184.
    Springer, Cham.
"""
function Rodas4P2Tableau(T, T2)
    gamma = 0.25
    b = T[RODAS4P2A[6, 1], RODAS4P2A[6, 2], RODAS4P2A[6, 3], RODAS4P2A[6, 4], RODAS4P2A[6, 5], one(T)]
    btilde = T[zero(T), zero(T), zero(T), zero(T), zero(T), one(T)]
    return RodasTableau{T, T2, Vector{T}}(RODAS4P2A, RODAS4P2C, gamma, RODAS4P2c, RODAS4P2d, RODAS4P2H, b, btilde)
end

################################################################################
# Rodas5 family (8-stage)
################################################################################

const RODAS5A = [
    0 0 0 0 0 0 0 0
    2.0 0 0 0 0 0 0 0
    3.040894194418781 1.041747909077569 0 0 0 0 0 0
    2.576417536461461 1.62208306077664 -0.9089668560264532 0 0 0 0 0
    2.760842080225597 1.446624659844071 -0.3036980084553738 0.2877498600325443 0 0 0 0
    -14.09640773051259 6.925207756232704 -41.47510893210728 2.343771018586405 24.13215229196062 0 0 0
    -14.09640773051259 6.925207756232704 -41.47510893210728 2.343771018586405 24.13215229196062 1 0 0
    -14.09640773051259 6.925207756232704 -41.47510893210728 2.343771018586405 24.13215229196062 1 1 0
]
const RODAS5C = [
    0 0 0 0 0 0 0
    -10.31323885133993 0 0 0 0 0 0
    -21.04823117650003 -7.234992135176716 0 0 0 0 0
    32.22751541853323 -4.943732386540191 19.44922031041879 0 0 0 0
    -20.69865579590063 -8.816374604402768 1.260436877740897 -0.7495647613787146 0 0 0
    -46.22004352711257 -17.49534862857472 -289.6389582892057 93.60855400400906 318.3822534212147 0 0
    34.20013733472935 -14.1553540271769 57.823356409884 25.83362985412365 1.408950972071624 -6.551835421242162 0
    42.57076742291101 -13.80770672017997 93.98938432427124 18.77919633714503 -31.5835918722337 -6.685968952921985 -5.810979938412932
]
const RODAS5c = [0, 0.38, 0.3878509998321533, 0.483971893787384, 0.457047700881958, 1, 1, 1]
const RODAS5d = [
    0.19, -0.1823079225333714636, -0.319231832186874912,
    0.3449828624725343, -0.377417564392089818, 0, 0, 0,
]

const RODAS5H = [
    27.354592673333357 -6.925207756232857 26.40037733258859 0.5635230501052979 -4.699151156849391 -1.6008677469422725 -1.5306074446748028 -1.3929872940716344
    44.19024239501722 1.3677947663381929e-13 202.93261852171622 -35.5669339789154 -181.91095152160645 3.4116351403665033 2.5793540257308067 2.2435122582734066
    -44.0988150021747 -5.755396159656812e-13 -181.26175034586677 56.99302194811676 183.21182741427398 -7.480257918273637 -5.792426076169686 -5.32503859794143
]

"""
    Rodas5Tableau(T, T2)

The tableau for the 5th order L-stable Rosenbrock method Rodas5.
It is an 8-stage method designed for high-accuracy stiff integration.
Reference: Di Marzo, G. (1993). Rodas5(4) -- Méthodes de Rosenbrock d'ordre 5(4)
    adaptées aux problemes différentiels-algébriques. Master's thesis,
    University of Geneva.
"""
function Rodas5Tableau(T, T2)

function ROS3PRodasTableau(T, T2)
    gamma = convert(T2, 1 / 2 + sqrt(3) / 6)
    igamma = inv(gamma)
    a21 = convert(T, igamma)
    a31 = convert(T, igamma)
    C21 = convert(T, -igamma^2)
    tmp = -igamma * (convert(T, 2) - convert(T, 1 / 2) * igamma)
    C31 = -igamma * (convert(T, 1) - tmp)
    C32 = tmp
    tmp2 = igamma * (convert(T, 2 / 3) - convert(T, 1 / 6) * igamma)
    b1 = igamma * (convert(T, 1) + tmp2)
    b2 = tmp2
    b3 = convert(T, 1 / 3) * igamma
    btilde1 = b1 - convert(T, 2.113248654051871)
    btilde2 = b2 - convert(T, 1.0)
    btilde3 = b3 - convert(T, 0.4226497308103742)
    d1 = convert(T, 0.7886751345948129)
    d2 = convert(T, -0.2113248654051871)
    d3 = convert(T, -1.077350269189626)

    A = zeros(T, 3, 3)
    A[2, 1] = a21
    A[3, 1] = a31

    C = zeros(T, 3, 3)
    C[2, 1] = C21
    C[3, 1] = C31
    C[3, 2] = C32

    c = T2[convert(T2, 0), convert(T2, 1), convert(T2, 1)]
    d = T[d1, d2, d3]
    b = T[b1, b2, b3]
    btilde = T[btilde1, btilde2, btilde3]
    H = zeros(T, 0, 3)

    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

################################################################################
# Rodas3 (4-stage, hand-written)
################################################################################

"""
    Rodas3RodasTableau(T, T2)

A 3rd order Rosenbrock method with 4 stages.
Reference: Sandu, A., et al. (1997). Benchmarking stiff ode solvers for atmospheric chemistry problems-I. implicit vs explicit. Atmospheric Environment, 31(19), 3151-3166.
"""
function Rodas3RodasTableau(T, T2)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 0)
    A[3, 1] = convert(T, 2)
    A[4, 1] = convert(T, 2)
    A[4, 3] = convert(T, 1)

    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, 4)
    C[3, 1] = convert(T, 1)
    C[3, 2] = convert(T, -1)
    C[4, 1] = convert(T, 1)
    C[4, 2] = convert(T, -1)
    C[4, 3] = convert(T, -8 // 3)

    gamma = convert(T2, 1 // 2)
    c = T2[convert(T2, 0), convert(T2, 0), convert(T2, 1), convert(T2, 1)]
    d = T[convert(T, 1 // 2), convert(T, 3 // 2), zero(T), zero(T)]
    b = T[convert(T, 2), zero(T), convert(T, 1), convert(T, 1)]
    btilde = T[zero(T), zero(T), zero(T), one(T)]
    H = zeros(T, 0, 4)

    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

################################################################################
# Rodas3P (5-stage, with H matrix for dense output)
################################################################################

"""
    Rodas3PRodasTableau(T, T2)

A 3rd order Rosenbrock method with 5 stages, including a dense output matrix H.
Reference: Steinebach, G. (2024). Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications. Proceedings of the JuliaCon Conferences.
"""
function Rodas3PRodasTableau(T, T2)

function RosShamp4RodasTableau(T, T2)
    gamma = convert(T2, 1 // 2)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 1.92)
    A[3, 2] = convert(T, 0.24)
    A[4, 1] = convert(T, 1.92)
    A[4, 2] = convert(T, 0.24)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -8 // 1)
    C[3, 1] = convert(T, 372 // 25)
    C[3, 2] = convert(T, 12 // 5)
    C[4, 1] = convert(T, -112 // 125)
    C[4, 2] = convert(T, -54 // 125)
    C[4, 3] = convert(T, -2 // 5)
    c = T2[0.0, 1.0, 0.6, 0.6]
    d = T[1 // 2, -3 // 2, 121 // 50, 29 // 250]
    b = T[19 // 9, 1 // 2, 25 // 108, 125 // 108]
    btilde = T[17 // 54, 7 // 36, 0 // 1, 125 // 108]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    Veldd4RodasTableau(T, T2)

A 4th order Rosenbrock method by van Veldhuizen.
"""
function Veldd4RodasTableau(T, T2)
    gamma = convert(T2, 0.2257081148225682)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 4.812234362695436)
    A[3, 2] = convert(T, 4.578146956747842)
    A[4, 1] = convert(T, 4.812234362695436)
    A[4, 2] = convert(T, 4.578146956747842)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -5.333333333333331)
    C[3, 1] = convert(T, 6.100529678848254)
    C[3, 2] = convert(T, 1.804736797378427)
    C[4, 1] = convert(T, -2.540515456634749)
    C[4, 2] = convert(T, -9.443746328915205)
    C[4, 3] = convert(T, -1.988471753215993)
    c = T2[0.0, 0.4514162296451364, 0.8755928946018455, 0.8755928946018455]
    d = T[0.2257081148225682, -0.04599403502680582, 0.5177590504944076, -0.03805623938054428]
    b = T[4.289339254654537, 5.036098482851414, 0.6085736420673917, 1.355958941201148]
    btilde = T[2.175672787531755, 2.950911222575741, -0.785974454488743, -1.355958941201148]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

function Velds4RodasTableau(T, T2)
    gamma = convert(T2, 1 // 2)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 1.75)
    A[3, 2] = convert(T, 0.25)
    A[4, 1] = convert(T, 1.75)
    A[4, 2] = convert(T, 0.25)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -8 // 1)
    C[3, 1] = convert(T, -8 // 1)
    C[3, 2] = convert(T, -1 // 1)
    C[4, 1] = convert(T, 1 // 2)
    C[4, 2] = convert(T, -1 // 2)
    C[4, 3] = convert(T, 2 // 1)
    c = T2[0.0, 1.0, 0.5, 0.5]
    d = T[1 // 2, -3 // 2, -3 // 4, 1 // 4]
    b = T[4 // 3, 2 // 3, -4 // 3, 4 // 3]
    btilde = T[-1 // 3, -1 // 3, 0 // 1, -4 // 3]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    GRK4TRodasTableau(T, T2)

A 4th order Generalized Runge-Kutta (Rosenbrock) method by Kaps and Rentrop.
"""
function GRK4TRodasTableau(T, T2)
    gamma = convert(T2, 0.231)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 4.524708207373116)
    A[3, 2] = convert(T, 4.163528788597648)
    A[4, 1] = convert(T, 4.524708207373116)
    A[4, 2] = convert(T, 4.163528788597648)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -5.071675338776316)
    C[3, 1] = convert(T, 6.020152728650786)
    C[3, 2] = convert(T, 0.1597506846727117)
    C[4, 1] = convert(T, -1.856343618686113)
    C[4, 2] = convert(T, -8.505380858179826)
    C[4, 3] = convert(T, -2.084075136023187)
    c = T2[0.0, 0.462, 0.8802083333333334, 0.8802083333333334]
    d = T[0.231, -0.03962966775244303, 0.5507789395789127, -0.05535098457052764]
    b = T[3.957503746640777, 4.624892388363313, 0.6174772638750108, 1.282612945269037]
    btilde = T[2.302155402932996, 3.073634485392623, -0.8732808018045032, -1.282612945269037]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

function GRK4ARodasTableau(T, T2)
    gamma = convert(T2, 0.395)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 1.108860759493671)
    A[3, 1] = convert(T, 2.37708526198336)
    A[3, 2] = convert(T, 0.1850114988899692)
    A[4, 1] = convert(T, 2.37708526198336)
    A[4, 2] = convert(T, 0.1850114988899692)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -4.920188402397641)
    C[3, 1] = convert(T, 1.055588686048583)
    C[3, 2] = convert(T, 3.351817267668938)
    C[4, 1] = convert(T, 3.846869007049313)
    C[4, 2] = convert(T, 3.42710924126818)
    C[4, 3] = convert(T, -2.162408848753263)
    c = T2[0.0, 0.438, 0.87, 0.87]
    d = T[0.395, -0.372672395484092, 0.06629196544571492, 0.4340946962568634]
    b = T[1.84568324040584, 0.1369796894360503, 0.7129097783291559, 0.6329113924050632]
    btilde = T[0.04831870177201765, -0.6471108651049505, 0.218687666050024, -0.6329113924050632]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    Ros4LStabRodasTableau(T, T2)

A 4th order L-stable Rosenbrock method.
"""
function Ros4LStabRodasTableau(T, T2)
    gamma = convert(T2, 0.57282)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 1.867943637803922)
    A[3, 2] = convert(T, 0.2344449711399156)
    A[4, 1] = convert(T, 1.867943637803922)
    A[4, 2] = convert(T, 0.2344449711399156)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -7.13761503641231)
    C[3, 1] = convert(T, 2.580708087951457)
    C[3, 2] = convert(T, 0.6515950076447975)
    C[4, 1] = convert(T, -2.137148994382534)
    C[4, 2] = convert(T, -0.3214669691237626)
    C[4, 3] = convert(T, -0.6949742501781779)
    c = T2[0.0, 1.14564, 0.65521686381559, 0.65521686381559]
    d = T[0.57282, -1.769193891319233, 0.7592633437920482, -0.104902108710045]
    b = T[2.255570073418735, 0.2870493262186792, 0.435317943184018, 1.093502252409163]
    btilde = T[-0.2815431932141155, -0.0727619912493892, -0.1082196201495311, -1.093502252409163]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS2RodasTableau(T, T2)

A 2nd order Rosenbrock method.
"""
function ROS2RodasTableau(T, T2)
    gamma = convert(T2, 1.7071067811865475)
    A = zeros(T, 2, 2)
    A[2, 1] = convert(T, 0.585786437626905)
    C = zeros(T, 2, 2)
    C[2, 1] = convert(T, -1.17157287525381)
    c = T2[0.0, 1.0]
    d = T[1.7071067811865475, -1.7071067811865475]
    b = T[0.8786796564403574, 0.2928932188134525]
    btilde = T[0.2928932188134525, 0.2928932188134525]
    H = zeros(T, 0, 2)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

function ROS2PRRodasTableau(T, T2)
    gamma = convert(T2, 0.228155493653962)
    A = zeros(T, 3, 3)
    A[2, 1] = convert(T, 4.382975767906234)
    A[3, 1] = convert(T, 4.382975767906234)
    A[3, 2] = convert(T, 4.382975767906234)
    C = zeros(T, 3, 3)
    C[2, 1] = convert(T, -4.382975767906234)
    C[3, 1] = convert(T, -4.382975767906234)
    C[3, 2] = convert(T, -16.827500814147)
    c = T2[0.0, 1.0, 1.0]
    d = T[0.228155493653962, 0.0, -2.7755575615628914e-17]
    b = T[4.382975767906234, 4.382975767906234, 1.0]
    btilde = T[-9.968705307220848e-18, 3.3829757679062333, 1.0]
    H = zeros(T, 0, 3)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS2SRodasTableau(T, T2)

A 2nd order Rosenbrock method with specific stability properties.
"""
function ROS2SRodasTableau(T, T2)
    gamma = convert(T2, 0.292893218813452)
    A = zeros(T, 3, 3)
    A[2, 1] = convert(T, 2.0000000000000036)
    A[3, 1] = convert(T, 6.828427124746214)
    A[3, 2] = convert(T, 3.4142135623731007)
    C = zeros(T, 3, 3)
    C[2, 1] = convert(T, -6.828427124746214)
    C[3, 1] = convert(T, -10.949747468305889)
    C[3, 2] = convert(T, -7.535533905932761)
    c = T2[0.0, 0.585786437626905, 1.0]
    d = T[0.292893218813452, -0.292893218813453, -5.551115123125783e-17]
    b = T[6.828427124746214, 3.414213562373101, 1.0]
    btilde = T[-0.23570226039551292, -0.23570226039551567, -0.13807118745769906]
    H = zeros(T, 0, 3)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS3RodasTableau(T, T2)

A 3rd order Rosenbrock method.
"""
function ROS3RodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 3, 3)
    A[2, 1] = convert(T, 1.0)
    A[3, 1] = convert(T, 1.0)
    C = zeros(T, 3, 3)
    C[2, 1] = convert(T, -1.0156171083877703)
    C[3, 1] = convert(T, 4.07599564525377)
    C[3, 2] = convert(T, 9.20767942983308)
    c = T2[0.0, 0.435866521508459, 0.435866521508459]
    d = T[0.435866521508459, 0.24291996454816805, 2.185138002766406]
    b = T[1.0000000000000002, 6.1697947043828245, -0.42772256543218573]
    btilde = T[0.49999999999999983, -2.907955871680547, 0.22354069897811568]
    H = zeros(T, 0, 3)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS3PRRodasTableau(T, T2)

A 3rd order Rosenbrock method.
"""
function ROS3PRRodasTableau(T, T2)
    gamma = convert(T2, 0.788675134594813)
    A = zeros(T, 3, 3)
    A[2, 1] = convert(T, 3.0000000000000018)
    A[3, 1] = convert(T, 3.80384757729337)
    A[3, 2] = convert(T, 1.2679491924311226)
    C = zeros(T, 3, 3)
    C[2, 1] = convert(T, -3.80384757729337)
    C[3, 1] = convert(T, -5.673079295488928)
    C[3, 2] = convert(T, -1.7384634363911504)
    c = T2[0.0, 2.36602540378444, 1.0]
    d = T[0.788675134594813, -1.577350269189627, -0.577350269189621]
    b = T[4.5358983848622305, 1.2679491924311173, 1.0]
    btilde = T[-0.4598572229918937, -0.22992861149594657, 0.0]
    H = zeros(T, 0, 3)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    Scholz4_7RodasTableau(T, T2)

A 4th order Rosenbrock method with 7 stages by Scholz.
Reference: Scholz, S. (1989).
"""
function Scholz4_7RodasTableau(T, T2)
    gamma = convert(T2, 0.788675134594813)
    A = zeros(T, 3, 3)
    A[2, 1] = convert(T, 3.0000000000000018)
    A[3, 1] = convert(T, 4.120834875401151)
    A[3, 2] = convert(T, 1.2679491924311226)
    C = zeros(T, 3, 3)
    C[2, 1] = convert(T, -3.80384757729337)
    C[3, 1] = convert(T, -6.310062633644914)
    C[3, 2] = convert(T, -1.7746264440079669)
    c = T2[0.0, 2.36602540378444, 1.25]
    d = T[0.788675134594813, -1.577350269189627, -0.928571905524962]
    b = T[4.096775673951863, 0.953252779460065, 0.7833659121534696]
    btilde = T[0.3028225394953986, -0.06093909935296366, 0.36071618134309585]
    H = zeros(T, 0, 3)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS34PW1aRodasTableau(T, T2)

A Rosenbrock-W method of order (3)4.
Reference: Rang, J., & Angermann, L. (2005). New Rosenbrock-W methods of order 3 and 4 for stiff problems.
"""
function ROS34PW1aRodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 5.0905205106702045)
    A[4, 1] = convert(T, 4.005173696367865)
    A[4, 2] = convert(T, 0.19316470237944158)
    A[4, 3] = convert(T, 1.147140180139521)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -11.679081231228288)
    C[3, 1] = convert(T, -0.7100952636543062)
    C[3, 2] = convert(T, -0.04165460771675499)
    C[4, 1] = convert(T, -11.979557762226603)
    C[4, 2] = convert(T, -0.48054400523894975)
    C[4, 3] = convert(T, 1.4350493021655284)
    c = T2[0.0, 2.218787467653286, 0.0, 1.7837037931914073]
    d = T[0.435866521508459, -1.7829209461448272, 0.33333333333333337, -1.258070496147625]
    b = T[6.1538321465310215, -0.8364233759732359, -0.8614792120957679, 2.294280360279042]
    btilde = T[-5.429679341539398, -1.3273810331413745, 0.0, 0.0]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

function ROS34PW1bRodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 5.0905205106702045)
    A[3, 1] = convert(T, 5.0905205106702045)
    A[4, 1] = convert(T, 4.976281110107875)
    A[4, 2] = convert(T, 0.027726816471584953)
    A[4, 3] = convert(T, 0.22942803602790418)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -11.679081231228288)
    C[3, 1] = convert(T, -16.40573264673668)
    C[3, 2] = convert(T, -0.27726816471584953)
    C[4, 1] = convert(T, -8.38103960500476)
    C[4, 2] = convert(T, -0.8483284091993433)
    C[4, 3] = convert(T, 0.28700986043310556)
    c = T2[0.0, 2.218787467653286, 2.218787467653286, 1.553923375357884]
    d = T[0.435866521508459, -1.7829209461448272, -2.4654190049693425, -0.8055299979063697]
    b = T[5.2258276123309395, -0.5569711481541647, 0.35797946935364533, 1.7233739852106407]
    btilde = T[-5.168452127840395, -1.2635194260384186, 0.0, 0.0]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS34PW2RodasTableau(T, T2)

A Rosenbrock-W method of order (3)4.
Reference: Rang, J., & Angermann, L. (2005).
"""
function ROS34PW2RodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 1.4192173174557647)
    A[3, 2] = convert(T, -0.2592322116729697)
    A[4, 1] = convert(T, 4.18476048231916)
    A[4, 2] = convert(T, -0.28519201735549593)
    A[4, 3] = convert(T, 2.294280360279042)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -4.588560720558084)
    C[3, 1] = convert(T, -4.18476048231916)
    C[3, 2] = convert(T, 0.28519201735549593)
    C[4, 1] = convert(T, -6.368179200128359)
    C[4, 2] = convert(T, -6.795620944466837)
    C[4, 3] = convert(T, 2.8700986043310563)
    c = T2[0.0, 0.871733043016918, 0.7315799577888524, 1.0]
    d = T[0.435866521508459, -0.435866521508459, -0.4133333762338865, -5.551115123125783e-17]
    b = T[4.1847604823191595, -0.28519201735549565, 2.2942803602790414, 1.0]
    btilde = T[0.2777499476479681, -1.4032398951759992, 1.7726301276675507, 0.5]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS34PRwRodasTableau(T, T2)

A 3rd order Rosenbrock-W method.
"""
function ROS34PRwRodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 2.0)
    A[3, 1] = convert(T, 1.9166355646921893)
    A[3, 2] = convert(T, -0.7305046154473316)
    A[4, 1] = convert(T, 3.7075384385487764)
    A[4, 2] = convert(T, 1.984721005641544)
    A[4, 3] = convert(T, -0.7228174329072325)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -4.588560720558084)
    C[3, 1] = convert(T, -1.4496008611374558)
    C[3, 2] = convert(T, 2.6585485498967283)
    C[4, 1] = convert(T, -0.8142320398640468)
    C[4, 2] = convert(T, 2.1949369533270104)
    C[4, 3] = convert(T, -0.9042300763629808)
    c = T2[0.0, 0.871733043016918, 1.1537997822626886, 0.9999999999999999]
    d = T[0.435866521508459, -0.435866521508459, -0.34459816128502135, 5.551115123125783e-17]
    b = T[3.7075384385487764, 1.9847210056415439, -0.7228174329072324, 1.0]
    btilde = T[-0.08016142700721947, 0.15059517863671545, -0.29187352202361583, 0.26131506383377556]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS3PRLRodasTableau(T, T2)

A 3rd order low-storage Rosenbrock method.
"""
function ROS3PRLRodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 1.147140180139521)
    A[3, 1] = convert(T, 2.4630707730300534)
    A[3, 2] = convert(T, 1.147140180139521)
    A[4, 1] = convert(T, 2.4630707730300534)
    A[4, 2] = convert(T, 1.147140180139521)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -2.631861185781065)
    C[3, 1] = convert(T, -2.038451402734394)
    C[3, 2] = convert(T, 1.8551577240019121)
    C[4, 1] = convert(T, -1.8050630466729911)
    C[4, 2] = convert(T, 3.411439279441918)
    C[4, 3] = convert(T, -1.7057196397209593)
    c = T2[0.0, 0.5, 1.0, 1.0]
    d = T[0.435866521508459, -0.064133478491541, -0.0032561147686690495, 0.0]
    b = T[2.4630707730300534, 1.1471401801395211, 0.0, 1.0]
    btilde = T[0.14188781262114447, 0.9677841576948438, -0.06855321332716582, 0.26131506383377634]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    ROS3PRL2RodasTableau(T, T2)

A 3rd order low-storage Rosenbrock method.
"""
function ROS3PRL2RodasTableau(T, T2)
    gamma = convert(T2, 0.435866521508459)
    A = zeros(T, 4, 4)
    A[2, 1] = convert(T, 3.000000000000007)
    A[3, 1] = convert(T, 4.588560720558092)
    A[3, 2] = convert(T, 1.147140180139521)
    A[4, 1] = convert(T, 4.588560720558092)
    A[4, 2] = convert(T, 1.147140180139521)
    C = zeros(T, 4, 4)
    C[2, 1] = convert(T, -6.882841080837141)
    C[3, 1] = convert(T, -12.579179703104524)
    C[3, 2] = convert(T, -2.9475127180857186)
    C[4, 1] = convert(T, 3.556592450580385)
    C[4, 2] = convert(T, -0.4663124315128337)
    C[4, 3] = convert(T, 3.545260255335102)
    c = T2[0.0, 1.30759956452538, 1.0, 1.0]
    d = T[0.435866521508459, -0.871733043016921, -0.8339865967040412, -5.551115123125783e-17]
    b = T[4.5885607205580925, 1.1471401801395211, 0.0, 1.0]
    btilde = T[0.8808636262903982, 0.3041167493114213, 0.14248471842891264, 0.26131506383377634]
    H = zeros(T, 0, 4)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end

"""
    RosenbrockW6S4OSRodasTableau(T, T2)

A 6-stage 4th order Rosenbrock-W method.
"""
function RosenbrockW6S4OSRodasTableau(T, T2)
    gamma = convert(T2, 0.25)
    A = zeros(T, 6, 6)
    A[2, 1] = convert(T, 0.5812383407115008)
    A[3, 1] = convert(T, 0.903962441371467)
    A[3, 2] = convert(T, 1.861519155534501)
    A[4, 1] = convert(T, 2.076579719675)
    A[4, 2] = convert(T, 0.1884255381414796)
    A[4, 3] = convert(T, 1.870158967491032)
    A[5, 1] = convert(T, 4.435550638484312)
    A[5, 2] = convert(T, 5.457181798610189)
    A[5, 3] = convert(T, 4.61635078806893)
    A[5, 4] = convert(T, 3.118111952402361)
    A[6, 1] = convert(T, 10.79170169848326)
    A[6, 2] = convert(T, -10.05691522584131)
    A[6, 3] = convert(T, 14.99564485428419)
    A[6, 4] = convert(T, 5.274339954390943)
    A[6, 5] = convert(T, 1.42973087126119)
    C = zeros(T, 6, 6)
    C[2, 1] = convert(T, -2.661294105131369)
    C[3, 1] = convert(T, -3.128450202373838)
    C[4, 1] = convert(T, -6.920335474535658)
    C[4, 2] = convert(T, -1.202675288266817)
    C[4, 3] = convert(T, -9.73356181141362)
    C[5, 1] = convert(T, -28.09530629102695)
    C[5, 2] = convert(T, 20.37126295479377)
    C[5, 3] = convert(T, -41.04375275302869)
    C[5, 4] = convert(T, -19.66373175620895)
    C[6, 1] = convert(T, 9.7998186780974)
    C[6, 2] = convert(T, 11.93579288660318)
    C[6, 3] = convert(T, 3.673874929013201)
    C[6, 4] = convert(T, 14.8078285410955)
    C[6, 5] = convert(T, 0.831858399869068)
    c = T2[0.0, 0.1453095851778752, 0.3817422770256738, 0.6367813704374599, 0.7560744496323561, 0.927104723987567]
    d = T[0.25, 0.0836691184292894, 0.0544718623516351, -0.3402289722355864, 0.0337651588339529, -0.090307426761854]
    b = T[6.456217074653235, -4.853141317768053, 9.76531833406926, 2.081084177278723, 0.6603936866352417, 0.6]
    btilde = nothing
    H = zeros(T, 0, 6)
    return RodasTableau(A, C, gamma, c, d, H, b, btilde)
end


# Tsit5DA - 12-stage order 5(4) hybrid explicit/linear-implicit method for DAEs
# Reference: Steinebach (2025), arXiv:2511.21252
