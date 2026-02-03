"""
constructDormandPrince()

Constructs the tableau object for the Dormand-Prince Order 4/5 method.
"""
function constructDormandPrince(T::Type = Float64)
    A = [
        0 0 0 0 0 0 0
        1 // 5 0 0 0 0 0 0
        3 // 40 9 // 40 0 0 0 0 0
        44 // 45 -56 // 15 32 // 9 0 0 0 0
        19372 // 6561 -25360 // 2187 64448 // 6561 -212 // 729 0 0 0
        9017 // 3168 -355 // 33 46732 // 5247 49 // 176 -5103 // 18656 0 0
        35 // 384 0 500 // 1113 125 // 192 -2187 // 6784 11 // 84 0
    ]
    c = [0; 1 // 5; 3 // 10; 4 // 5; 8 // 9; 1; 1]
    α = [35 // 384; 0; 500 // 1113; 125 // 192; -2187 // 6784; 11 // 84; 0]
    αEEst = [
        5179 // 57600,
        0,
        7571 // 16695,
        393 // 640,
        -92097 // 339200,
        187 // 2100,
        1 // 40,
    ]
    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
            fsal = true, stability_size = 3.3066
        )
    )
end

# ============================================================================
# Tsit5 Interpolation Coefficients in Matrix Form
# ============================================================================

"""
    construct_tsit5_interp_matrix(::Type{T}) where {T <: Union{Float32, Float64}}

Constructs the interpolation coefficient matrix for Tsit5 method.
The matrix B_interp has dimensions (7, 5) where:
- Row i contains coefficients for stage i's interpolation polynomial
- Column j contains coefficients for Θ^(j-1) term
"""
function construct_tsit5_interp_matrix(::Type{T}) where {T <: Union{Float32, Float64}}
    r11 = convert(T, 1.0)
    r12 = convert(T, -2.763706197274826)
    r13 = convert(T, 2.9132554618219126)
    r14 = convert(T, -1.0530884977290216)

    r22 = convert(T, 0.13169999999999998)
    r23 = convert(T, -0.2234)
    r24 = convert(T, 0.1017)

    r32 = convert(T, 3.9302962368947516)
    r33 = convert(T, -5.941033872131505)
    r34 = convert(T, 2.490627285651253)

    r42 = convert(T, -12.411077166933676)
    r43 = convert(T, 30.33818863028232)
    r44 = convert(T, -16.548102889244902)

    r52 = convert(T, 37.50931341651104)
    r53 = convert(T, -88.1789048947664)
    r54 = convert(T, 47.37952196281928)

    r62 = convert(T, -27.896526289197286)
    r63 = convert(T, 65.09189467479366)
    r64 = convert(T, -34.87065786149661)

    r72 = convert(T, 1.5)
    r73 = convert(T, -4.0)
    r74 = convert(T, 2.5)

    B_interp = zeros(T, 7, 5)
    B_interp[1, :] = [0, r11, r12, r13, r14]
    B_interp[2, :] = [0, 0, r22, r23, r24]
    B_interp[3, :] = [0, 0, r32, r33, r34]
    B_interp[4, :] = [0, 0, r42, r43, r44]
    B_interp[5, :] = [0, 0, r52, r53, r54]
    B_interp[6, :] = [0, 0, r62, r63, r64]
    B_interp[7, :] = [0, 0, r72, r73, r74]

    return B_interp
end

construct_tsit5_interp_matrix() = construct_tsit5_interp_matrix(Float64)

"""
    construct_tsit5_interp_matrix(::Type{T}) where T

High-precision version for BigFloat and other arbitrary-precision types.
"""
function construct_tsit5_interp_matrix(::Type{T}) where {T}
    r11 = convert(T, big"0.999999999999999974283372471559910888475488471328")
    r12 = convert(T, big"-2.763706197274825911336735930481400260916070804192")
    r13 = convert(T, big"2.91325546182191274375068099306808")
    r14 = convert(T, -1.0530884977290216)

    r22 = convert(T, big"0.13169999999999999727")
    r23 = convert(T, big"-0.22339999999999999818")
    r24 = convert(T, 0.1017)

    r32 = convert(T, big"3.93029623689475152850687446709813398")
    r33 = convert(T, big"-5.94103387213150473470249202589458001")
    r34 = convert(T, big"2.490627285651252793")

    r42 = convert(T, big"-12.411077166933676983734381540685453484102414134010752")
    r43 = convert(T, big"30.3381886302823215981729903691836576")
    r44 = convert(T, big"-16.54810288924490272")

    r52 = convert(T, big"37.50931341651103919496903965334519631242339792120440212")
    r53 = convert(T, big"-88.1789048947664011014276693541209817")
    r54 = convert(T, big"47.37952196281928122")

    r62 = convert(T, big"-27.896526289197287805948263144598643896")
    r63 = convert(T, big"65.09189467479367152629021928716553658")
    r64 = convert(T, big"-34.87065786149660974")

    r72 = convert(T, 1.5)
    r73 = convert(T, -4.0)
    r74 = convert(T, 2.5)

    B_interp = zeros(T, 7, 5)
    B_interp[1, :] = [0, r11, r12, r13, r14]
    B_interp[2, :] = [0, 0, r22, r23, r24]
    B_interp[3, :] = [0, 0, r32, r33, r34]
    B_interp[4, :] = [0, 0, r42, r43, r44]
    B_interp[5, :] = [0, 0, r52, r53, r54]
    B_interp[6, :] = [0, 0, r62, r63, r64]
    B_interp[7, :] = [0, 0, r72, r73, r74]

    return B_interp
end

# ============================================================================
# Tsit5 in ExplicitRK tableau format
# ============================================================================

"""
    constructTsit5ExplicitRK(::Type{T}) where {T <: Union{Float32, Float64}}

Constructs the Tsitouras 5/4 method in ExplicitRK tableau format with compiled float coefficients.
"""
function constructTsit5ExplicitRK(::Type{T}) where {T <: Union{Float32, Float64}}
    A = [
        0.0        0.0         0.0         0.0          0.0          0.0       0.0
        0.161      0.0         0.0         0.0          0.0          0.0       0.0
        -0.008480655492356989  0.335480655492357  0.0  0.0  0.0  0.0  0.0
        2.8971530571054935   -6.359448489975075   4.3622954328695815  0.0  0.0  0.0  0.0
        5.325864828439257   -11.748883564062828   7.4955393428898365  -0.09249506636175525  0.0  0.0  0.0
        5.86145544294642    -12.92096931784711    8.159367898576159   -0.071584973281401   -0.028269050394068383  0.0  0.0
        0.09646076681806523   0.01  0.4798896504144996  1.379008574103742  -3.290069515436081  2.324710524099774  0.0
    ]

    c = [0.0, 0.161, 0.327, 0.9, 0.9800255409045097, 1.0, 1.0]

    α = [
        0.09468075576583945, 0.009183565540343254, 0.4877705284247616,
        1.234297566930479, -2.7077123499835256, 1.866628418170587,
        0.015151515151515152,
    ]

    btilde = [
        -0.00178001105222577714, -0.0008164344596567469, 0.007880878010261995,
        -0.1447110071732629, 0.5823571654525552, -0.45808210592918697,
        0.015151515151515152,
    ]

    αEEst = α .- btilde

    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    B_interp = construct_tsit5_interp_matrix(T)

    return DiffEqBase.ExplicitRKTableau(
        A, c, α, 5,
        αEEst = αEEst,
        adaptiveorder = 4,
        fsal = true,
        stability_size = 2.9,
        B_interp = B_interp
    )
end

constructTsit5ExplicitRK() = constructTsit5ExplicitRK(Float64)

"""
    constructTsit5ExplicitRK(::Type{T}) where T

High-precision version for BigFloat and other arbitrary-precision types.
"""
function constructTsit5ExplicitRK(::Type{T}) where {T}
    A = [
        0 0 0 0 0 0 0
        14 // 87 0 0 0 0 0 0
        -1 // 117 50 // 149 0 0 0 0 0
        310 // 107 -407 // 64 301 // 69 0 0 0 0
        474 // 89 -2479 // 211 817 // 109 -5 // 54 0 0 0
        381 // 65 -491 // 38 563 // 69 -19 // 265 -3 // 106 0 0
        8 // 83 1 // 100 107 // 223 131 // 95 -329 // 100 179 // 77 0
    ]

    c = [
        0; 161 // 1000; 327 // 1000; 9 // 10;
        big".9800255409045096857298102862870245954942137979563024768854764293221195950761080302604";
        1; 1
    ]

    α = [
        big".9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1",
        big".9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2",
        big".4877705284247615707855642599631228241516691959761363774365216240304071651579571959813",
        big"1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761",
        big"-2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702",
        big"1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089",
        1 // 66,
    ]

    btilde = [
        big"-1.780011052225771443378550607539534775944678804333659557637450799792588061629796e-03",
        big"-8.164344596567469032236360633546862401862537590159047610940604670770447527463931e-04",
        big"7.880878010261996010314727672526304238628733777103128603258129604952959142646516e-03",
        big"-1.44711007173262907537165147972635116720922712343167677619514233896760819649515e-01",
        big"5.823571654525552250199376106520421794260781239567387797673045438803694038950012e-01",
        big"-4.580821059291869466616365188325542974428047279788398179474684434732070620889539e-01",
        1 // 66,
    ]

    αEEst = α .- btilde

    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    B_interp = construct_tsit5_interp_matrix(T)

    return DiffEqBase.ExplicitRKTableau(
        A, c, α, 5,
        αEEst = αEEst,
        adaptiveorder = 4,
        fsal = true,
        stability_size = 2.9,
        B_interp = B_interp
    )
end

"""
ODE_DEFAULT_TABLEAU

Sets the default tableau for the ODE solver. Currently Dormand-Prince 4/5.
"""
const ODE_DEFAULT_TABLEAU = constructDormandPrince()

@doc raw"""
    ExplicitRK(; tableau = ODE_DEFAULT_TABLEAU)

A generic explicit Runge-Kutta method that allows you to define a custom tableau.
The default tableau is Dormand-Prince 4/5. This solver is primarily for research
purposes or when you need a specific tableau not already implemented.

# Parameters
- `tableau`: A `DiffEqBase.ExplicitRKTableau` object defining the Runge-Kutta tableau.

For most applications, prefer the named methods like `DP5()`, `Tsit5()`, etc.
"""
struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
    tableau::TabType
end
ExplicitRK(; tableau = ODE_DEFAULT_TABLEAU) = ExplicitRK(tableau)

@truncate_stacktrace ExplicitRK
