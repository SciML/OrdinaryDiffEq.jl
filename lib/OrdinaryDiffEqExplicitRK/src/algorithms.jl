"""
constructDormandPrince()

Constructs the tableau object for the Dormand-Prince Order 4/5 method.
"""
function constructDormandPrince(::Type{T} = Float64) where T
    A = [0 0 0 0 0 0 0
         1//5 0 0 0 0 0 0
         3//40 9//40 0 0 0 0 0
         44//45 -56//15 32//9 0 0 0 0
         19372//6561 -25360//2187 64448//6561 -212//729 0 0 0
         9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
         35//384 0 500//1113 125//192 -2187//6784 11//84 0]
    c = [0; 1 // 5; 3 // 10; 4 // 5; 8 // 9; 1; 1]
    α = [35 // 384; 0; 500 // 1113; 125 // 192; -2187 // 6784; 11 // 84; 0]
    αEEst = [
        5179 // 57600,
        0,
        7571 // 16695,
        393 // 640,
        -92097 // 339200,
        187 // 2100,
        1 // 40
    ]
    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
        fsal = true, stability_size = 3.3066))
end

# ============================================================================
# Tsit5 Interpolation Coefficients in Matrix Form
# ============================================================================

"""
    construct_tsit5_interp_matrix(T::Type = Float64)

Constructs the interpolation coefficient matrix for Tsit5 method.
This converts the polynomial coefficients from the original Tsit5 implementation
into a matrix format for generic interpolation.

The matrix B_interp has dimensions (7, 5) where:
- Row i contains coefficients for stage i's interpolation polynomial
- Column j contains coefficients for Θ^(j-1) term

Each polynomial bᵢ(Θ) is defined as:
    bᵢ(Θ) = bᵢ₀ + bᵢ₁*Θ + bᵢ₂*Θ² + bᵢ₃*Θ³ + bᵢ₄*Θ⁴

For Tsit5, the original formulation was:
    b₁(Θ) = Θ * (r11 + r12*Θ + r13*Θ² + r14*Θ³)
          = 0 + r11*Θ + r12*Θ² + r13*Θ³ + r14*Θ⁴

    b₂(Θ) = Θ² * (r22 + r23*Θ + r24*Θ²)
          = 0 + 0*Θ + r22*Θ² + r23*Θ³ + r24*Θ⁴

    ... and so on for all 7 stages
"""
function construct_tsit5_interp_matrix(::Type{T} = Float64) where T
    # Original Tsit5 interpolation coefficients
    # From OrdinaryDiffEqTsit5/src/tsit_tableaus.jl

    # Stage 1: b₁(Θ) = Θ * (r11 + r12*Θ + r13*Θ² + r14*Θ³)
    r11 = convert(T, 1.0)
    r12 = convert(T, -2.763706197274826)
    r13 = convert(T, 2.9132554618219126)
    r14 = convert(T, -1.0530884977290216)

    # Stage 2: b₂(Θ) = Θ² * (r22 + r23*Θ + r24*Θ²)
    r22 = convert(T, 0.13169999999999998)
    r23 = convert(T, -0.2234)
    r24 = convert(T, 0.1017)

    # Stage 3: b₃(Θ) = Θ² * (r32 + r33*Θ + r34*Θ²)
    r32 = convert(T, 3.9302962368947516)
    r33 = convert(T, -5.941033872131505)
    r34 = convert(T, 2.490627285651253)

    # Stage 4: b₄(Θ) = Θ² * (r42 + r43*Θ + r44*Θ²)
    r42 = convert(T, -12.411077166933676)
    r43 = convert(T, 30.33818863028232)
    r44 = convert(T, -16.548102889244902)

    # Stage 5: b₅(Θ) = Θ² * (r52 + r53*Θ + r54*Θ²)
    r52 = convert(T, 37.50931341651104)
    r53 = convert(T, -88.1789048947664)
    r54 = convert(T, 47.37952196281928)

    # Stage 6: b₆(Θ) = Θ² * (r62 + r63*Θ + r64*Θ²)
    r62 = convert(T, -27.896526289197286)
    r63 = convert(T, 65.09189467479366)
    r64 = convert(T, -34.87065786149661)

    # Stage 7: b₇(Θ) = Θ² * (r72 + r73*Θ + r74*Θ²)
    r72 = convert(T, 1.5)
    r73 = convert(T, -4.0)
    r74 = convert(T, 2.5)

    # Construct the interpolation matrix
    # B_interp[i, j] = coefficient of Θ^(j-1) in bᵢ(Θ)
    B_interp = zeros(T, 7, 5)

    # Stage 1: bᵢ(Θ) = 0 + r11*Θ + r12*Θ² + r13*Θ³ + r14*Θ⁴
    B_interp[1, :] = [0, r11, r12, r13, r14]

    # Stages 2-7: bᵢ(Θ) = 0 + 0*Θ + ri2*Θ² + ri3*Θ³ + ri4*Θ⁴
    B_interp[2, :] = [0, 0, r22, r23, r24]
    B_interp[3, :] = [0, 0, r32, r33, r34]
    B_interp[4, :] = [0, 0, r42, r43, r44]
    B_interp[5, :] = [0, 0, r52, r53, r54]
    B_interp[6, :] = [0, 0, r62, r63, r64]
    B_interp[7, :] = [0, 0, r72, r73, r74]

    return B_interp
end

"""
    construct_tsit5_interp_matrix_highprecision(T::Type)

High-precision version for BigFloat and other arbitrary-precision types.
We have not tested this
"""
function construct_tsit5_interp_matrix_highprecision(::Type{T}) where T
    # Stage 1
    r11 = convert(T, big"0.999999999999999974283372471559910888475488471328")
    r12 = convert(T, big"-2.763706197274825911336735930481400260916070804192")
    r13 = convert(T, big"2.91325546182191274375068099306808")
    r14 = convert(T, -1.0530884977290216)

    # Stage 2
    r22 = convert(T, big"0.13169999999999999727")
    r23 = convert(T, big"-0.22339999999999999818")
    r24 = convert(T, 0.1017)

    # Stage 3
    r32 = convert(T, big"3.93029623689475152850687446709813398")
    r33 = convert(T, big"-5.94103387213150473470249202589458001")
    r34 = convert(T, big"2.490627285651252793")

    # Stage 4
    r42 = convert(T, big"-12.411077166933676983734381540685453484102414134010752")
    r43 = convert(T, big"30.3381886302823215981729903691836576")
    r44 = convert(T, big"-16.54810288924490272")

    # Stage 5
    r52 = convert(T, big"37.50931341651103919496903965334519631242339792120440212")
    r53 = convert(T, big"-88.1789048947664011014276693541209817")
    r54 = convert(T, big"47.37952196281928122")

    # Stage 6
    r62 = convert(T, big"-27.896526289197287805948263144598643896")
    r63 = convert(T, big"65.09189467479367152629021928716553658")
    r64 = convert(T, big"-34.87065786149660974")

    # Stage 7
    r72 = convert(T, 1.5)
    r73 = convert(T, -4.0)
    r74 = convert(T, 2.5)

    # Construct matrix
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

"""
    construct_tsit5_interp_matrix_auto(T::Type)

Automatically selects appropriate precision based on type.
"""
construct_tsit5_interp_matrix_auto(::Type{T}) where {T <: Union{Float32, Float64}} = construct_tsit5_interp_matrix(T)
construct_tsit5_interp_matrix_auto(::Type{T}) where T = construct_tsit5_interp_matrix_highprecision(T)

# Convert Tsit5 tableau to ExplicitRK format

"""
    constructTsit5ExplicitRK(T::Type = Float64)

Constructs the Tsitouras 5/4 method in ExplicitRK tableau format.
This allows using Tsit5 with the generic ExplicitRK solver.

Tsit5 is a 7-stage, 5th-order method with 4th-order embedded error estimate.
"""
function constructTsit5ExplicitRK(::Type{T} = Float64) where T
    # Build the A matrix (Butcher tableau coefficients)
    # 7 stages, lower triangular (explicit method)
    A=[0 0 0 0 0 0 0
       14//87 0 0 0 0 0 0
       -1//117 50//149 0 0 0 0 0
       310//107 -407//64 301//69 0 0 0 0
       474//89 -2479//211 817//109 -5//54 0 0 0
       381//65 -491//38 563//69 -19//265 -3//106 0 0
       8//83 1//100 107//223 131//95 -329//100 179//77 0]
    # A = Float8.(A)

    # Time nodes (c vector)
    c = [0; 161//1000; 327//1000; 9//10;
         big".9800255409045096857298102862870245954942137979563024768854764293221195950761080302604";
         1; 1]


    # Solution weights (b vector) - 5th order
    α = [
        big".9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1",
        big".9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2",
        big".4877705284247615707855642599631228241516691959761363774365216240304071651579571959813",
        big"1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761",
        big"-2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702",
        big"1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089",
        1//66  # = 0.015151515151515152
    ]
    # Error estimate weights (b̂ vector) - 4th order
    # Note: In Tsit5, btilde = b - b̂, so b̂ = b - btilde
    btilde = [
        big"-1.780011052225771443378550607539534775944678804333659557637450799792588061629796e-03",
        big"-8.164344596567469032236360633546862401862537590159047610940604670770447527463931e-04",
        big"7.880878010261996010314727672526304238628733777103128603258129604952959142646516e-03",
        big"-1.44711007173262907537165147972635116720922712343167677619514233896760819649515e-01",
        big"5.823571654525552250199376106520421794260781239567387797673045438803694038950012e-01",
        big"-4.580821059291869466616365188325542974428047279788398179474684434732070620889539e-01",
        1//66
    ]

    # Calculate b̂ = b - btilde for the embedded 4th-order method
    αEEst = α .- btilde

    # Convert to requested type
    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    B_interp=construct_tsit5_interp_matrix_highprecision(T)

    return DiffEqBase.ExplicitRKTableau(A, c, α, 5,
        αEEst = αEEst,
        adaptiveorder = 4,
        fsal = true,
        stability_size = 2.9,  # Approximate stability region size
        B_interp = B_interp)
end

"""
    constructTsit5ExplicitRKSimple(T::Type = Float64)

Simplified version using rational and decimal approximations.
Faster to construct but slightly less accurate than the full precision version.
"""
function constructTsit5ExplicitRKSimple(::Type{T} = Float64) where T
    #Tested a few more variants and leaving them commented out here for future reference
    # Build the A matrix with simpler rationals/decimals
    # A = [0          0          0          0          0          0       0
    #      0.161      0          0          0          0          0       0
    #      -0.00848   0.3355     0          0          0          0       0
    #      2.8972     -6.3594    4.3623     0          0          0       0
    #      5.3259     -11.7489   7.4955     -0.0925    0          0       0
    #      5.8615     -12.9210   8.1594     -0.0716    -0.0283    0       0
    #      0.09646    0.01       0.4799     1.3790     -3.2901    2.3247  0]
#   A = [0         0         0         0         0         0      0
#          161//1000 0         0         0         0         0      0
#          -8480655492356989//1000000000000000000 335480655492357//1000000000000000 0 0 0 0 0
#          2897153057105493//1000000000000000 -6359448489975075//1000000000000000 4362295432869582//1000000000000000 0 0 0 0
#          5325864828439257//1000000000000000 -11748883564062828//10000000000000000 7495539342889836//1000000000000000 -92495066361755//1000000000000000 0 0 0
#          5861455442946420//1000000000000000 -12920969317847109//1000000000000000 8159367898576159//1000000000000000 -71584973281401//1000000000000000 -28269050394068//1000000000000000 0 0
#          96460766818065//1000000000000000 1//100 479889650414500//1000000000000000 1379008574103742//1000000000000000 -3290069515436081//1000000000000000 2324710524099774//1000000000000000 0]

#    A = Float64.(A)
A=[0 0 0 0 0 0 0
   14//87 0 0 0 0 0 0
   -1//117 50//149 0 0 0 0 0
   310//107 -407//64 301//69 0 0 0 0
   474//89 -2479//211 817//109 -5//54 0 0 0
   381//65 -491//38 563//69 -19//265 -3//106 0 0
   8//83 1//100 107//223 131//95 -329//100 179//77 0]
# A=[0.0 0.0 0.0 0.0 0.0 0.0 0.0
#    0.161 0.0 0.0 0.0 0.0 0.0 0.0
#    -0.008484 0.3354 0.0 0.0 0.0 0.0 0.0
#    2.896 -6.36 4.363 0.0 0.0 0.0 0.0
#    5.324 -1.175 7.496 -0.09247 0.0 0.0 0.0
#    5.863 -12.92 8.16 -0.0716 -0.02827 0.0 0.0
#    0.09644 0.01 0.48 1.379 -3.291 2.324 0.0]

    # Time nodes
    c = [0, 0.161, 0.327, 0.9, 0.9800255409045097, 1.0, 1.0]


    # Solution weights (5th order)
    α = [0.09468075576583945, 0.009183565540343254, 0.4877705284247616,
         1.234297566930479, -2.7077123499835256, 1.866628418170587,
         0.015151515151515152]

    # Error estimate - computed from btilde
    btilde = [-0.00178001105222577714, -0.0008164344596567469, 0.007880878010261995,
              -0.1447110071732629, 0.5823571654525552, -0.45808210592918697,
              0.015151515151515152]

    αEEst = α .- btilde

    # Convert to requested type
    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    B_interp=construct_tsit5_interp_matrix(T)

    return DiffEqBase.ExplicitRKTableau(A, c, α, 5,
        αEEst = αEEst,
        adaptiveorder = 4,
        fsal = true,
        stability_size = 2.9,
        B_interp = B_interp)
end

# Example usage:
# tableau = constructTsit5ExplicitRK()
# solve(prob, ExplicitRK(tableau = tableau))

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
