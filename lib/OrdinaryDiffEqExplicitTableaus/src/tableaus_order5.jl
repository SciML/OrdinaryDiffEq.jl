"""
    DormandPrince()

Constructs the tableau object for the Dormand-Prince Order 4/5 method (DP5).
"""
function DormandPrince(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
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
    c = map(T_time, c)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
            fsal = true, stability_size = 3.3066
        )
    )
end

"""
Runge-Kutta-Fehlberg Order 4/5 method.
"""
function RKF5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0 0 0
        1 // 4 0 0 0 0 0
        3 // 32 9 // 32 0 0 0 0
        1932 // 2197 -7200 // 2197 7296 // 2197 0 0 0
        439 // 216 -8 3680 // 513 -845 // 4104 0 0
        -8 // 27 2 -3544 // 2565 1859 // 4104 -11 // 40 0
    ]
    c = [0; 1 // 4; 3 // 8; 12 // 13; 1; 1 // 2]
    α = [16 // 135; 0; 6656 // 12825; 28561 // 56430; -9 // 50; 2 // 55]
    αEEst = [25 // 216; 0; 1408 // 2565; 2197 // 4104; -1 // 5; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
            stability_size = 3.6777066213218945
        )
    )
end

"""

Runge's First Order 5 method

"""
function RungeFirst5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 6, 6)
    c = zeros(T_time, 6)
    α = zeros(T, 6)
    c[2] = 1 // 5
    c[3] = 2 // 5
    c[4] = 1
    c[5] = 3 // 5
    c[6] = 4 // 5
    A[2, 1] = 1 // 5
    A[3, 1] = 0
    A[3, 2] = 2 // 5
    A[4, 1] = 9 // 4
    A[4, 2] = -5
    A[4, 3] = 15 // 4
    A[5, 1] = -63 // 100
    A[5, 2] = 9 // 5
    A[5, 3] = -13 // 20
    A[5, 4] = 2 // 25
    A[6, 1] = -6 // 25
    A[6, 2] = 4 // 5
    A[6, 3] = 2 // 15
    A[6, 4] = 8 // 75
    A[6, 5] = 0
    α[1] = 17 // 144
    α[2] = 0
    α[3] = 25 // 36
    α[4] = 1 // 72
    α[5] = -25 // 72
    α[6] = 25 // 48
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, stability_size = 3.2170478666401054))
end

"""

Cassity's Order 5 method
"""
function Cassity5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 6, 6)
    c = zeros(T_time, 6)
    α = zeros(T, 6)
    c[2] = 1 // 7
    c[3] = 5 // 14
    c[4] = 9 // 14
    c[5] = 6 // 7
    c[6] = 1
    A[2, 1] = 1 // 7
    A[3, 1] = -367 // 4088
    A[3, 2] = 261 // 584
    A[4, 1] = 41991 // 2044
    A[4, 2] = -2493 // 73
    A[4, 3] = 57 // 4
    A[5, 1] = -108413 // 196224
    A[5, 2] = 58865 // 65408
    A[5, 3] = 5 // 16
    A[5, 4] = 265 // 1344
    A[6, 1] = -204419 // 58984
    A[6, 2] = 143829 // 58984
    A[6, 3] = 171 // 202
    A[6, 4] = 2205 // 404
    A[6, 5] = -432 // 101
    α[1] = 1 // 9
    α[2] = 7 // 2700
    α[3] = 413 // 810
    α[4] = 7 // 450
    α[5] = 28 // 75
    α[6] = -101 // 8100
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, stability_size = 2.1686021275264866))
end

"""

Lawson's 5th order scheme

An Order Five Runge Kutta Process with Extended Region of Stability, J. Douglas Lawson,
 Siam Journal on Numerical Analysis, Vol. 3, No. 4, (Dec., 1966) pages 593-597

"""
function Lawson5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 6, 6)
    c = zeros(T_time, 6)
    α = zeros(T, 6)
    c[2] = 1 // 12
    c[3] = 1 // 4
    c[4] = 1 // 2
    c[5] = 3 // 4
    c[6] = 1
    A[2, 1] = 1 // 12
    A[3, 1] = -1 // 8
    A[3, 2] = 3 // 8
    A[4, 1] = 3 // 5
    A[4, 2] = -9 // 10
    A[4, 3] = 4 // 5
    A[5, 1] = 39 // 80
    A[5, 2] = -9 // 20
    A[5, 3] = 3 // 20
    A[5, 4] = 9 // 16
    A[6, 1] = -59 // 35
    A[6, 2] = 66 // 35
    A[6, 3] = 48 // 35
    A[6, 4] = -12 // 7
    A[6, 5] = 8 // 7
    α[1] = 7 // 90
    α[2] = 0
    α[3] = 16 // 45
    α[4] = 2 // 15
    α[5] = 16 // 45
    α[6] = 7 // 90
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, stability_size = 3.7343596072347225))
end

"""

Luther and Konen's First Order 5
Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,
 Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

"""
function LutherKonen5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 6, 6)
    c = zeros(T_time, 6)
    α = zeros(T, 6)
    c[2] = 1 / 2
    c[3] = 1 / 2 - 1 / 10 * 5^(1 / 2)
    c[4] = 1 / 2
    c[5] = 1 / 2 + 1 / 10 * 5^(1 / 2)
    c[6] = 1
    A[2, 1] = 1 / 2
    A[3, 1] = 1 / 5
    A[3, 2] = 3 / 10 - 1 / 10 * 5^(1 / 2)
    A[4, 1] = 1 / 4
    A[4, 2] = 1 / 4
    A[4, 3] = 0
    A[5, 1] = 1 / 20 - 1 / 20 * 5^(1 / 2)
    A[5, 2] = -1 / 5
    A[5, 3] = 1 / 4 + 3 / 20 * 5^(1 / 2)
    A[5, 4] = 2 / 5
    A[6, 1] = 1 / 4 * 5^(1 / 2) - 1 / 4
    A[6, 2] = 1 / 2 * 5^(1 / 2) - 1 / 2
    A[6, 3] = 5 / 4 - 1 / 4 * 5^(1 / 2)
    A[6, 4] = -2
    A[6, 5] = 5 / 2 - 1 / 2 * 5^(1 / 2)
    α[1] = 1 / 12
    α[2] = 0
    α[3] = 5 / 12
    α[4] = 0
    α[5] = 5 / 12
    α[6] = 1 / 12
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, stability_size = 3.217047866640106))
end

"""

Luther and Konen's Second Order 5
Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,
 Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

"""
function LutherKonen52(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 6, 6)
    c = zeros(T_time, 6)
    α = zeros(T, 6)
    c[2] = 2 / 5
    c[3] = 1 / 2
    c[4] = 1
    c[5] = 1 / 2 - 1 / 10 * 15^(1 / 2)
    c[6] = 1 / 2 + 1 / 10 * 15^(1 / 2)
    A[2, 1] = 2 / 5
    A[3, 1] = 3 / 16
    A[3, 2] = 5 / 16
    A[4, 1] = 1 / 4
    A[4, 2] = -5 / 4
    A[4, 3] = 2
    A[5, 1] = 3 / 20 - 1 / 100 * 15^(1 / 2)
    A[5, 2] = -1 / 4
    A[5, 3] = 3 / 5 - 2 / 25 * 15^(1 / 2)
    A[5, 4] = -1 / 100 * 15^(1 / 2)
    A[6, 1] = -3 / 20 - 1 / 20 * 15^(1 / 2)
    A[6, 2] = -1 / 4
    A[6, 3] = 3 / 5
    A[6, 4] = 3 / 10 - 1 / 20 * 15^(1 / 2)
    A[6, 5] = 1 / 5 * 15^(1 / 2)
    α[1] = 0
    α[2] = 0
    α[3] = 4 / 9
    α[4] = 0
    α[5] = 5 / 18
    α[6] = 5 / 18
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, stability_size = 2.6515956444339794))
end

"""

Luther and Konen's Third Order 5
Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,
 Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

"""
function LutherKonen53(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 6, 6)
    c = zeros(T_time, 6)
    α = zeros(T, 6)
    c[2] = 3 / 25
    c[3] = 5 / 18
    c[4] = 45 / 89
    c[5] = 3 / 5 - 1 / 10 * 6^(1 / 2)
    c[6] = 3 / 5 + 1 / 10 * 6^(1 / 2)
    A[2, 1] = 3 / 25
    A[3, 1] = -85 / 1944
    A[3, 2] = 625 / 1944
    A[4, 1] = 610425 / 1409938
    A[4, 2] = -961875 / 1409938
    A[4, 3] = 532170 / 704969
    A[5, 1] = 7411 / 37500 - 673 / 18750 * 6^(1 / 2)
    A[5, 2] = 0
    A[5, 3] = 27621 / 228125 * 6^(1 / 2) - 6561 / 228125
    A[5, 4] = 1180229 / 2737500 - 126736 / 684375 * 6^(1 / 2)
    A[6, 1] = -5351 / 62500 - 7087 / 281250 * 6^(1 / 2)
    A[6, 2] = 0
    A[6, 3] = 2736423 / 1140625 + 786753 / 1140625 * 6^(1 / 2)
    A[6, 4] = 73736589 / 86687500 + 101816534 / 195046875 * 6^(1 / 2)
    A[6, 5] = -30448 / 11875 - 12903 / 11875 * 6^(1 / 2)
    α[1] = 1 / 9
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 4 / 9 + 1 / 36 * 6^(1 / 2)
    α[6] = 4 / 9 - 1 / 36 * 6^(1 / 2)
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, stability_size = 3.679935798458283))
end

"""

 S.N. Papakostas and G. PapaGeorgiou higher error more stable

 A Family of Fifth-order Runge-Kutta Pairs, by S.N. Papakostas and G. PapaGeorgiou,
 Mathematics of Computation,Volume 65, Number 215, July 1996, Pages 1165-1181.

"""
function PapakostasPapaGeorgiou5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 7, 7)
    c = zeros(T_time, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)
    c[2] = 64 // 315
    c[3] = 115 // 381
    c[4] = 762 // 935
    c[5] = 25 // 28
    c[6] = 1
    c[7] = 1
    A[2, 1] = 64 // 315
    A[3, 1] = 480815 // 6193536
    A[3, 2] = 462875 // 2064512
    A[4, 1] = 344904825219069 // 345923838700000
    A[4, 2] = -2360077908867 // 601606676000
    A[4, 3] = 40439332863108 // 10810119959375
    A[5, 1] = 12078745127989699 // 5009699157786624
    A[5, 2] = -791781731775 // 81669668864
    A[5, 3] = 39297175833216951 // 4694461413969824
    A[5, 4] = -10508413393960625 // 54097233987826176
    A[6, 1] = 2251421737440863 // 828701767536000
    A[6, 2] = -39895842357 // 3782730880
    A[6, 3] = 34564628685305112534 // 3916944972468643375
    A[6, 4] = 12051135495733565 // 36492943960723917
    A[6, 5] = -26808346215168 // 82592376030125
    A[7, 1] = 2405713 // 26289000
    A[7, 2] = 0
    A[7, 3] = 63896466577779 // 141024193000600
    A[7, 4] = 454128848141375 // 589615117674696
    A[7, 5] = -1359311744 // 2892576375
    A[7, 6] = 256979 // 1656648
    α[1] = 2405713 // 26289000
    α[2] = 0
    α[3] = 63896466577779 // 141024193000600
    α[4] = 454128848141375 // 589615117674696
    α[5] = -1359311744 // 2892576375
    α[6] = 256979 // 1656648
    αEEst[1] = 1818563883019 // 20194131951000
    αEEst[2] = 0
    αEEst[3] = 5513862498202899713 // 12036555896794210600
    αEEst[4] = 324806515311046773125 // 452918159177876804664
    αEEst[5] = -126112324722496 // 317422653663375
    αEEst[6] = 137695258717 // 1272569071032
    αEEst[7] = 1 // 42
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4, fsal = true,
            stability_size = 3.306567892634947
        )
    )
end

"""

 S.N. Papakostas and G. PapaGeorgiou less stable lower error
 Strictly better than DP5

 A Family of Fifth-order Runge-Kutta Pairs, by S.N. Papakostas and G. PapaGeorgiou,
 Mathematics of Computation,Volume 65, Number 215, July 1996, Pages 1165-1181.

"""
function PapakostasPapaGeorgiou52(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 7, 7)
    c = zeros(T_time, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)

    c[2] = 35 // 159
    c[3] = 42 // 131
    c[4] = 131 // 143
    c[5] = 21 // 22
    c[6] = 1
    c[7] = 1
    A[2, 1] = 35 // 159
    A[3, 1] = 7476 // 85805
    A[3, 2] = 20034 // 85805
    A[4, 1] = 2438549411 // 1983961980
    A[4, 2] = -3707256508 // 716430715
    A[4, 3] = 25077455105 // 5158301148
    A[5, 1] = 105337889067 // 64388030080
    A[5, 2] = -1698584121 // 245755840
    A[5, 3] = 6869523776931 // 1096562558080
    A[5, 4] = -927215289 // 26981535520
    A[6, 1] = 67512025387 // 32454592380
    A[6, 2] = -20051384 // 2293935
    A[6, 3] = 10587214001321 // 1373901639516
    A[6, 4] = 731293420 // 8209319229
    A[6, 5] = -144610048 // 1077690663
    A[7, 1] = 669707 // 6932520
    A[7, 2] = 0
    A[7, 3] = 2215522905683 // 4570867891800
    A[7, 4] = 349043981 // 116904400
    A[7, 5] = -2234144 // 575505
    A[7, 6] = 9363 // 7120
    α[1] = 669707 // 6932520
    α[2] = 0
    α[3] = 2215522905683 // 4570867891800
    α[4] = 349043981 // 116904400
    α[5] = -2234144 // 575505
    α[6] = 9363 // 7120
    αEEst[1] = 2243660497 // 23535905400
    αEEst[2] = 0
    αEEst[3] = 7589131232781673 // 15518096492661000
    αEEst[4] = 1104461697911 // 396890438000
    αEEst[5] = -6925033984 // 1953839475
    αEEst[6] = 3529851 // 3021550
    αEEst[7] = 1 // 112

    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4, fsal = true,
            stability_size = 3.4272251630453394
        )
    )
end

"""
Runge–Kutta pairs of orders 5(4) using the minimal set of simplifying assumptions,
 by Ch. Tsitouras, TEI of Chalkis, Dept. of Applied Sciences, GR34400, Psahna, Greece.
"""
function Tsitouras5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 7, 7)
    c = zeros(T_time, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)
    c[2] = convert(T, 161 // 1000)
    c[3] = convert(T, 327 // 1000)
    c[4] = convert(T, 9 // 10)
    c[5] = convert(
        T,
        big".9800255409045096857298102862870245954942137979563024768854764293221195950761080302604"
    )
    c[6] = convert(T, 1)
    c[7] = convert(T, 1)
    A[2, 1] = convert(T, 161 // 1000)
    A[3, 1] = convert(
        T,
        big"-.8480655492356988544426874250230774675121177393430391537369234245294192976164141156943e-2"
    )
    A[3, 2] = convert(
        T,
        big".3354806554923569885444268742502307746751211773934303915373692342452941929761641411569"
    )
    A[4, 1] = convert(
        T,
        big"2.897153057105493432130432594192938764924887287701866490314866693455023795137503079289"
    )
    A[4, 2] = convert(
        T,
        big"-6.359448489975074843148159912383825625952700647415626703305928850207288721235210244366"
    )
    A[4, 3] = convert(
        T,
        big"4.362295432869581411017727318190886861027813359713760212991062156752264926097707165077"
    )
    A[5, 1] = convert(
        T,
        big"5.325864828439256604428877920840511317836476253097040101202360397727981648835607691791"
    )
    A[5, 2] = convert(
        T,
        big"-11.74888356406282787774717033978577296188744178259862899288666928009020615663593781589"
    )
    A[5, 3] = convert(
        T,
        big"7.495539342889836208304604784564358155658679161518186721010132816213648793440552049753"
    )
    A[5, 4] = convert(
        T,
        big"-.9249506636175524925650207933207191611349983406029535244034750452930469056411389539635e-1"
    )
    A[6, 1] = convert(
        T,
        big"5.861455442946420028659251486982647890394337666164814434818157239052507339770711679748"
    )
    A[6, 2] = convert(
        T,
        big"-12.92096931784710929170611868178335939541780751955743459166312250439928519268343184452"
    )
    A[6, 3] = convert(
        T,
        big"8.159367898576158643180400794539253485181918321135053305748355423955009222648673734986"
    )
    A[6, 4] = convert(
        T,
        big"-.7158497328140099722453054252582973869127213147363544882721139659546372402303777878835e-1"
    )
    A[6, 5] = convert(
        T,
        big"-.2826905039406838290900305721271224146717633626879770007617876201276764571291579142206e-1"
    )
    A[7, 1] = convert(
        T,
        big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1"
    )
    A[7, 2] = convert(T, 1 // 100)
    A[7, 3] = convert(
        T,
        big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509"
    )
    A[7, 4] = convert(
        T,
        big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331"
    )
    A[7, 5] = convert(
        T,
        big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677"
    )
    A[7, 6] = convert(
        T,
        big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841"
    )
    α[1] = convert(
        T,
        big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1"
    )
    α[2] = convert(T, 1 // 100)
    α[3] = convert(
        T,
        big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509"
    )
    α[4] = convert(
        T,
        big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331"
    )
    α[5] = convert(
        T,
        big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677"
    )
    α[6] = convert(
        T,
        big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841"
    )
    αEEst[1] = convert(
        T,
        big".9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1"
    )
    αEEst[2] = convert(
        T,
        big".9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2"
    )
    αEEst[3] = convert(
        T,
        big".4877705284247615707855642599631228241516691959761363774365216240304071651579571959813"
    )
    αEEst[4] = convert(
        T,
        big"1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761"
    )
    αEEst[5] = convert(
        T,
        big"-2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702"
    )
    αEEst[6] = convert(
        T,
        big"1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089"
    )
    αEEst[7] = convert(T, 1 // 66)
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4, fsal = true,
            stability_size = 3.5068469938049547
        )
    )
end

"""

An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28
"""
function BogakiShampine5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 8, 8)
    c = zeros(T_time, 8)
    α = zeros(T, 8)
    αEEst = zeros(T, 8)
    αEEst2 = zeros(T, 8)
    c[2] = 1 // 6
    c[3] = 2 // 9
    c[4] = 3 // 7
    c[5] = 2 // 3
    c[6] = 3 // 4
    c[7] = 1
    c[8] = 1
    A[2, 1] = 1 // 6
    A[3, 1] = 2 // 27
    A[3, 2] = 4 // 27
    A[4, 1] = 183 // 1372
    A[4, 2] = -162 // 343
    A[4, 3] = 1053 // 1372
    A[5, 1] = 68 // 297
    A[5, 2] = -4 // 11
    A[5, 3] = 42 // 143
    A[5, 4] = 1960 // 3861
    A[6, 1] = 597 // 22528
    A[6, 2] = 81 // 352
    A[6, 3] = 63099 // 585728
    A[6, 4] = 58653 // 366080
    A[6, 5] = 4617 // 20480
    A[7, 1] = 174197 // 959244
    A[7, 2] = -30942 // 79937
    A[7, 3] = 8152137 // 19744439
    A[7, 4] = 666106 // 1039181
    A[7, 5] = -29421 // 29068
    A[7, 6] = 482048 // 414219
    A[8, 1] = 587 // 8064
    A[8, 2] = 0
    A[8, 3] = 4440339 // 15491840
    A[8, 4] = 24353 // 124800
    A[8, 5] = 387 // 44800
    A[8, 6] = 2152 // 5985
    A[8, 7] = 7267 // 94080
    α[1] = 587 // 8064
    α[2] = 0
    α[3] = 4440339 // 15491840
    α[4] = 24353 // 124800
    α[5] = 387 // 44800
    α[6] = 2152 // 5985
    α[7] = 7267 // 94080
    α[8] = 0
    αEEst[1] = 6059 // 80640
    αEEst[2] = 0
    αEEst[3] = 8559189 // 30983680
    αEEst[4] = 26411 // 124800
    αEEst[5] = -927 // 89600
    αEEst[6] = 443 // 1197
    αEEst[7] = 7267 // 94080
    αEEst2[1] = 2479 // 34992
    αEEst2[2] = 0
    αEEst2[3] = 123 // 416
    αEEst2[4] = 612941 // 3411720
    αEEst2[5] = 43 // 1440
    αEEst2[6] = 2272 // 6561
    αEEst2[7] = 79937 // 1113912
    αEEst2[8] = 3293 // 556956
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
            stability_size = 3.9879271987261333
        )
    )
end

"""
Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum, by P.W.Sharp and E.Smart,
 Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.

"""
function SharpSmart5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = zeros(T, 7, 7)
    c = zeros(T_time, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)
    αEEst2 = zeros(T, 7)
    c[2] = 16 // 105
    c[3] = 8 // 35
    c[4] = 9 // 20
    c[5] = 2 // 3
    c[6] = 7 // 9
    c[7] = 1
    A[2, 1] = 16 // 105
    A[3, 1] = 2 // 35
    A[3, 2] = 6 // 35
    A[4, 1] = 8793 // 40960
    A[4, 2] = -5103 // 8192
    A[4, 3] = 17577 // 20480
    A[5, 1] = 347 // 1458
    A[5, 2] = -7 // 20
    A[5, 3] = 3395 // 10044
    A[5, 4] = 49792 // 112995
    A[6, 1] = -1223224109959 // 9199771214400
    A[6, 2] = 1234787701 // 2523942720
    A[6, 3] = 568994101921 // 3168810084960
    A[6, 4] = -105209683888 // 891227836395
    A[6, 5] = 9 // 25
    A[7, 1] = 2462504862877 // 8306031988800
    A[7, 2] = -123991 // 287040
    A[7, 3] = 106522578491 // 408709510560
    A[7, 4] = 590616498832 // 804646848915
    A[7, 5] = -319138726 // 534081275
    A[7, 6] = 52758 // 71449
    α[1] = 1093 // 15120
    α[2] = 0
    α[3] = 60025 // 190992
    α[4] = 3200 // 20709
    α[5] = 1611 // 11960
    α[6] = 712233 // 2857960
    α[7] = 3 // 40
    αEEst[1] = 84018211 // 991368000
    αEEst[2] = 0
    αEEst[3] = 92098979 // 357791680
    αEEst[4] = 17606944 // 67891005
    αEEst[5] = 3142101 // 235253200
    αEEst[6] = 22004596809 // 70270091500
    αEEst[7] = 9 // 125
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
            stability_size = 3.9156746135081772
        )
    )
end

"""
BogakiShampine3()

Constructs the tableau object for the Bogakai-Shampine Order 2/3 method.
"""
function BogakiShampine3(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0
        1 // 2 0 0 0
        0 3 // 4 0 0
        2 // 9 1 // 3 4 // 9 0
    ]
    c = [0; 1 // 2; 3 // 4; 1]
    α = [2 // 9; 1 // 3; 4 // 9; 0]
    αEEst = [7 // 24; 1 // 4; 1 // 3; 1 // 8]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 3, αEEst = αEEst, adaptiveorder = 2,
            stability_size = 2.5127453266183286
        )
    )
end

"""
CashKarp()

Constructs the tableau object for the Cash-Karp Order 4/5 method.
"""
function CashKarp(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0 0 0
        1 // 5 0 0 0 0 0
        3 // 40 9 // 40 0 0 0 0
        3 // 10 -9 // 10 6 // 5 0 0 0
        -11 // 54 5 // 2 -70 // 27 35 // 27 0 0
        1631 // 55296 175 // 512 575 // 13824 44275 // 110592 253 // 4096 0
    ]
    c = [0; 1 // 5; 3 // 10; 3 // 5; 1; 7 // 8]
    α = [37 // 378; 0; 250 // 621; 125 // 594; 0; 512 // 1771]
    αEEst = [2825 // 27648; 0; 18575 // 48384; 13525 // 55296; 277 // 14336; 1 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
            stability_size = 3.7343596072347225
        )
    )
end

"""
Runge-Kutta-Fehberg Order 4/3
"""
function RKF4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    c = [0; 1 // 4; 4 // 9; 6 // 7; 1]
    A = [
        0 0 0 0 0
        1 // 4 0 0 0 0
        4 // 81 32 // 81 0 0 0
        57 // 98 -432 // 343 1053 // 686 0 0
        1 // 6 0 27 // 52 49 // 156 0
    ]
    α = [43 // 288; 0; 243 // 416; 343 // 1872; 1 // 12]
    αEEst = [1 // 6; 0; 27 // 52; 49 // 156; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 4, αEEst = αEEst, adaptiveorder = 3,
            stability_size = 4.109222736949077
        )
    )
end

"""
Butcher's Third Order 6

On Runge-Kutta Processes of High Order, by J. C. Butcher,
 Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194

"""
