"""
    RoesslerRI

Holds the Butcher tableaus for a Roessler RI method. (high weak order)
"""
struct RoesslerRI{T, T2} <: Tableau
    c₀::Vector{T2}
    c₁::Vector{T2}
    c₂::Vector{T2}
    A₀::Matrix{T}
    A₁::Matrix{T}
    A₂::Matrix{T}
    B₀::Matrix{T}
    B₁::Matrix{T}
    B₂::Matrix{T}
    α::Vector{T}
    β₁::Vector{T}
    β₂::Vector{T}
    β₃::Vector{T}
    β₄::Vector{T}
    order::Rational{Int}
    quantile::T
end

function constructDRI1(T = Float64, T2 = Float64)
    c₀ = [0; 1 // 2; 1]
    c₁ = [0; 342 // 491; 342 // 491]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 // 2 0 0
        -1 2 0
    ]
    A₁ = [
        0 0 0
        342 // 491 0 0
        342 // 491 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        (6 - sqrt(6)) / 10 0 0
        (3 + 2 * sqrt(6)) / 5 0 0
    ]
    B₁ = [
        0 0 0
        3 * sqrt(38 // 491) 0 0
        -3 * sqrt(38 / 491) 0 0
    ]
    B₂ = [
        0 0 0
        -214 // 513 * sqrt(1105 // 991) -491 // 513 * sqrt(221 // 4955) -491 // 513 * sqrt(221 // 4955)
        214 // 513 * sqrt(1105 // 991) 491 // 513 * sqrt(221 // 4955) 491 // 513 * sqrt(221 // 4955)
    ]
    α = [1 // 6; 2 // 3; 1 // 6]

    β₁ = [193 // 684; 491 // 1368; 491 // 1368]
    β₂ = [0; 1 // 6 * sqrt(491 // 38); -1 // 6 * sqrt(491 // 38)]
    β₃ = [-4955 // 7072; 4955 // 14144; 4955 // 14144]
    β₄ = [0; -1 // 8 * sqrt(4955 // 221); 1 // 8 * sqrt(4955 // 221)]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRI1(T = Float64, T2 = Float64)
    c₀ = [0; 2 // 3; 2 // 3]
    c₁ = [0; 1; 1]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        2 // 3 0 0
        -1 // 3 1 0
    ]
    A₁ = [
        0 0 0
        1 0 0
        1 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        1 0 0
        0 0 0
    ]
    B₁ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    B₂ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    α = [1 // 4; 1 // 2; 1 // 4]

    β₁ = [1 // 2; 1 // 4; 1 // 4]
    β₂ = [0; 1 // 2; -1 // 2]
    β₃ = [-1 // 2; 1 // 4; 1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRI3(T = Float64, T2 = Float64)
    c₀ = [0; 1; 1 // 2]
    c₁ = [0; 1; 1]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 0 0
        1 // 4 1 // 4 0
    ]
    A₁ = [
        0 0 0
        1 0 0
        1 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        (3 - 2 * sqrt(6)) / 5 0 0
        (6 + sqrt(6)) / 10 0 0
    ]
    B₁ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    B₂ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    α = [1 // 6; 1 // 6; 2 // 3]

    β₁ = [1 // 2; 1 // 4; 1 // 4]
    β₂ = [0; 1 // 2; -1 // 2]
    β₃ = [-1 // 2; 1 // 4; 1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRI5(T = Float64, T2 = Float64)
    c₀ = [0; 1; 5 // 12]
    c₁ = [0; 1 // 4; 1 // 4]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 0 0
        25 // 144 35 // 144 0
    ]
    A₁ = [
        0 0 0
        1 // 4 0 0
        1 // 4 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        1 // 3 0 0
        -5 // 6 0 0
    ]
    B₁ = [
        0 0 0
        1 // 2 0 0
        -1 // 2 0 0
    ]
    B₂ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    α = [1 // 10; 3 // 14; 24 // 35]

    β₁ = [1; -1; -1]
    β₂ = [0; 1; -1]
    β₃ = [1 // 2; -1 // 4; -1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRI6(T = Float64, T2 = Float64)
    c₀ = [0; 1; 0]
    c₁ = [0; 1; 1]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 0 0
        0 0 0
    ]
    A₁ = [
        0 0 0
        1 0 0
        1 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        1 0 0
        0 0 0
    ]
    B₁ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    B₂ = [
        0 0 0
        1 0 0
        -1 0 0
    ]
    α = [1 // 2; 1 // 2; 0]

    β₁ = [1 // 2; 1 // 4; 1 // 4]
    β₂ = [0; 1 // 2; -1 // 2]
    β₃ = [-1 // 2; 1 // 4; 1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRDI1WM(T = Float64, T2 = Float64)
    c₀ = [0; 2 // 3]
    c₁ = [0; 0]
    c₂ = [0; 0]
    A₀ = [
        0 0
        2 // 3 0
    ]
    A₁ = [
        0 0
        0 0
    ]
    A₂ = [
        0 0
        0 0
    ]
    B₀ = [
        0 0
        2 // 3 0
    ]
    B₁ = [
        0 0
        0 0
    ]
    B₂ = [
        0 0
        0 0
    ]
    α = [1 // 4; 3 // 4]

    β₁ = [1; 0]
    β₂ = [0; 0]
    β₃ = [0; 0]
    β₄ = [0; 0]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRDI2WM(T = Float64, T2 = Float64)
    c₀ = [0; 1; 0]
    c₁ = [0; 2 // 3; 2 // 3]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 0 0
        0 0 0
    ]
    A₁ = [
        0 0 0
        2 // 3 0 0
        2 // 3 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        1 0 0
        0 0 0
    ]
    B₁ = [
        0 0 0
        sqrt(2 // 3) 0 0
        -sqrt(2 // 3) 0 0
    ]
    B₂ = [
        0 0 0
        sqrt(2) 0 0
        -sqrt(2) 0 0
    ]
    α = [1 // 2; 1 // 2; 0]

    β₁ = [1 // 4; 3 // 8; 3 // 8]
    β₂ = [0; sqrt(6) / 4; -sqrt(6) / 4]
    β₃ = [-1 // 4; 1 // 8; 1 // 8]
    β₄ = [0; sqrt(2) / 4; -sqrt(2) / 4]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRDI3WM(T = Float64, T2 = Float64)
    c₀ = [0; 1 // 2; 3 // 4]
    c₁ = [0; 2 // 3; 2 // 3]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 // 2 0 0
        0 3 // 4 0
    ]
    A₁ = [
        0 0 0
        2 // 3 0 0
        2 // 3 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        (9 - 2 * sqrt(15)) / 14 0 0
        (18 + 3 * sqrt(15)) / 28 0 0
    ]
    B₁ = [
        0 0 0
        sqrt(2 // 3) 0 0
        -sqrt(2 // 3) 0 0
    ]
    B₂ = [
        0 0 0
        sqrt(2) 0 0
        -sqrt(2) 0 0
    ]
    α = [2 // 9; 1 // 3; 4 // 9]

    β₁ = [1 // 4; 3 // 8; 3 // 8]
    β₂ = [0; sqrt(6) / 4; -sqrt(6) / 4]
    β₃ = [-1 // 4; 1 // 8; 1 // 8]
    β₄ = [0; sqrt(2) / 4; -sqrt(2) / 4]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function constructRDI4WM(T = Float64, T2 = Float64)
    c₀ = [0; 1 // 2; 1]
    c₁ = [0; 2 // 3; 2 // 3]
    c₂ = [0; 0; 0]
    A₀ = [
        0 0 0
        1 // 2 0 0
        -1 2 0
    ]
    A₁ = [
        0 0 0
        2 // 3 0 0
        2 // 3 0 0
    ]
    A₂ = [
        0 0 0
        0 0 0
        0 0 0
    ]
    B₀ = [
        0 0 0
        (6 - sqrt(6)) / 10 0 0
        (3 + 2 * sqrt(6)) / 5 0 0
    ]
    B₁ = [
        0 0 0
        sqrt(2 // 3) 0 0
        -sqrt(2 // 3) 0 0
    ]
    B₂ = [
        0 0 0
        sqrt(2) 0 0
        -sqrt(2) 0 0
    ]
    α = [1 // 6; 2 // 3; 1 // 6]

    β₁ = [1 // 4; 3 // 8; 3 // 8]
    β₂ = [0; sqrt(6) / 4; -sqrt(6) / 4]
    β₃ = [-1 // 4; 1 // 8; 1 // 8]
    β₄ = [0; sqrt(2) / 4; -sqrt(2) / 4]
    return RoesslerRI(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014
    )
end

function checkRIOrder(RI; tol = 1.0e-6, ps = 2)
    (; c₀, c₁, c₂, A₀, A₁, A₂, B₀, B₁, B₂, α, β₁, β₂, β₃, β₄) = RI
    e = ones(size(α))
    if ps == 2
        conditions = Vector{Bool}(undef, 59) # 9 conditions for first order, 59 in total
    elseif ps == 1
        conditions = Vector{Bool}(undef, 9)
    else
        error("Only the conditions for first and second order weak convergence are implemented. Choose ps=1 or ps=2.")
    end
    #first order weak sense
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₄, e) - 0) < tol
    conditions[3] = abs(dot(β₃, e) - 0) < tol
    conditions[4] = abs(dot(β₁, e)^2 - 1) < tol
    conditions[5] = abs(dot(β₂, e) - 0) < tol
    conditions[6] = abs(dot(β₁, (B₁ * e)) - 0) < tol
    conditions[7] = abs(dot(β₄, (A₂ * e)) - 0) < tol
    conditions[8] = abs(dot(β₃, (B₂ * e)) - 0) < tol
    conditions[9] = abs(dot(β₄, (B₂ * e) .^ 2) - 0) < tol

    #second order weak sense
    if ps == 2
        conditions[10] = abs(dot(α, (A₀ * e)) - 0.5) < tol
        conditions[11] = abs(dot(α, (B₀ * e) .^ 2) - 0.5) < tol
        conditions[12] = abs(dot(β₁, e) * dot(α, (B₀ * e)) - 0.5) < tol
        conditions[13] = abs(dot(β₁, e) * dot(β₁, (A₁ * e)) - 0.5) < tol
        conditions[14] = abs(dot(β₃, (A₂ * e)) - 0) < tol
        conditions[15] = abs(dot(β₂, (B₁ * e)) - 1) < tol
        conditions[16] = abs(dot(β₄, (B₂ * e)) - 1) < tol
        conditions[17] = abs(dot(β₁, e) * dot(β₁, (B₁ * e) .^ 2) - 0.5) < tol
        conditions[18] = abs(dot(β₁, e) * dot(β₃, (B₂ * e) .^ 2) - 0.5) < tol
        conditions[19] = abs(dot(β₁, B₁ * (B₁ * e)) - 0) < tol
        conditions[20] = abs(dot(β₃, B₂ * (B₁ * e)) - 0) < tol
        conditions[21] = abs(dot(β₃, B₂ * (B₁ * (B₁ * e))) - 0) < tol
        conditions[22] = abs(dot(β₁, A₁ * (B₀ * e)) - 0) < tol
        conditions[23] = abs(dot(β₃, A₂ * (B₀ * e)) - 0) < tol
        conditions[24] = abs(dot(β₄, (A₂ * e) .^ 2) - 0) < tol
        conditions[25] = abs(dot(β₄, A₂ * (A₀ * e)) - 0) < tol
        conditions[26] = abs(dot(α, B₀ * (B₁ * e)) - 0) < tol
        conditions[27] = abs(dot(β₂, A₁ * e) - 0) < tol
        conditions[28] = abs(dot(β₁, (A₁ * e) .* (B₁ * e)) - 0) < tol
        conditions[29] = abs(dot(β₃, (A₂ * e) .* (B₂ * e)) - 0) < tol
        conditions[30] = abs(dot(β₄, A₂ * (B₀ * e)) - 0) < tol
        conditions[31] = abs(dot(β₂, A₁ * (B₀ * e)) - 0) < tol
        conditions[32] = abs(dot(β₄, ((B₂ * e) .^ 2) .* (A₂ * e)) - 0) < tol
        conditions[33] = abs(dot(β₄, A₂ * (B₀ * e) .^ 2) - 0) < tol
        conditions[34] = abs(dot(β₂, A₁ * (B₀ * e) .^ 2) - 0) < tol
        conditions[35] = abs(dot(β₁, B₁ * (A₁ * e)) - 0) < tol
        conditions[36] = abs(dot(β₃, B₂ * (A₁ * e)) - 0) < tol
        conditions[37] = abs(dot(β₂, (B₁ * e) .^ 2) - 0) < tol
        conditions[38] = abs(dot(β₄, B₂ * (B₁ * e)) - 0) < tol
        conditions[39] = abs(dot(β₂, B₁ * (B₁ * e)) - 0) < tol
        conditions[40] = abs(dot(β₁, (B₁ * e) .^ 3) - 0) < tol
        conditions[41] = abs(dot(β₃, (B₂ * e) .^ 3) - 0) < tol
        conditions[42] = abs(dot(β₁, B₁ * ((B₁ * e) .^ 2)) - 0) < tol
        conditions[43] = abs(dot(β₃, B₂ * ((B₁ * e) .^ 2)) - 0) < tol
        conditions[44] = abs(dot(β₄, (B₂ * e) .^ 4) - 0) < tol
        conditions[45] = abs(dot(β₄, (B₂ * (B₁ * e)) .^ 2) - 0) < tol
        conditions[46] = abs(dot(β₄, (B₂ * e) .* (B₂ * (B₁ * e))) - 0) < tol
        conditions[47] = abs(dot(α, (B₀ * e) .* (B₀ * (B₁ * e))) - 0) < tol
        conditions[48] = abs(dot(β₁, (A₁ * (B₀ * e)) .* (B₁ * e)) - 0) < tol
        conditions[49] = abs(dot(β₃, (A₂ * (B₀ * e)) .* (B₂ * e)) - 0) < tol
        conditions[50] = abs(dot(β₁, A₁ * (B₀ * (B₁ * e))) - 0) < tol
        conditions[51] = abs(dot(β₃, A₂ * (B₀ * (B₁ * e))) - 0) < tol
        conditions[52] = abs(dot(β₄, (B₂ * (A₁ * e)) .* (B₂ * e)) - 0) < tol
        conditions[53] = abs(dot(β₁, B₁ * (A₁ * (B₀ * e))) - 0) < tol
        conditions[54] = abs(dot(β₃, B₂ * (A₁ * (B₀ * e))) - 0) < tol
        conditions[55] = abs(dot(β₁, (B₁ * e) .* (B₁ * (B₁ * e))) - 0) < tol
        conditions[56] = abs(dot(β₃, (B₂ * e) .* (B₂ * (B₁ * e))) - 0) < tol
        conditions[57] = abs(dot(β₁, B₁ * (B₁ * (B₁ * e))) - 0) < tol
        conditions[58] = abs(dot(β₄, (B₂ * e) .* (B₂ * ((B₁ * e) .^ 2))) - 0) < tol
        conditions[59] = abs(dot(β₄, (B₂ * e) .* (B₂ * (B₁ * (B₁ * e)))) - 0) < tol
    end
    return (conditions)
end

"""
    RoesslerRS

Holds the Butcher tableaus for a Roessler RS method. (high weak order in Stratonovich sense)
"""
struct RoesslerRS{T, T2} <: Tableau
    c₀::Vector{T2}
    c₁::Vector{T2}
    c₂::Vector{T2}
    A₀::Matrix{T}
    A₁::Matrix{T}
    A₂::Matrix{T}
    B₀::Matrix{T}
    B₁::Matrix{T}
    B₂::Matrix{T}
    B₃::Matrix{T}
    α::Vector{T}
    β₁::Vector{T}
    β₂::Vector{T}
    order::Rational{Int}
    quantile::T
end

function constructRS1(T = Float64, T2 = Float64)
    c₀ = [0; 0; 1; 0]
    c₁ = [0; 0; 1; 1]
    c₂ = [0; 0; 0; 0]
    A₀ = [
        0 0 0 0
        0 0 0 0
        1 0 0 0
        0 0 0 0
    ]
    A₁ = [
        0 0 0 0
        0 0 0 0
        1 0 0 0
        1 0 0 0
    ]
    A₂ = [
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    B₀ = [
        0 0 0 0
        0 0 0 0
        1 // 4 3 // 4 0 0
        0 0 0 0
    ]
    B₁ = [
        0 0 0 0
        2 // 3 0 0 0
        1 // 12 1 // 4 0 0
        -5 // 4 1 // 4 2 0
    ]
    B₂ = [
        0 0 0 0
        1 0 0 0
        -1 0 0 0
        0 0 0 0
    ]
    B₃ = [
        0 0 0 0
        0 0 0 0
        1 // 4 3 // 4 0 0
        1 // 4 3 // 4 0 0
    ]
    α = [0; 0; 1 // 2; 1 // 2]

    β₁ = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    β₂ = [0; -1 // 4; 1 // 4; 0]

    return RoesslerRS(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂), map(T, B₃),
        map(T, α), map(T, β₁), map(T, β₂),
        1 // 1, -0.9674215661017014
    )
end

function constructRS2(T = Float64, T2 = Float64)
    c₀ = [0; 2 // 3; 2 // 3; 0]
    c₁ = [0; 0; 1; 1]
    c₂ = [0; 0; 0; 0]
    A₀ = [
        0 0 0 0
        2 // 3 0 0 0
        1 // 6 1 // 2 0 0
        0 0 0 0
    ]
    A₁ = [
        0 0 0 0
        0 0 0 0
        1 0 0 0
        1 0 0 0
    ]
    A₂ = [
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    B₀ = [
        0 0 0 0
        0 0 0 0
        1 // 4 3 // 4 0 0
        0 0 0 0
    ]
    B₁ = [
        0 0 0 0
        2 // 3 0 0 0
        1 // 12 1 // 4 0 0
        -5 // 4 1 // 4 2 0
    ]
    B₂ = [
        0 0 0 0
        1 0 0 0
        -1 0 0 0
        0 0 0 0
    ]
    B₃ = [
        0 0 0 0
        0 0 0 0
        1 // 4 3 // 4 0 0
        1 // 4 3 // 4 0 0
    ]
    α = [1 // 4; 1 // 4; 1 // 2; 0]

    β₁ = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    β₂ = [0; -1 // 4; 1 // 4; 0]

    return RoesslerRS(
        map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂), map(T, B₃),
        map(T, α), map(T, β₁), map(T, β₂),
        1 // 1, -0.9674215661017014
    )
end

function checkRSOrder(RS; tol = 1.0e-6, ps = 2)
    (; c₀, c₁, c₂, A₀, A₁, A₂, B₀, B₁, B₂, B₂, B₃, α, β₁, β₂) = RS
    e = ones(size(α))
    if ps == 2
        conditions = Vector{Bool}(undef, 55) # 6 conditions for first order, 55 in total
    elseif ps == 1
        conditions = Vector{Bool}(undef, 6)
    else
        error("Only the conditions for first and second order weak convergence are implemented. Choose ps=1 or ps=2.")
    end
    #first order weak sense
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₁, e)^2 - 1) < tol
    conditions[3] = abs(dot(β₂, e) - 0) < tol
    conditions[4] = abs(dot(β₁, (B₁ * e)) - 0.5) < tol
    conditions[5] = abs(dot(β₂, (A₂ * e)) - 0) < tol
    conditions[6] = abs(dot(β₂, (B₂ * e) .^ 2) - 0) < tol

    #second order weak sense
    if ps == 2
        conditions[7] = abs(dot(α, (A₀ * e)) - 0.5) < tol
        conditions[8] = abs(dot(α, B₀ * (B₁ * e)) - 0.25) < tol
        conditions[9] = abs(dot(α, (B₀ * e) .^ 2) - 0.5) < tol
        conditions[10] = abs(dot(β₁, e) * dot(α, (B₀ * e)) - 0.5) < tol
        conditions[11] = abs(dot(β₁, e) * dot(β₁, (A₁ * e)) - 0.5) < tol
        conditions[12] = abs(dot(β₁, B₁ * (A₁ * e)) - 0.25) < tol
        conditions[13] = abs(dot(β₁, (B₁ * e) .* (A₁ * e)) - 0.25) < tol
        conditions[14] = abs(dot(β₁, (B₃ * e)) - 0.5) < tol
        conditions[15] = abs(dot(β₁, e) * dot(β₁, (B₁ * e) .^ 2) - 1 / 3) < tol
        conditions[16] = abs(dot(β₁, e) * dot(β₁, (B₃ * e) .^ 2) - 0.5) < tol
        conditions[17] = abs(dot(β₁, B₃ * (B₃ * e)) - 0) < tol
        conditions[18] = abs(dot(β₂, (B₂ * e))^2 - 0.25) < tol
        conditions[19] = abs(dot(β₁, (B₁ * e) .^ 3) - 0.25) < tol
        conditions[20] = abs(dot(β₁, B₁ * (B₁ * e) .^ 2) - 1 / 12) < tol
        conditions[21] = abs(dot(β₁, B₁ * (B₃ * e) .^ 2) - 0.25) < tol
        conditions[22] = abs(dot(β₁, A₁ * (B₀ * e)) - 0) < tol
        conditions[23] = abs(dot(β₂, (A₂ * e) .^ 2) - 0) < tol
        conditions[24] = abs(dot(β₂, A₂ * (A₀ * e)) - 0) < tol
        conditions[25] = abs(dot(β₁, B₁ * B₁ * (B₁ * e)) - 1 / 24) < tol
        conditions[26] = abs(dot(β₂, A₂ * (B₀ * e)) - 0) < tol
        conditions[27] = abs(dot(β₂, A₂ * (B₀ * e) .^ 2) - 0) < tol
        conditions[28] = abs(dot(β₂, (B₂ * e) .^ 4) - 0) < tol
        conditions[29] = abs(dot(β₂, (B₂ * (B₁ * e)) .^ 2) - 0) < tol
        conditions[30] = abs(dot(β₂, (B₂ * (B₃ * e)) .^ 2) - 0) < tol
        conditions[31] = abs(dot(β₁, (B₁ * e) .* (B₃ * e) .^ 2) - 0.25) < tol
        conditions[32] = abs(dot(β₂, (A₂ * e) .* (B₂ * e) .^ 2) - 0) < tol
        conditions[33] = abs(dot(β₁, B₁ * B₃ * (B₁ * e)) - 1 // 8) < tol
        conditions[34] = abs(dot(β₁, B₃ * B₃ * (B₃ * e)) - 0) < tol
        conditions[35] = abs(dot(β₁, B₃ * B₁ * (B₃ * e)) - 0) < tol
        conditions[36] = abs(dot(β₂, A₂ * B₀ * (B₁ * e)) - 0) < tol
        conditions[37] = abs(dot(β₁, e) * dot(β₁, (B₃ * e) .* (B₁ * e)) - 0.25) < tol
        conditions[38] = abs(dot(β₁, e) * dot(β₁, B₁ * (B₁ * e)) - 1 / 6) < tol
        conditions[39] = abs(dot(β₁, e) * dot(β₁, B₃ * (B₁ * e)) - 0.25) < tol
        conditions[40] = abs(dot(β₁, e) * dot(β₁, B₁ * (B₃ * e)) - 0.25) < tol
        conditions[41] = abs(dot(β₁, (B₁ * e) .* (B₁ * (B₁ * e))) - 1 / 8) < tol
        conditions[42] = abs(dot(β₁, (B₁ * e) .* (B₃ * (B₁ * e))) - 1 / 8) < tol
        conditions[43] = abs(dot(β₁, (B₃ * e) .* (B₁ * (B₃ * e))) - 1 / 4) < tol
        conditions[44] = abs(dot(β₁, (B₃ * e) .* (B₃ * (B₃ * e))) - 0) < tol
        conditions[45] = abs(dot(β₁, B₃ * ((B₃ * e) .* (B₁ * e))) - 0) < tol
        conditions[46] = abs(dot(β₂, (B₂ * (A₁ * e)) .* (B₂ * e)) - 0) < tol
        conditions[47] = abs(dot(β₂, (B₂ * e) .* (B₂ * (B₁ * e))) - 0) < tol
        conditions[48] = abs(dot(β₂, (B₂ * e) .* (B₂ * (B₃ * e))) - 0) < tol
        conditions[49] = abs(dot(β₂, (B₂ * e) .* (B₂ * ((B₁ * e) .^ 2))) - 0) < tol
        conditions[50] = abs(dot(β₂, (B₂ * e) .* (B₂ * ((B₃ * e) .^ 2))) - 0) < tol
        conditions[51] = abs(dot(β₂, (B₂ * e) .* (B₂ * ((B₁ * e) .* (B₃ * e)))) - 0) < tol
        conditions[52] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₁ * (B₁ * e))) - 0) < tol
        conditions[53] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₃ * (B₁ * e))) - 0) < tol
        conditions[54] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₃ * (B₃ * e))) - 0) < tol
        conditions[55] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₁ * (B₃ * e))) - 0) < tol
    end
    return (conditions)
end

"""
    KomoriNON

Holds the Butcher tableau for Komori's NON method. (second weak order in Stratonovich sense)
"""
struct KomoriNON{T} <: Tableau
    c0::Vector{T}
    cj::Vector{T}
    cjl::Vector{T}
    clj::Vector{T}

    α00::Matrix{T}
    α0j::Matrix{T}
    αj0::Matrix{T}
    αjj::Matrix{T}
    αjl::Matrix{T}

    αjljj::Matrix{T}
    αljjl::Matrix{T}

    order::Rational{Int}
    quantile::T
end

function constructNON(T = Float64)
    c0 = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    cj = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    cjl = [0; 1 // 2; -1 // 2; 0]
    clj = [0; -1 // 2; 1 // 2; 0]

    α00 = [
        0 0 0 0
        1 // 2 0 0 0
        0 1 // 2 0 0
        0 0 1 0
    ]

    α0j = [
        0 0 0 0
        1 0 0 0
        -9 // 8 9 // 8 0 0
        1 0 0 0
    ]

    αj0 = [
        0 0 0 0
        2 0 0 0
        0 0 0 0
        -2 0 0 0
    ]

    αjj = [
        0 0 0 0
        2 // 3 0 0 0
        1 // 12 1 // 4 0 0
        -5 // 4 1 // 4 2 0
    ]

    αjl = [
        0 0 0 0
        0 0 0 0
        1 // 4 3 // 4 0 0
        1 // 4 3 // 4 0 0
    ]

    αjljj = [
        0 0 0 0
        0 0 0 0
        1 0 0 0
        0 0 0 0
    ]

    αljjl = [
        0 0 0 0
        -1 // 2 0 0 0
        1 // 2 0 0 0
        0 0 0 0
    ]

    return KomoriNON(
        map(T, c0), map(T, cj), map(T, cjl),
        map(T, clj), map(T, α00), map(T, α0j),
        map(T, αj0), map(T, αjj), map(T, αjl), map(T, αjljj), map(T, αljjl),
        1 // 1, -0.9674215661017014
    )
end

function checkNONOrder(NON; tol = 1.0e-6)
    if NON isa KomoriNON
        (; c0, cj, cjl, clj, α00, α0j, αj0, αjj, αjl, αjljj, αljjl) = NON
    elseif NON isa KomoriNON2
        (; c0, cj, ckj, α00, α0j, αj0, αjj, αjl, αkjjl) = NON
    else
        (; c0, cj, α00, α0j, αj0, αjj, αjl) = NON
    end
    e = ones(size(c0))

    conditions = Vector{Bool}(undef, 44) # 38 conditions for second order, 6 extra conditions for fourth deterministic order, 44 in total
    conditions[1] = abs(dot(c0, e) - 1) < tol
    conditions[2] = abs(dot(cj, e) - 1) < tol
    conditions[3] = abs(dot(cj, αjj * e) - 1 / 2) < tol
    conditions[4] = abs(dot(cj, α0j * e) - 1 / 2) < tol
    conditions[5] = abs(dot(cj, αj0 * e) - 1 / 2) < tol
    conditions[6] = abs(dot(c0, α00 * e) - 1 / 2) < tol
    conditions[7] = abs(dot(cj, αjj * (αj0 * e)) - 1 / 4) < tol

    conditions[8] = abs(dot(c0, α0j * (αjj * e)) - 1 / 4) < tol
    conditions[9] = abs(dot(cj, αj0 * (α0j * e)) - 0) < tol
    conditions[10] = abs(dot(c0, (α0j * e) .^ 2) - 1 / 2) < tol
    conditions[11] = abs(dot(cj, (αj0 * e) .* (αjj * e)) - 1 / 4) < tol
    conditions[12] = abs(dot(cj, αjj * αjj * (αjj * e)) - 1 / 24) < tol

    conditions[13] = abs(dot(cj, αjj * (αjj * e) .^ 2) - 1 / 12) < tol
    conditions[14] = abs(dot(cj .* (αjj * e), αjj * (αjj * e)) - 1 / 8) < tol
    conditions[15] = abs(dot(cj, (αjj * e) .^ 3) - 1 / 4) < tol
    conditions[16] = abs(dot(cj, αjj * (αjj * e)) - 1 / 6) < tol
    conditions[17] = abs(dot(cj, (αjj * e) .^ 2) - 1 / 3) < tol
    conditions[18] = abs(dot(cj, (αjl * e)) - 1 / 2) < tol
    conditions[19] = abs(dot(cj .* (αjl * e), αjl * (αjl * e)) - 0) < tol
    conditions[20] = abs(dot(cj, (αjl * e) .^ 2) - 1 / 2) < tol
    conditions[21] = abs(dot(cj, αjj * αjl * (αjj * e)) - 1 / 8) < tol
    conditions[22] = abs(dot(cj, αjl * αjl * (αjl * e)) - 0) < tol
    conditions[23] = abs(dot(cj, αjl * αjj * (αjl * e)) - 0) < tol
    conditions[24] = abs(dot(cj, αjj * (αjl * e) .^ 2) - 1 / 4) < tol
    conditions[25] = abs(dot(cj, αjl * ((αjj * e) .* (αjl * e))) - 0) < tol
    conditions[26] = abs(dot(cj .* (αjj * e), αjl * (αjj * e)) - 1 / 8) < tol
    conditions[27] = abs(dot(cj .* (αjl * e), αjj * (αjl * e)) - 1 / 4) < tol
    conditions[28] = abs(dot(cj .* (αjj * e), (αjl * e) .^ 2) - 1 / 4) < tol
    conditions[29] = abs(dot(cj, αjj * (αjl * e)) - 1 / 4) < tol
    conditions[30] = abs(dot(cj, αjl * (αjj * e)) - 1 / 4) < tol
    conditions[31] = abs(dot(cj, αjl * (αjl * e)) - 0) < tol
    conditions[32] = abs(dot(cj, (αjj * e) .* (αjl * e)) - 1 / 4) < tol

    if NON isa KomoriNON
        conditions[33] = abs(dot(clj, e) - 0) < tol
        conditions[34] = abs(dot(cjl, e) - 0) < tol
        conditions[35] = abs(dot(clj, αljjl * e) - 1 / 2) < tol
        conditions[36] = abs(dot(cjl, αjljj * e) + 1 / 2) < tol
        conditions[37] = abs(dot(clj .* (αljjl * e), αljjl * (αjljj * e)) - 0) < tol
        conditions[38] = abs(dot(clj, (αljjl * e) .^ 2) - 0) < tol
    elseif NON isa KomoriNON2
        conditions[33] = abs(ckj[4] + ckj[3] - 0) < tol
        conditions[34] = abs(ckj[2] - 0) < tol
        conditions[35] = abs(αkjjl[4, 3] - 0) < tol
        conditions[36] = true # we set alpha^(k(j),j,0,0) = 0
        conditions[37] = abs(ckj[3] * αkjjl[3, 2]^2 + ckj[4] * αkjjl[4, 2]^2 - 0) < tol
        conditions[38] = abs(ckj[3] * αkjjl[3, 2] + ckj[4] * αkjjl[4, 2] - 1 / 2) < tol
    end
    # deterministic fourth order
    conditions[39] = abs(dot(c0, α00 * α00 * (α00 * e)) - 1 / 24) < tol
    conditions[40] = abs(dot(c0, α00 * (α00 * e) .^ 2) - 1 / 12) < tol
    conditions[41] = abs(dot(c0 .* (α00 * e), α00 * (α00 * e)) - 1 / 8) < tol
    conditions[42] = abs(dot(c0, (α00 * e) .^ 3) - 1 / 4) < tol
    conditions[43] = abs(dot(c0, α00 * (α00 * e)) - 1 / 6) < tol
    conditions[44] = abs(dot(c0, (α00 * e) .^ 2) - 1 / 3) < tol

    return (conditions)
end

struct KomoriNON2{T} <: Tableau
    c0::Vector{T}
    cj::Vector{T}
    ckj::Vector{T}

    α00::Matrix{T}
    α0j::Matrix{T}
    αj0::Matrix{T}
    αjj::Matrix{T}
    αjl::Matrix{T}
    αkjjl::Matrix{T}

    order::Rational{Int}
    quantile::T
end

function constructNON2(T = Float64)
    # gamma is a free parameter
    γ = 1
    c0 = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    cj = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    ckj = [0; 0; γ; -γ]

    α00 = [
        0 0 0 0
        1 // 2 0 0 0
        0 1 // 2 0 0
        0 0 1 0
    ]

    α0j = [
        0 0 0 0
        1 0 0 0
        -9 // 8 9 // 8 0 0
        1 0 0 0
    ]

    αj0 = [
        0 0 0 0
        2 0 0 0
        0 0 0 0
        -2 0 0 0
    ]

    αjj = [
        0 0 0 0
        2 // 3 0 0 0
        1 // 12 1 // 4 0 0
        -5 // 4 1 // 4 2 0
    ]

    αjl = [
        0 0 0 0
        0 0 0 0
        1 // 4 3 // 4 0 0
        1 // 4 3 // 4 0 0
    ]

    αkjjl = [
        0 0 0 0
        0 0 0 0
        0 1 / (4 * γ) 0 0
        0 -1 / (4 * γ) 0 0
    ]

    return KomoriNON2(
        map(T, c0), map(T, cj), map(T, ckj),
        map(T, α00), map(T, α0j),
        map(T, αj0), map(T, αjj), map(T, αjl), map(T, αkjjl),
        1 // 1, -0.9674215661017014
    )
end
