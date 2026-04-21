"""
Gauss-Legendre Order 2.
"""
function GL2(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    c = [1 // 2]
    A = fill(1 // 2, 1, 1)
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 2))
end

"""
Gauss-Legendre Order 4.
"""
function GL4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    c = [(3 - sqrt(3)) / 6; (3 + sqrt(3)) / 6]
    A = [
        1 / 4 (3 - 2 * sqrt(3)) / 12;
        (3 + 2 * sqrt(3)) / 12 1 / 4
    ]
    α = [1 / 2; 1 / 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 4))
end

"""
Gauss-Legendre Order 6.
"""
function GL6(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    c = [(5 - sqrt(15)) / 10; 1 / 2; (5 + sqrt(15)) / 10]
    A = [
        5 / 36 (10 - 3 * sqrt(15)) / 45 (25 - 6 * sqrt(15)) / 180;
        (10 + 3 * sqrt(15)) / 72 2 / 9 (10 - 3 * sqrt(15)) / 72;
        (25 + 6 * sqrt(15)) / 180 (10 + 3 * sqrt(15)) / 45 5 / 36
    ]
    α = [5 / 18; 4 / 9; 5 / 18]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 6))
end

"""
Implicit Euler Method
"""
function ImplicitEuler(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = Matrix{T}(undef, 1, 1)
    A[1] = 1
    c = [1]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 1))
end

"""
Order 2 Midpoint Method
"""
function MidpointRule(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = Matrix{T}(undef, 1, 1)
    A[1] = 1 // 2
    c = [1 // 2]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 2))
end

"""
Order 2 Trapezoidal Rule (LobattoIIIA2)
"""
function TrapezoidalRule(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0
        1 // 2 1 // 2
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1))
end

"""
LobattoIIIA Order 4 method
"""
function LobattoIIIA4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0
        5 // 24 1 // 3 -1 // 24
        1 // 6 2 // 3 1 // 6
    ]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIIB Order 2 method
"""
function LobattoIIIB2(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 2 0
        1 // 2 0
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1))
end

"""
LobattoIIIB Order 4 method

"""
function LobattoIIIB4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 6 -1 // 6 0
        1 // 6 1 // 3 0
        1 // 6 5 // 6 0
    ]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIIC Order 2 method

"""
function LobattoIIIC2(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 2 -1 // 2
        1 // 2 1 // 2
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1))
end

"""
LobattoIIIC Order 4 method

"""
function LobattoIIIC4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 6 -1 // 3 1 // 6
        1 // 6 5 // 12 -1 // 12
        1 // 6 2 // 3 1 // 6
    ]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIIC* Order 4 method

"""
function LobattoIIICStar4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0
        1 // 4 1 // 4 0
        0 1 0
    ]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIID Order 2 method

"""
function LobattoIIID2(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 2 1 // 2
        -1 // 2 1 // 2
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 2))
end

"""
LobattoIIID Order 4 method

"""
function LobattoIIID4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 6 0 -1 // 6
        1 // 12 5 // 12 0
        1 // 2 1 // 3 1 // 6
    ]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
RadauIA Order 3 method

"""
function RadauIA3(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 4 -1 // 4
        1 // 4 5 // 12
    ]
    c = [0; 2 // 3]
    α = [1 // 4; 3 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 3))
end

"""
RadauIA Order 5 method

"""
function RadauIA5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        1 // 9 (-1 - sqrt(6)) / 18 (-1 + sqrt(6)) / 18
        1 // 9 11 // 45 + 7 * sqrt(6) / 360 11 // 45 - 43 * sqrt(6) / 360
        1 // 9 11 // 45 + 43 * sqrt(6) / 360 11 // 45 - 7 * sqrt(6) / 360
    ]
    c = [0; 3 // 5 - sqrt(6) / 10; 3 // 5 + sqrt(6) / 10]
    α = [1 // 9; 4 // 9 + sqrt(6) / 36; 4 // 9 - sqrt(6) / 36]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 5))
end

"""
RadauIIA Order 3 method

"""
function RadauIIA3(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        5 // 12 -1 // 12
        3 // 4 1 // 4
    ]
    c = [1 // 3; 1]
    α = [3 // 4; 1 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 3))
end

"""
RadauIIA Order 5 method

"""
function RadauIIA5(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    sq6 = sqrt(6)
    A = [
        11 // 45 - 7sq6 / 360 37 // 225 - 169sq6 / 1800 -2 // 225 + sq6 / 75
        37 // 225 + 169sq6 / 1800 11 // 45 + 7sq6 / 360 -2 // 225 - sq6 / 75
        4 // 9 - sq6 / 36 4 // 9 + sq6 / 36 1 // 9
    ]
    c = [2 // 5 - sq6 / 10; 2 / 5 + sq6 / 10; 1]
    α = [4 // 9 - sq6 / 36; 4 // 9 + sq6 / 36; 1 // 9]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ImplicitRKTableau(A, c, α, 5))
end
