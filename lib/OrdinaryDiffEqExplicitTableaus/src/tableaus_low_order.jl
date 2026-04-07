"""
Heun's Order 2 method.
"""
function Heun(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0
        1 0
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    αEEst = map(T, αEEst)
    return (
        DiffEqBase.ExplicitRKTableau(
            A, c, α, 2, αEEst = αEEst, adaptiveorder = 1,
            stability_size = 2.0
        )
    )
end

"""
Ralston's Order 2 method.
"""
function Ralston(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0
        2 // 3 0
    ]
    c = [0; 2 // 3]
    α = [1 // 4; 3 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 2, stability_size = 2.0))
end

"""
Euler's method.
"""
function Euler(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = Matrix{T}(undef, 1, 1)
    A[1] = 0
    c = [0]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 1, stability_size = 2.0))
end

"""
Kutta's Order 3 method.
"""
function Kutta3(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0
        1 // 2 0 0
        -1 2 0
    ]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 3, stability_size = 2.5127453266183286))
end

"""
Classic RK4 method.
"""
function RK4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0
        1 // 2 0 0 0
        0 1 // 2 0 0
        0 0 1 0
    ]
    c = [0; 1 // 2; 1 // 2; 1]
    α = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 4, stability_size = 2.785293563405282))
end

"""
Classic RK4 3/8's rule method.
"""
function RK438Rule(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0
        1 // 3 0 0 0
        -1 // 3 1 0 0
        1 -1 1 0
    ]
    c = [0; 1 // 3; 2 // 3; 1]
    α = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 4, stability_size = 2.785293563405282))
end

"""
Ralston's Order 4 method with minimum truncation error.
"""
function Ralston4(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    sqrt5 = sqrt(convert(T, 5))
    a21 = 4 // 10
    a31 = (-2889 + 1428 * sqrt5) / 1024
    a32 = (3785 - 1620 * sqrt5) / 1024
    a41 = (-3365 + 2094 * sqrt5) / 6040
    a42 = (-975 - 3046 * sqrt5) / 2552
    a43 = (467040 + 203968 * sqrt5) / 240845
    A = [
        0 0 0 0
        a21 0 0 0
        a31 a32 0 0
        a41 a42 a43 0
    ]
    b2 = 4 // 10
    b3 = (14 - 3 * sqrt5) / 16
    c = [0, b2, b3, 1]
    b1 = (263 + 24 * sqrt5) / 1812
    b2 = (125 - 1000 * sqrt5) / 3828
    b3 = 1024 * (3346 + 1623 * sqrt5) / 5924787
    b4 = (30 - 4 * sqrt5) / 123
    α = [b1, b2, b3, b4]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return DiffEqBase.ExplicitRKTableau(A, c, α, 4, stability_size = 2.7852935634052822)
end

"""
Explicit SSP method of order 2 using 2 stages.
"""
function SSPRK22(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0
        1 0
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 2, stability_size = -2.0))
end

"""
Explicit SSP method of order 3 using 3 stages.
"""
function SSPRK33(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0
        1 0 0
        1 // 4 1 // 4 0
    ]
    c = [0; 1; 1 // 2]
    α = [1 // 6; 1 // 6; 2 // 3]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 3, stability_size = 2.5127453266183286))
end

"""
Explicit SSP method of order 3 using 4 stages.
"""
function SSPRK43(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0
        1 // 2 0 0 0
        1 // 2 1 // 2 0 0
        1 // 6 1 // 6 1 // 6 0
    ]
    c = [0; 1 // 2; 1; 1 // 2]
    α = [1 // 6; 1 // 6; 1 // 6; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 3, stability_size = 5.149486147774043))
end

"""
Explicit SSP method of order 4 using 10 stages.
"""
function SSPRK104(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0 0 0 0 0 0 0 0 0
        1 // 6 0 0 0 0 0 0 0 0 0
        1 // 6 1 // 6 0 0 0 0 0 0 0 0
        1 // 6 1 // 6 1 // 6 0 0 0 0 0 0 0
        1 // 6 1 // 6 1 // 6 1 // 6 0 0 0 0 0 0
        1 // 15 1 // 15 1 // 15 1 // 15 1 // 15 0 0 0 0 0
        1 // 15 1 // 15 1 // 15 1 // 15 1 // 15 1 // 6 0 0 0 0
        1 // 15 1 // 15 1 // 15 1 // 15 1 // 15 1 // 6 1 // 6 0 0 0
        1 // 15 1 // 15 1 // 15 1 // 15 1 // 15 1 // 6 1 // 6 1 // 6 0 0
        1 // 15 1 // 15 1 // 15 1 // 15 1 // 15 1 // 6 1 // 6 1 // 6 1 // 6 0
    ]
    c = [0; 1 // 6; 1 // 3; 1 // 2; 2 // 3; 1 // 3; 1 // 2; 2 // 3; 5 // 6; 1]
    α = [
        1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10;
        1 // 10
    ]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 4, stability_size = 13.917047464637367))
end

"""
LobattoIIIC* Order 2 method

"""
function LobattoIIICStar2(::Type{T} = Float64, ::Type{T_time} = T) where {T, T_time}
    A = [
        0 0
        1 0
    ]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T_time, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 2, stability_size = -2.0))
end
