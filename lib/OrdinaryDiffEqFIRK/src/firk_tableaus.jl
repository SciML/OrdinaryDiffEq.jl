struct RadauIIA3Tableau{T, T2}
    T11::T
    T12::T
    T21::T
    T22::T
    TI11::T
    TI12::T
    TI21::T
    TI22::T
    c1::T2
    c2::T2
    α::T
    β::T
    e1::T
    e2::T
end

function RadauIIA3Tableau(T, T2)
    T11 = T(0.10540925533894596)
    T12 = T(-0.29814239699997197)
    T21 = T(0.9486832980505138)
    T22 = T(0.0)
    TI11 = T(0.0)
    TI12 = T(1.0540925533894598)
    TI21 = T(-3.3541019662496843)
    TI22 = T(0.3726779962499649)
    c1 = T2(1 / 3)
    c2 = T2(1.0)
    α = T(2.0)
    β = T(-sqrt(2))
    e1 = T(1 / 4)
    e2 = T(-1 / 4)
    RadauIIA3Tableau{T, T2}(T11, T12, T21, T22,
        TI11, TI12, TI21, TI22,
        c1, c2, α, β, e1, e2)
end

struct RadauIIA5Tableau{T, T2}
    T11::T
    T12::T
    T13::T
    T21::T
    T22::T
    T23::T
    T31::T
    #T32::T = 1
    #T33::T = 0
    TI11::T
    TI12::T
    TI13::T
    TI21::T
    TI22::T
    TI23::T
    TI31::T
    TI32::T
    TI33::T
    c1::T2
    c2::T2
    #c3::T2 = 1
    γ::T
    α::T
    β::T
    e1::T
    e2::T
    e3::T
end

# inv(T) * inv(A) * T = [γ  0  0
#                        0  α -β
#                        0  β  α]
function RadauIIA5Tableau(T, T2)
    T11 = convert(T, 9.1232394870892942792e-02)
    T12 = convert(T, -0.14125529502095420843e0)
    T13 = convert(T, -3.0029194105147424492e-02)
    T21 = convert(T, 0.24171793270710701896e0)
    T22 = convert(T, 0.20412935229379993199e0)
    T23 = convert(T, 0.38294211275726193779e0)
    T31 = convert(T, 0.96604818261509293619e0)
    TI11 = convert(T, 4.3255798900631553510e0)
    TI12 = convert(T, 0.33919925181580986954e0)
    TI13 = convert(T, 0.54177053993587487119e0)
    TI21 = convert(T, -4.1787185915519047273e0)
    TI22 = convert(T, -0.32768282076106238708e0)
    TI23 = convert(T, 0.47662355450055045196e0)
    TI31 = convert(T, -0.50287263494578687595e0)
    TI32 = convert(T, 2.5719269498556054292e0)
    TI33 = convert(T, -0.59603920482822492497e0)

    sqrt6 = sqrt(6)

    c1 = convert(T2, (4 - sqrt6) / 10)
    c2 = convert(T2, (4 + sqrt6) / 10)

    cbrt9 = cbrt(9)

    γ′ = convert(T, (6.0 + cbrt9 * (cbrt9 - 1)) / 30) # eigval of `A`
    α′ = convert(T, (12.0 - cbrt9 * (cbrt9 - 1)) / 60) # eigval of `A`
    β′ = convert(T, cbrt9 * (cbrt9 + 1) * sqrt(3) / 60) # eigval of `A`
    scale = α′^2 + β′^2
    γ = inv(γ′)  # eigval of `inv(A)`
    α = α′ / scale # eigval of `inv(A)`
    β = β′ / scale # eigval of `inv(A)`

    e1 = convert(T, -(13 + 7 * sqrt6) / 3)
    e2 = convert(T, (-13 + 7 * sqrt6) / 3)
    e3 = convert(T, -1 / 3)
    RadauIIA5Tableau{T, T2}(T11, T12, T13, T21, T22, T23, T31, #= T33 = 0 =#
        TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33,
        c1, c2, #= c3 = 1 =#
        γ, α, β,
        e1, e2, e3)
end