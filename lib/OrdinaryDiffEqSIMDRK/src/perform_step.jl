function initialize!(integrator,
        cache::Union{MER5v2ConstantCache, MER6v2ConstantCache, RK6v4ConstantCache})
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

using VectorizationBase, SLEEFPirates

#=
function foo(x, n = x ÷ 2 - 1)
    s = """
    c_$(x)_$(x+1) = Vec(c_$(x), c_$(x+1))
    a$(x)_$(x+1)_1 = Vec(a$(x)_1, a$(x+1)_1)
    """

    for i in 1:n
        s *= """
        a$(x)_$(2i)_$(x+1)_$(2i+1) = Vec(a$(x)_$(2i), a$(x+1)_$(2i+1))
        a$(x)_$(2i+1)_$(x+1)_$(2i) = Vec(a$(x)_$(2i+1), a$(x+1)_$(2i))
        """
    end
    s *= "k_$(x-1)_$(x-2) = VectorizationBase.shufflevector.(k_$(x-2)_$(x-1), Val{(1,0)}())"

    s *= "\nk_$(x)_$(x+1) = f(uprev +
    dt * (a$(x)_$(x+1)_1 * k_1"
    for i in 1:n
        s *= " + a$(x)_$(2i)_$(x+1)_$(2i+1) * k_$(2i)_$(2i+1)"
        s *= " + a$(x)_$(2i+1)_$(x+1)_$(2i) * k_$(2i+1)_$(2i)"
    end
    s *= "), p, t + c_$(x)_$(x+1) * dt)"
    return s
end
=#

# pirate
@inline function (f::ODEFunction)(
        v::AbstractArray{<:VectorizationBase.AbstractSIMD}, args...)
    @inline f.f(v, args...)
end

@muladd function perform_step!(
        integrator, cache::MER5v2ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (;  c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10, c_11, c_12, c_13, c_14, a2_1, a3_1,
        a4_1, a5_1, a6_1, a7_1, a8_1, a9_1, a10_1, a11_1, a12_1, a13_1,
        a14_1, a4_2, a5_2, a6_2, a7_2, a8_2, a9_2, a10_2, a11_2, a12_2,
        a13_2, a14_2, a4_3, a5_3, a6_3, a7_3, a8_3, a9_3, a10_3, a11_3,
        a12_3, a13_3, a14_3, a6_4, a7_4, a8_4, a9_4, a10_4, a11_4, a12_4,
        a13_4, a14_4, a6_5, a7_5, a8_5, a9_5, a10_5, a11_5, a12_5, a13_5,
        a14_5, a8_6, a9_6, a10_6, a11_6, a12_6, a13_6, a14_6, a8_7, a9_7,
        a10_7, a11_7, a12_7, a13_7, a14_7, a10_8, a11_8, a12_8, a13_8,
        a14_8, a10_9, a11_9, a12_9, a13_9, a14_9, a12_10, a13_10, a14_10,
        a12_11, a13_11, a14_11, a14_12, a14_13, btilde_1, btilde_2, btilde_3,
        btilde_4, btilde_5, btilde_6, btilde_7, btilde_8, btilde_9, btilde_10,
        btilde_11, btilde_12, btilde_13, btilde_14) = cache

    k_1 = integrator.fsalfirst

    a_2_3_1 = Vec(a2_1, a3_1)
    c_2_3 = Vec(c_2, c_3)

    k_2_3 = @inline f(uprev + dt * a_2_3_1 * k_1, p, t + c_2_3 * dt)
    #k_2 = f(uprev + dt * a2_1 * k_1, p, t + c_2 * dt)
    #k_3 = f(uprev + dt * a3_1 * k_1, p, t + c_3 * dt)

    k_3_2 = VectorizationBase.shufflevector.(k_2_3, Val{(1, 0)}())
    #Base.Cartesian.@nexprs 3 i -> a_4_5_i = Vec(a4_i, a5_i)
    a_4_5_1 = Vec(a4_1, a5_1)
    a_4_2_5_3 = Vec(a4_2, a5_3)
    a_4_3_5_2 = Vec(a4_3, a5_2)
    c_4_5 = Vec(c_4, c_5)
    k_4_5 = @inline f(
        uprev + dt * (a_4_5_1 * k_1 + a_4_2_5_3 * k_2_3 + a_4_3_5_2 * k_3_2), p,
        t + c_4_5 * dt)
    #k_4 = f(uprev + dt * (a4_1 * k_1 + a4_2 * k_2 + a4_3 * k_3), p, t + c_4 * dt)
    #k_5 = f(uprev + dt * (a5_1 * k_1 + a5_2 * k_2 + a5_3 * k_3), p, t + c_5 * dt)

    k_5_4 = VectorizationBase.shufflevector.(k_4_5, Val{(1, 0)}())
    a6_7_1 = Vec(a6_1, a7_1)
    a6_2_7_3 = Vec(a6_2, a7_3)
    a6_3_7_2 = Vec(a6_3, a7_2)
    a6_4_7_5 = Vec(a6_4, a7_5)
    a6_5_7_4 = Vec(a6_5, a7_4)
    c_6_7 = Vec(c_6, c_7)

    #k_2 = map(x -> x(1), k_2_3)
    #k_3 = map(x -> x(2), k_2_3)
    #k_4 = map(x -> x(1), k_4_5)
    #k_5 = map(x -> x(2), k_4_5)

    k_6_7 = @inline f(
        uprev +
        dt *
        (a6_7_1 * k_1 + a6_2_7_3 * k_2_3 + a6_3_7_2 * k_3_2 + a6_4_7_5 * k_4_5 +
         a6_5_7_4 * k_5_4),
        p,
        t + c_6_7 * dt)
    k_7_6 = VectorizationBase.shufflevector.(k_6_7, Val{(1, 0)}())

    c_8_9 = Vec(c_8, c_9)
    a8_9_1 = Vec(a8_1, a9_1)
    a8_2_9_3 = Vec(a8_2, a9_3)
    a8_3_9_2 = Vec(a8_3, a9_2)
    a8_4_9_5 = Vec(a8_4, a9_5)
    a8_5_9_4 = Vec(a8_5, a9_4)
    a8_6_9_7 = Vec(a8_6, a9_7)
    a8_7_9_6 = Vec(a8_7, a9_6)
    k_7_6 = VectorizationBase.shufflevector.(k_6_7, Val{(1, 0)}())
    k_8_9 = @inline f(
        uprev +
        dt *
        (a8_9_1 * k_1 + a8_2_9_3 * k_2_3 + a8_3_9_2 * k_3_2 + a8_4_9_5 * k_4_5 +
         a8_5_9_4 * k_5_4 + a8_6_9_7 * k_6_7 + a8_7_9_6 * k_7_6),
        p,
        t + c_8_9 * dt)

    c_10_11 = Vec(c_10, c_11)
    a10_11_1 = Vec(a10_1, a11_1)
    a10_2_11_3 = Vec(a10_2, a11_3)
    a10_3_11_2 = Vec(a10_3, a11_2)
    a10_4_11_5 = Vec(a10_4, a11_5)
    a10_5_11_4 = Vec(a10_5, a11_4)
    a10_6_11_7 = Vec(a10_6, a11_7)
    a10_7_11_6 = Vec(a10_7, a11_6)
    a10_8_11_9 = Vec(a10_8, a11_9)
    a10_9_11_8 = Vec(a10_9, a11_8)
    k_9_8 = VectorizationBase.shufflevector.(k_8_9, Val{(1, 0)}())
    k_10_11 = @inline f(
        uprev +
        dt * (a10_11_1 * k_1 + a10_2_11_3 * k_2_3 + a10_3_11_2 * k_3_2 +
         a10_4_11_5 * k_4_5 + a10_5_11_4 * k_5_4 + a10_6_11_7 * k_6_7 +
         a10_7_11_6 * k_7_6 + a10_8_11_9 * k_8_9 + a10_9_11_8 * k_9_8),
        p,
        t + c_10_11 * dt)
    c_12_13 = Vec(c_12, c_13)
    a12_13_1 = Vec(a12_1, a13_1)
    a12_2_13_3 = Vec(a12_2, a13_3)
    a12_3_13_2 = Vec(a12_3, a13_2)
    a12_4_13_5 = Vec(a12_4, a13_5)
    a12_5_13_4 = Vec(a12_5, a13_4)
    a12_6_13_7 = Vec(a12_6, a13_7)
    a12_7_13_6 = Vec(a12_7, a13_6)
    a12_8_13_9 = Vec(a12_8, a13_9)
    a12_9_13_8 = Vec(a12_9, a13_8)
    a12_10_13_11 = Vec(a12_10, a13_11)
    a12_11_13_10 = Vec(a12_11, a13_10)
    k_11_10 = VectorizationBase.shufflevector.(k_10_11, Val{(1, 0)}())
    k_12_13 = @inline f(
        uprev +
        dt * (a12_13_1 * k_1 + a12_2_13_3 * k_2_3 + a12_3_13_2 * k_3_2 +
         a12_4_13_5 * k_4_5 + a12_5_13_4 * k_5_4 + a12_6_13_7 * k_6_7 +
         a12_7_13_6 * k_7_6 + a12_8_13_9 * k_8_9 + a12_9_13_8 * k_9_8 +
         a12_10_13_11 * k_10_11 + a12_11_13_10 * k_11_10),
        p,
        t + c_12_13 * dt)

    u = uprev +
        dt * (a14_1 * k_1 +
         VectorizationBase.vsum.(Vec(a14_2, a14_3) * k_2_3 + Vec(a14_4, a14_5) * k_4_5 +
                                 Vec(a14_6, a14_7) * k_6_7 + Vec(a14_8, a14_9) * k_8_9 +
                                 Vec(a14_10, a14_11) * k_10_11 +
                                 Vec(a14_12, a14_13) * k_12_13))

    k_14 = integrator.fsallast = @inline f(u, p, t + dt)

    integrator.stats.nf += 12
    if integrator.opts.adaptive
        utilde = dt * (VectorizationBase.vsum.(Vec(btilde_2, btilde_3) * k_2_3 +
                                          Vec(btilde_4, btilde_5) * k_4_5 +
                                          Vec(btilde_6, btilde_7) * k_6_7 +
                                          Vec(btilde_8, btilde_9) * k_8_9 +
                                          Vec(btilde_10, btilde_11) * k_10_11 +
                                          Vec(btilde_12, btilde_13) * k_12_13 +
                                          Vec(btilde_1, btilde_14) * Vec.(k_1, k_14)))
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.u = u
end

@muladd function perform_step!(
        integrator, cache::MER6v2ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10,
        c_11, c_12, c_13, c_14, c_15, a2_1, a3_1,
        a4_1, a5_1, a6_1, a7_1, a8_1, a9_1, a10_1, a11_1, a12_1, a13_1,
        a14_1, a15_1, a4_2, a5_2, a6_2, a7_2, a8_2, a9_2, a10_2, a11_2, a12_2,
        a13_2, a14_2, a15_2, a4_3, a5_3, a6_3, a7_3, a8_3, a9_3, a10_3, a11_3,
        a12_3, a13_3, a14_3, a15_3, a6_4, a7_4, a8_4, a9_4, a10_4, a11_4, a12_4,
        a13_4, a14_4, a15_4, a6_5, a7_5, a8_5, a9_5, a10_5, a11_5, a12_5, a13_5,
        a14_5, a15_5, a8_6, a9_6, a10_6, a11_6, a12_6, a13_6, a14_6, a15_6, a8_7, a9_7,
        a10_7, a11_7, a12_7, a13_7, a14_7, a15_7, a10_8, a11_8, a12_8, a13_8,
        a14_8, a15_8, a10_9, a11_9, a12_9, a13_9, a14_9, a15_9, a12_10, a13_10, a14_10, a15_10,
        a12_11, a13_11, a14_11, a15_11, a14_12, a15_12, a14_13,
        a15_13, a15_14, btilde_1, btilde_2, btilde_3,
        btilde_4, btilde_5, btilde_6, btilde_7, btilde_8, btilde_9, btilde_10,
        btilde_11, btilde_12, btilde_13, btilde_14) = cache

    k_1 = integrator.fsalfirst

    a_2_3_1 = Vec(a2_1, a3_1)
    c_2_3 = Vec(c_2, c_3)

    k_2_3 = @inline f(uprev + dt * a_2_3_1 * k_1, p, t + c_2_3 * dt)

    k_3_2 = VectorizationBase.shufflevector.(k_2_3, Val{(1, 0)}())
    a_4_5_1 = Vec(a4_1, a5_1)
    a_4_2_5_3 = Vec(a4_2, a5_3)
    a_4_3_5_2 = Vec(a4_3, a5_2)
    c_4_5 = Vec(c_4, c_5)
    k_4_5 = @inline f(
        uprev + dt * (a_4_5_1 * k_1 + a_4_2_5_3 * k_2_3 + a_4_3_5_2 * k_3_2), p,
        t + c_4_5 * dt)

    k_5_4 = VectorizationBase.shufflevector.(k_4_5, Val{(1, 0)}())
    a6_7_1 = Vec(a6_1, a7_1)
    a6_2_7_3 = Vec(a6_2, a7_3)
    a6_3_7_2 = Vec(a6_3, a7_2)
    a6_4_7_5 = Vec(a6_4, a7_5)
    a6_5_7_4 = Vec(a6_5, a7_4)
    c_6_7 = Vec(c_6, c_7)
    k_6_7 = @inline f(
        uprev +
        dt *
        (a6_7_1 * k_1 + a6_2_7_3 * k_2_3 + a6_3_7_2 * k_3_2 + a6_4_7_5 * k_4_5 +
         a6_5_7_4 * k_5_4),
        p,
        t + c_6_7 * dt)
    k_7_6 = VectorizationBase.shufflevector.(k_6_7, Val{(1, 0)}())

    c_8_9 = Vec(c_8, c_9)
    a8_9_1 = Vec(a8_1, a9_1)
    a8_2_9_3 = Vec(a8_2, a9_3)
    a8_3_9_2 = Vec(a8_3, a9_2)
    a8_4_9_5 = Vec(a8_4, a9_5)
    a8_5_9_4 = Vec(a8_5, a9_4)
    a8_6_9_7 = Vec(a8_6, a9_7)
    a8_7_9_6 = Vec(a8_7, a9_6)
    k_7_6 = VectorizationBase.shufflevector.(k_6_7, Val{(1, 0)}())
    k_8_9 = @inline f(
        uprev +
        dt *
        (a8_9_1 * k_1 + a8_2_9_3 * k_2_3 + a8_3_9_2 * k_3_2 + a8_4_9_5 * k_4_5 +
         a8_5_9_4 * k_5_4 + a8_6_9_7 * k_6_7 + a8_7_9_6 * k_7_6),
        p,
        t + c_8_9 * dt)

    c_10_11 = Vec(c_10, c_11)
    a10_11_1 = Vec(a10_1, a11_1)
    a10_2_11_3 = Vec(a10_2, a11_3)
    a10_3_11_2 = Vec(a10_3, a11_2)
    a10_4_11_5 = Vec(a10_4, a11_5)
    a10_5_11_4 = Vec(a10_5, a11_4)
    a10_6_11_7 = Vec(a10_6, a11_7)
    a10_7_11_6 = Vec(a10_7, a11_6)
    a10_8_11_9 = Vec(a10_8, a11_9)
    a10_9_11_8 = Vec(a10_9, a11_8)
    k_9_8 = VectorizationBase.shufflevector.(k_8_9, Val{(1, 0)}())
    k_10_11 = @inline f(
        uprev +
        dt * (a10_11_1 * k_1 + a10_2_11_3 * k_2_3 + a10_3_11_2 * k_3_2 +
         a10_4_11_5 * k_4_5 + a10_5_11_4 * k_5_4 + a10_6_11_7 * k_6_7 +
         a10_7_11_6 * k_7_6 + a10_8_11_9 * k_8_9 + a10_9_11_8 * k_9_8),
        p,
        t + c_10_11 * dt)
    c_12_13 = Vec(c_12, c_13)
    a12_13_1 = Vec(a12_1, a13_1)
    a12_2_13_3 = Vec(a12_2, a13_3)
    a12_3_13_2 = Vec(a12_3, a13_2)
    a12_4_13_5 = Vec(a12_4, a13_5)
    a12_5_13_4 = Vec(a12_5, a13_4)
    a12_6_13_7 = Vec(a12_6, a13_7)
    a12_7_13_6 = Vec(a12_7, a13_6)
    a12_8_13_9 = Vec(a12_8, a13_9)
    a12_9_13_8 = Vec(a12_9, a13_8)
    a12_10_13_11 = Vec(a12_10, a13_11)
    a12_11_13_10 = Vec(a12_11, a13_10)
    k_11_10 = VectorizationBase.shufflevector.(k_10_11, Val{(1, 0)}())
    k_12_13 = @inline f(
        uprev +
        dt * (a12_13_1 * k_1 + a12_2_13_3 * k_2_3 + a12_3_13_2 * k_3_2 +
         a12_4_13_5 * k_4_5 + a12_5_13_4 * k_5_4 + a12_6_13_7 * k_6_7 +
         a12_7_13_6 * k_7_6 + a12_8_13_9 * k_8_9 + a12_9_13_8 * k_9_8 +
         a12_10_13_11 * k_10_11 + a12_11_13_10 * k_11_10),
        p,
        t + c_12_13 * dt)

    a14_2_14_3 = Vec(a14_2, a14_3)
    a14_4_14_5 = Vec(a14_4, a14_5)
    a14_6_14_7 = Vec(a14_6, a14_7)
    a14_8_14_9 = Vec(a14_8, a14_9)
    a14_10_14_11 = Vec(a14_10, a14_11)
    a14_12_14_13 = Vec(a14_12, a14_13)

    k_14 = @inline f(
        uprev +
        dt * (a14_1 * k_1 +
         VectorizationBase.vsum.(a14_2_14_3 * k_2_3 + a14_4_14_5 * k_4_5
                                 + a14_6_14_7 * k_6_7 + a14_8_14_9 * k_8_9
                                 + a14_10_14_11 * k_10_11 + a14_12_14_13 * k_12_13)),
        p,
        t + c_14 * dt)

    u = uprev +
        dt *
        (VectorizationBase.vsum.(Vec(a15_2, a15_3) * k_2_3 + Vec(a15_4, a15_5) * k_4_5 +
                                 Vec(a15_6, a15_7) * k_6_7 + Vec(a15_8, a15_9) * k_8_9 +
                                 Vec(a15_10, a15_11) * k_10_11 +
                                 Vec(a15_12, a15_13) * k_12_13 +
                                 Vec(a15_1, a15_14) * Vec.(k_1, k_14)))

    integrator.fsallast = @inline f(u, p, t + dt)

    integrator.stats.nf += 14
    if integrator.opts.adaptive
        utilde = dt * (VectorizationBase.vsum.(Vec(btilde_2, btilde_3) * k_2_3 +
                                          Vec(btilde_4, btilde_5) * k_4_5 +
                                          Vec(btilde_6, btilde_7) * k_6_7 +
                                          Vec(btilde_8, btilde_9) * k_8_9 +
                                          Vec(btilde_10, btilde_11) * k_10_11 +
                                          Vec(btilde_12, btilde_13) * k_12_13 +
                                          Vec(btilde_1, btilde_14) * Vec.(k_1, k_14)))
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.u = u
end

@muladd function perform_step!(
        integrator, cache::RK6v4ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10, c_11, c_12, c_13, c_14, c_15,
    c_16, c_17, c_18, c_19, c_20,
    c_21, c_22, a2_1, a3_1, a4_1, a5_1,
    a6_1, a7_1, a8_1, a9_1, a10_1, a11_1, a12_1,
    a13_1, a14_1, a15_1, a16_1, a17_1, a18_1, a19_1, a20_1,
    a21_1, a22_1, a6_2, a7_2, a8_2, a9_2, a10_2, a11_2, a12_2,
    a13_2, a14_2, a15_2, a16_2, a17_2, a18_2, a19_2,
    a20_2, a21_2, a22_2, a6_3, a7_3, a8_3, a9_3, a10_3, a11_3,
    a12_3, a13_3, a14_3, a15_3, a16_3, a17_3, a18_3,
    a19_3, a20_3, a21_3, a22_3, a6_4, a7_4, a8_4, a9_4, a10_4,
    a11_4, a12_4, a13_4, a14_4, a15_4, a16_4, a17_4,
    a18_4, a19_4, a20_4, a21_4,
    a22_4, a6_5, a7_5, a8_5, a9_5, a10_5, a11_5, a12_5, a13_5,
    a14_5, a15_5, a16_5, a17_5, a18_5, a19_5, a20_5,
    a21_5, a22_5, a10_6, a11_6, a12_6, a13_6, a14_6, a15_6,
    a16_6, a17_6, a18_6, a19_6, a20_6, a21_6, a22_6, a10_7,
    a11_7, a12_7, a13_7, a14_7, a15_7, a16_7, a17_7, a18_7,
    a19_7, a20_7, a21_7, a22_7, a10_8, a11_8, a12_8,
    a13_8, a14_8, a15_8, a16_8, a17_8, a18_8, a19_8, a20_8,
    a21_8, a22_8, a10_9, a11_9, a12_9, a13_9, a14_9, a15_9,
    a16_9, a17_9, a18_9, a19_9, a20_9, a21_9, a22_9, a14_10,
    a15_10, a16_10, a17_10, a18_10, a19_10, a20_10, a21_10,
    a22_10, a14_11, a15_11, a16_11, a17_11, a18_11, a19_11,
    a20_11, a21_11, a22_11, a14_12, a15_12, a16_12, a17_12,
    a18_12, a19_12, a20_12, a21_12, a22_12, a14_13, a15_13,
    a16_13, a17_13, a18_13, a19_13, a20_13, a21_13, a22_13,
    a18_14, a19_14, a20_14, a21_14, a22_14, a18_15, a19_15,
    a20_15, a21_15, a22_15, a18_16, a19_16, a20_16, a21_16,
    a22_16, a18_17, a19_17, a20_17, a21_17, a22_17, a22_18,
    a22_19, a22_20, a22_21, btilde_1, btilde_2, btilde_3,
    btilde_4, btilde_5, btilde_6, btilde_7, btilde_8, btilde_9,
    btilde_10, btilde_11, btilde_12, btilde_13,
    btilde_14, btilde_15, btilde_16, btilde_17, btilde_18, btilde_19, btilde_20, btilde_21) = cache

    k_1 = integrator.fsalfirst

    c_2_3_4_5 = Vec(c_2, c_3, c_4, c_5)
    a_2_3_4_5_1 = Vec(a2_1, a3_1, a4_1, a5_1)

    k_2_3_4_5 = @inline f(uprev + dt * (a_2_3_4_5_1 * k_1), p, t +
                                                               c_2_3_4_5 * dt)

    c_6_7_8_9 = Vec(c_6, c_7, c_8, c_9)
    a_6_7_8_9_1 = Vec(a6_1, a7_1, a8_1, a9_1)
    a_6_2_7_3_8_4_9_5 = Vec(a6_2, a7_3, a8_4, a9_5)
    a_6_3_7_4_8_5_9_2 = Vec(a6_3, a7_4, a8_5, a9_2)
    a_6_4_7_5_8_2_9_3 = Vec(a6_4, a7_5, a8_2, a9_3)
    a_6_5_7_2_8_3_9_4 = Vec(a6_5, a7_2, a8_3, a9_4)

    k_3_4_5_2 = VectorizationBase.shufflevector.(k_2_3_4_5, Val{(1, 2, 3, 0)}())
    k_4_5_2_3 = VectorizationBase.shufflevector.(k_2_3_4_5, Val{(2, 3, 0, 1)}())
    k_5_2_3_4 = VectorizationBase.shufflevector.(k_2_3_4_5, Val{(3, 0, 1, 2)}())
    k_6_7_8_9 = @inline f(
        uprev +
        dt *
        (a_6_7_8_9_1 * k_1 + a_6_2_7_3_8_4_9_5 * k_2_3_4_5 + a_6_3_7_4_8_5_9_2 * k_3_4_5_2 +
         a_6_4_7_5_8_2_9_3 * k_4_5_2_3 + a_6_5_7_2_8_3_9_4 * k_5_2_3_4),
        p,
        t +
        c_6_7_8_9 * dt)

    c_10_11_12_13 = Vec(c_10, c_11, c_12, c_13)
    a_10_11_12_13_1 = Vec(a10_1, a11_1, a12_1, a13_1)
    a_10_2_11_3_12_4_13_5 = Vec(a10_2, a11_3, a12_4, a13_5)
    a_10_3_11_4_12_5_13_2 = Vec(a10_3, a11_4, a12_5, a13_2)
    a_10_4_11_5_12_2_13_3 = Vec(a10_4, a11_5, a12_2, a13_3)
    a_10_5_11_2_12_3_13_4 = Vec(a10_5, a11_2, a12_3, a13_4)
    a_10_6_11_7_12_8_13_9 = Vec(a10_6, a11_7, a12_8, a13_9)
    a_10_7_11_8_12_9_13_6 = Vec(a10_7, a11_8, a12_9, a13_6)
    a_10_8_11_9_12_6_13_7 = Vec(a10_8, a11_9, a12_6, a13_7)
    a_10_9_11_6_12_7_13_8 = Vec(a10_9, a11_6, a12_7, a13_8)

    k_7_8_9_6 = VectorizationBase.shufflevector.(k_6_7_8_9, Val{(1, 2, 3, 0)}())
    k_8_9_6_7 = VectorizationBase.shufflevector.(k_6_7_8_9, Val{(2, 3, 0, 1)}())
    k_9_6_7_8 = VectorizationBase.shufflevector.(k_6_7_8_9, Val{(3, 0, 1, 2)}())
    k_10_11_12_13 = f(
        uprev +
        dt * (a_10_11_12_13_1 * k_1 + a_10_2_11_3_12_4_13_5 * k_2_3_4_5 +
         a_10_3_11_4_12_5_13_2 * k_3_4_5_2 + a_10_4_11_5_12_2_13_3 * k_4_5_2_3 +
         a_10_5_11_2_12_3_13_4 * k_5_2_3_4 + a_10_6_11_7_12_8_13_9 * k_6_7_8_9 +
         a_10_7_11_8_12_9_13_6 * k_7_8_9_6 + a_10_8_11_9_12_6_13_7 * k_8_9_6_7 +
         a_10_9_11_6_12_7_13_8 * k_9_6_7_8),
        p,
        t +
        c_10_11_12_13 * dt)

    c_14_15_16_17 = Vec(c_14, c_15, c_16, c_17)
    a_14_15_16_17_1 = Vec(a14_1, a15_1, a16_1, a17_1)
    a_14_2_15_3_16_4_17_5 = Vec(a14_2, a15_3, a16_4, a17_5)
    a_14_3_15_4_16_5_17_2 = Vec(a14_3, a15_4, a16_5, a17_2)
    a_14_4_15_5_16_2_17_3 = Vec(a14_4, a15_5, a16_2, a17_3)
    a_14_5_15_2_16_3_17_4 = Vec(a14_5, a15_2, a16_3, a17_4)
    a_14_6_15_7_16_8_17_9 = Vec(a14_6, a15_7, a16_8, a17_9)
    a_14_7_15_8_16_9_17_6 = Vec(a14_7, a15_8, a16_9, a17_6)
    a_14_8_15_9_16_6_17_7 = Vec(a14_8, a15_9, a16_6, a17_7)
    a_14_9_15_6_16_7_17_8 = Vec(a14_9, a15_6, a16_7, a17_8)
    a_14_10_15_11_16_12_17_13 = Vec(a14_10, a15_11, a16_12, a17_13)
    a_14_11_15_12_16_13_17_10 = Vec(a14_11, a15_12, a16_13, a17_10)
    a_14_12_15_13_16_10_17_11 = Vec(a14_12, a15_13, a16_10, a17_11)
    a_14_13_15_10_16_11_17_12 = Vec(a14_13, a15_10, a16_11, a17_12)

    k_11_12_13_10 = VectorizationBase.shufflevector.(k_10_11_12_13, Val{(1, 2, 3, 0)}())
    k_12_13_10_11 = VectorizationBase.shufflevector.(k_10_11_12_13, Val{(2, 3, 0, 1)}())
    k_13_10_11_12 = VectorizationBase.shufflevector.(k_10_11_12_13, Val{(3, 0, 1, 2)}())
    k_14_15_16_17 = f(
        uprev +
        dt * (a_14_15_16_17_1 * k_1 + a_14_2_15_3_16_4_17_5 * k_2_3_4_5 +
         a_14_3_15_4_16_5_17_2 * k_3_4_5_2 + a_14_4_15_5_16_2_17_3 * k_4_5_2_3 +
         a_14_5_15_2_16_3_17_4 * k_5_2_3_4 + a_14_6_15_7_16_8_17_9 * k_6_7_8_9 +
         a_14_7_15_8_16_9_17_6 * k_7_8_9_6 + a_14_8_15_9_16_6_17_7 * k_8_9_6_7 +
         a_14_9_15_6_16_7_17_8 * k_9_6_7_8 + a_14_10_15_11_16_12_17_13 * k_10_11_12_13 +
         a_14_11_15_12_16_13_17_10 * k_11_12_13_10 +
         a_14_12_15_13_16_10_17_11 * k_12_13_10_11 +
         a_14_13_15_10_16_11_17_12 * k_13_10_11_12),
        p,
        t +
        c_14_15_16_17 * dt)

    c_18_19_20_21 = Vec(c_18, c_19, c_20, c_21)
    a_18_19_20_21_1 = Vec(a18_1, a19_1, a20_1, a21_1)
    a_18_2_19_3_20_4_21_5 = Vec(a18_2, a19_3, a20_4, a21_5)
    a_18_3_19_4_20_5_21_2 = Vec(a18_3, a19_4, a20_5, a21_2)
    a_18_4_19_5_20_2_21_3 = Vec(a18_4, a19_5, a20_2, a21_3)
    a_18_5_19_2_20_3_21_4 = Vec(a18_5, a19_2, a20_3, a21_4)
    a_18_6_19_7_20_8_21_9 = Vec(a18_6, a19_7, a20_8, a21_9)
    a_18_7_19_8_20_9_21_6 = Vec(a18_7, a19_8, a20_9, a21_6)
    a_18_8_19_9_20_6_21_7 = Vec(a18_8, a19_9, a20_6, a21_7)
    a_18_9_19_6_20_7_21_8 = Vec(a18_9, a19_6, a20_7, a21_8)
    a_18_10_19_11_20_12_21_13 = Vec(a18_10, a19_11, a20_12, a21_13)
    a_18_11_19_12_20_13_21_10 = Vec(a18_11, a19_12, a20_13, a21_10)
    a_18_12_19_13_20_10_21_11 = Vec(a18_12, a19_13, a20_10, a21_11)
    a_18_13_19_10_20_11_21_12 = Vec(a18_13, a19_10, a20_11, a21_12)
    a_18_14_19_15_20_16_21_17 = Vec(a18_14, a19_15, a20_16, a21_17)
    a_18_15_19_16_20_17_21_14 = Vec(a18_15, a19_16, a20_17, a21_14)
    a_18_16_19_17_20_14_21_15 = Vec(a18_16, a19_17, a20_14, a21_15)
    a_18_17_19_14_20_15_21_16 = Vec(a18_17, a19_14, a20_15, a21_16)

    k_15_16_17_14 = VectorizationBase.shufflevector.(k_14_15_16_17, Val{(1, 2, 3, 0)}())
    k_16_17_14_15 = VectorizationBase.shufflevector.(k_14_15_16_17, Val{(2, 3, 0, 1)}())
    k_17_14_15_16 = VectorizationBase.shufflevector.(k_14_15_16_17, Val{(3, 0, 1, 2)}())
    k_18_19_20_21 = f(
        uprev +
        dt * (a_18_19_20_21_1 * k_1 + a_18_2_19_3_20_4_21_5 * k_2_3_4_5 +
         a_18_3_19_4_20_5_21_2 * k_3_4_5_2 + a_18_4_19_5_20_2_21_3 * k_4_5_2_3 +
         a_18_5_19_2_20_3_21_4 * k_5_2_3_4 + a_18_6_19_7_20_8_21_9 * k_6_7_8_9 +
         a_18_7_19_8_20_9_21_6 * k_7_8_9_6 + a_18_8_19_9_20_6_21_7 * k_8_9_6_7 +
         a_18_9_19_6_20_7_21_8 * k_9_6_7_8 + a_18_10_19_11_20_12_21_13 * k_10_11_12_13 +
         a_18_11_19_12_20_13_21_10 * k_11_12_13_10 +
         a_18_12_19_13_20_10_21_11 * k_12_13_10_11 +
         a_18_13_19_10_20_11_21_12 * k_13_10_11_12 +
         a_18_14_19_15_20_16_21_17 * k_14_15_16_17 +
         a_18_15_19_16_20_17_21_14 * k_15_16_17_14 +
         a_18_16_19_17_20_14_21_15 * k_16_17_14_15 +
         a_18_17_19_14_20_15_21_16 * k_17_14_15_16),
        p,
        t +
        c_18_19_20_21 * dt)

    u = uprev +
        dt * (a22_1 * k_1 +
         VectorizationBase.vsum.(Vec(a22_2, a22_3, a22_4, a22_5) * k_2_3_4_5 +
                                 Vec(a22_6, a22_7, a22_8, a22_9) * k_6_7_8_9 +
                                 Vec(a22_10, a22_11, a22_12, a22_13) * k_10_11_12_13 +
                                 Vec(a22_14, a22_15, a22_16, a22_17) * k_14_15_16_17 +
                                 Vec(a22_18, a22_19, a22_20, a22_21) * k_18_19_20_21))

    integrator.fsallast = f(u, p, t + dt)

    integrator.stats.nf += 21
    if integrator.opts.adaptive
        utilde = dt * (btilde_1 * k_1 +
                  VectorizationBase.vsum.(Vec(btilde_2, btilde_3, btilde_4, btilde_5) *
                                          k_2_3_4_5 +
                                          Vec(btilde_6, btilde_7, btilde_8, btilde_9) *
                                          k_6_7_8_9 +
                                          Vec(btilde_10, btilde_11, btilde_12, btilde_13) *
                                          k_10_11_12_13 +
                                          Vec(btilde_14, btilde_15, btilde_16, btilde_17) *
                                          k_14_15_16_17 +
                                          Vec(btilde_18, btilde_19, btilde_20, btilde_21) *
                                          k_18_19_20_21))
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.u = u
end
