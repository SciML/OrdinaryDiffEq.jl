function initialize!(integrator, cache::MagnusMidpointCache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusMidpointCache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, u, p, t + dt / 2)

    if alg.krylov
        u .= expv(dt, L, u; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        A = convert(AbstractMatrix, L)
        u .= exponential!(dt * A, exp_method, exp_cache) * u
    end

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end
function initialize!(integrator, cache::LieRK4Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LieRK4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f

    Y1 = uprev
    update_coefficients!(L, Y1, p, t)
    k1 = dt * convert(AbstractMatrix, L)

    Y2 = exponential!(k1 / 2, exp_method, exp_cache) * uprev
    update_coefficients!(L, Y2, p, t)
    k2 = dt * convert(AbstractMatrix, L)

    Y3 = exponential!(k2 / 2, exp_method, exp_cache) * uprev
    update_coefficients!(L, Y3, p, t)
    k3 = dt * convert(AbstractMatrix, L)

    Y4 = exponential!(k3 - k1 / 2, exp_method, exp_cache) * Y2
    update_coefficients!(L, Y4, p, t)
    k4 = dt * convert(AbstractMatrix, L)

    y1_2 = exponential!((3 * k1 + 2 * k2 + 2 * k3 - k4) / 12, exp_method, exp_cache) * uprev

    if alg.krylov
        u .= expv((1 / 12), (-k1 + 2 * k2 + 2 * k3 + 3 * k4), y1_2;
            m = min(alg.m, size(L, 1)), opnorm = integrator.opts.internalopnorm,
            iop = alg.iop)
    else
        u .= exponential!((1 / 12) * (-k1 + 2 * k2 + 2 * k3 + 3 * k4),
            exp_method,
            exp_cache) *
             y1_2
    end

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::RKMK4Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::RKMK4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    k1 = dt * convert(AbstractMatrix, L)
    update_coefficients!(L, exponential!(k1 / 2, exp_method, exp_cache) * uprev, p, t)
    k2 = dt * convert(AbstractMatrix, L)
    update_coefficients!(L,
        exponential!(k1 / 2 - (k1 * k2 - k2 * k1) / 8, exp_method, exp_cache) * uprev,
        p,
        t)
    k3 = dt * convert(AbstractMatrix, L)
    update_coefficients!(L, exponential!(k3, exp_method, exp_cache) * uprev, p, t)
    k4 = dt * convert(AbstractMatrix, L)
    if alg.krylov
        u .= expv(1 / 6, (k1 + 2 * k2 + 2 * k3 + k4 - (k1 * k4 - k4 * k1) / 2), uprev;
            m = min(alg.m, size(L, 1)), opnorm = integrator.opts.internalopnorm,
            iop = alg.iop)
    else
        u .= exponential!((1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4 - (k1 * k4 - k4 * k1) / 2),
            exp_method,
            exp_cache) * uprev
    end

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::RKMK2Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::RKMK2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    k1 = dt * convert(AbstractMatrix, L)
    update_coefficients!(L, exponential!(k1, exp_method, exp_cache) * uprev, p, t)
    k2 = dt * convert(AbstractMatrix, L)
    if alg.krylov
        u .= expv(1 / 2, (k1 + k2), uprev; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!((1 / 2) * (k1 + k2), exp_method, exp_cache) * uprev
    end

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::CG3Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::CG3Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    A = deepcopy(convert(AbstractMatrix, L))
    v2 = exponential!((3 / 4) * dt * A, exp_method, exp_cache) * uprev
    update_coefficients!(L, v2, p, t + (3 * dt / 4))
    B = deepcopy(convert(AbstractMatrix, L))
    v3 = exponential!((119 / 216) * dt * B, exp_method, exp_cache) *
         exponential!((17 / 108) * dt * A, exp_method, exp_cache) * uprev
    update_coefficients!(L, v3, p, t + (17 * dt / 24))
    C = convert(AbstractMatrix, L)
    u .= (exponential!(dt * (24 / 17) * C, exp_method, exp_cache) *
          exponential!(dt * (-2 / 3) * B, exp_method, exp_cache) *
          exponential!(dt * (13 / 51) * A, exp_method, exp_cache)) *
         uprev

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::CG2Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::CG2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    k1 = dt * convert(AbstractMatrix, L)
    update_coefficients!(L, exponential!(k1, exp_method, exp_cache) * uprev, p, t)
    k2 = dt * convert(AbstractMatrix, L)
    u .= (exponential!((1 / 2) * (k1), exp_method, exp_cache) *
          (exponential!((1 / 2) * k2, exp_method, exp_cache))) * uprev

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::CG4aCache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::CG4aCache, repeat_step = false)
    @unpack t, dt, uprev, u, p, alg = integrator
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    A = deepcopy(convert(AbstractMatrix, L))
    v2 = exponential!((0.8177227988124852) * dt * A, exp_method, exp_cache) * uprev
    update_coefficients!(L, v2, p, t + (0.8177227988124852 * dt))
    B = deepcopy(convert(AbstractMatrix, L))
    v3 = exponential!((0.3199876375476427) * dt * B, exp_method, exp_cache) *
         exponential!((0.0659864263556022) * dt * A, exp_method, exp_cache) * uprev
    update_coefficients!(L, v3, p, t + (0.3859740639032449 * dt))
    C = deepcopy(convert(AbstractMatrix, L))
    v4 = exponential!((0.9214417194464946) * dt * C, exp_method, exp_cache) *
         exponential!((0.4997857776773573) * dt * B, exp_method, exp_cache) *
         exponential!((-1.0969984448371582) * dt * A, exp_method, exp_cache) * uprev
    update_coefficients!(L, v4, p, t + (0.3242290522866937 * dt))
    D = deepcopy(convert(AbstractMatrix, L))
    v5 = exponential!((0.3552358559023322) * dt * D, exp_method, exp_cache) *
         exponential!((0.2390958372307326) * dt * C, exp_method, exp_cache) *
         exponential!((1.3918565724203246) * dt * B, exp_method, exp_cache) *
         exponential!((-1.1092979392113465) * dt * A, exp_method, exp_cache) * uprev
    update_coefficients!(L, v5, p, t + (0.8768903263420429 * dt))
    E = convert(AbstractMatrix, L)
    u .= (exponential!(dt * (0.3322195591068374) * E, exp_method, exp_cache) *
          exponential!(dt * (-0.1907142565505889) * D, exp_method, exp_cache) *
          exponential!(dt * (0.7397813985370780) * C, exp_method, exp_cache) *
          exponential!(dt * (-0.0183698531564020) * B, exp_method, exp_cache) *
          exponential!(dt * (0.1370831520630755) * A, exp_method, exp_cache)) * uprev

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusAdapt4Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusAdapt4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack W, k, tmp, utilde, atmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    A0 = deepcopy(convert(AbstractMatrix, L))
    k1 = dt * A0
    Q1 = k1

    y2 = (1 / 2) * Q1
    update_coefficients!(L, exponential!(y2, exp_method, exp_cache) * uprev, p, t + dt / 2)
    A1 = deepcopy(convert(AbstractMatrix, L))
    k2 = dt * A1
    Q2 = k2 - k1

    y3 = (1 / 2) * Q1 + (1 / 4) * Q2
    update_coefficients!(L, exponential!(y3, exp_method, exp_cache) * uprev, p, t + dt / 2)
    A2 = deepcopy(convert(AbstractMatrix, L))
    k3 = dt * A2
    Q3 = k3 - k2

    y4 = Q1 + Q2
    update_coefficients!(L, exponential!(y4, exp_method, exp_cache) * uprev, p, t + dt)
    A3 = deepcopy(convert(AbstractMatrix, L))
    k4 = dt * A3
    Q4 = k4 - 2 * k2 + k1

    y5 = (1 / 2) * Q1 + (1 / 4) * Q2 + (1 / 3) * Q3 - (1 / 24) * Q4 -
         (1 / 48) * (Q1 * Q2 - Q2 * Q1)
    update_coefficients!(L, exponential!(y5, exp_method, exp_cache) * uprev, p, t + dt / 2)
    A4 = deepcopy(convert(AbstractMatrix, L))
    k5 = dt * A4
    Q5 = k5 - k2

    y6 = Q1 + Q2 + (2 / 3) * Q3 + (1 / 6) * Q4 - (1 / 6) * (Q1 * Q2 - Q2 * Q1)
    update_coefficients!(L, exponential!(y6, exp_method, exp_cache) * uprev, p, t + dt)
    A5 = convert(AbstractMatrix, L)
    k6 = dt * A5
    Q6 = k6 - 2 * k2 + k1

    v4 = Q1 + Q2 + (2 / 3) * Q5 + (1 / 6) * Q6 -
         (1 / 6) * (Q1 * (Q2 - Q3 + Q5 + (1 / 2) * Q6) - (Q2 - Q3 + Q5 + (1 / 2) * Q6) * Q1)

    u .= exponential!(v4, exp_method, exp_cache) * uprev

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    if integrator.opts.adaptive
        utilde = u - exponential!(y6, exp_method, exp_cache) * uprev
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::MagnusNC8Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusNC8Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f

    update_coefficients!(L, u, p, t)
    A0 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, u, p, t + dt / 6)
    A1 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, u, p, t + 2 * dt / 6)
    A2 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, u, p, t + 3 * dt / 6)
    A3 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, u, p, t + 4 * dt / 6)
    A4 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, u, p, t + 5 * dt / 6)
    A5 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, u, p, t + dt)
    A6 = convert(AbstractMatrix, L)

    S1 = A0 + A6
    S2 = A1 + A5
    S3 = A2 + A4
    S4 = A3
    R1 = A6 - A0
    R2 = A5 - A1
    R3 = A4 - A2

    B0 = (1 / 840) * (41 * S1 + 216 * S2 + 27 * S3 + 272 * S4)
    B2 = (1 / 840) * (41 * S1 / 4 + 216 * S2 / 9 + 27 * S3 / 36)
    B1 = (1 / 840) * (41 * R1 / 2 + 216 * R2 / 3 + 27 * R3 / 6)
    B3 = (1 / 840) * (41 * R1 / 8 + 216 * R2 / 27 + 27 * R3 / 216)

    Q1 = (-38 * B0 / 5 + 24 * B2) * B3 - B3 * (-38 * B0 / 5 + 24 * B2)
    Q2 = (63 * B0 / 5 - 84 * B2) * (-5 * B1 / 28 + B3) -
         (-5 * B1 / 28 + B3) * (63 * B0 / 5 - 84 * B2)
    Q3 = (19 * B0 / 28 - 15 * B2 / 7) * (B0 * (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) -
          (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) * B0) -
         (B0 * (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) -
          (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) * B0) * (19 * B0 / 28 - 15 * B2 / 7)
    Q4 = B3 * (20 * Q1 / 7 + 10 * Q2) - (20 * Q1 / 7 + 10 * Q2) * B3
    Q5 = (-6025 * B0 / 4116 + 2875 * B2 / 343) * (B2 * Q1 - Q1 * B2) -
         (B2 * Q1 - Q1 * B2) * (-6025 * B0 / 4116 + 2875 * B2 / 343)
    Q6 = B3 * (20 * (Q3 + Q4) / 7 + 820 * dt * Q5 / 189) -
         (20 * (Q3 + Q4) / 7 + 820 * dt * Q5 / 189) * B3
    Q7 = (-1 / 42) * (B0 * (B0 * (Q3 - Q4 / 3 + dt * Q5) - (Q3 - Q4 / 3 + dt * Q5) * B0) -
          (B0 * (Q3 - Q4 / 3 + dt * Q5) - (Q3 - Q4 / 3 + dt * Q5) * B0) * B0)

    Ω1 = dt * B0
    Ω2 = (dt^2) * (Q1 + Q2)
    Ω3_4_5_6 = (dt^3) * (Q3 + Q4) + (dt^4) * (Q5 + Q6) + (dt^5) * Q7
    if alg.krylov
        u .= expv(1.0, Ω1 + Ω2 + Ω3_4_5_6, uprev; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!(Ω1 + Ω2 + Ω3_4_5_6, exp_method, exp_cache) * uprev
    end
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusGL4Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusGL4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t + dt * (1 / 2 - sqrt(3) / 6))
    A1 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt * (1 / 2 + sqrt(3) / 6))
    A2 = convert(AbstractMatrix, L)

    Ω = (dt / 2) * (A1 + A2) - (dt^2) * (sqrt(3) / 12) * (A1 * A2 - A2 * A1)

    if alg.krylov
        u .= expv(1.0, Ω, uprev; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!(Ω, exp_method, exp_cache) * uprev
    end
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusGL8Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusGL8Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    v1 = (1 / 2) * sqrt((3 + 2 * sqrt(6 / 5)) / 7)
    v2 = (1 / 2) * sqrt((3 - 2 * sqrt(6 / 5)) / 7)
    update_coefficients!(L, uprev, p, t - dt * v1 + dt / 2)
    A1 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt * v1 + dt / 2)
    A4 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t - dt * v2 + dt / 2)
    A2 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt * v2 + dt / 2)
    A3 = convert(AbstractMatrix, L)
    w1 = (1 / 2) - (1 / 6) * sqrt(5 / 6)
    w2 = (1 / 2) + (1 / 6) * sqrt(5 / 6)
    S1 = A1 + A4
    S2 = A2 + A3
    R1 = A4 - A1
    R2 = A3 - A2

    B0 = (1 / 2) * (w1 * S1 + w2 * S2)
    B2 = (1 / 2) * ((v1^2) * w1 * S1 + (v2^2) * w2 * S2)
    B1 = (1 / 2) * (v1 * w1 * R1 + v2 * w2 * R2)
    B3 = (1 / 2) * ((v1^3) * w1 * R1 + (v2^3) * w2 * R2)
    Q1 = (-38 * B0 / 5 + 24 * B2) * B3 - B3 * (-38 * B0 / 5 + 24 * B2)
    Q2 = (63 * B0 / 5 - 84 * B2) * (-5 * B1 / 28 + B3) -
         (-5 * B1 / 28 + B3) * (63 * B0 / 5 - 84 * B2)
    Q3 = (19 * B0 / 28 - 15 * B2 / 7) * (B0 * (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) -
          (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) * B0) -
         (B0 * (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) -
          (B2 + dt * (61 * Q1 / 588 - Q2 / 12)) * B0) * (19 * B0 / 28 - 15 * B2 / 7)
    Q4 = B3 * (20 * Q1 / 7 + 10 * Q2) - (20 * Q1 / 7 + 10 * Q2) * B3
    Q5 = (-6025 * B0 / 4116 + 2875 * B2 / 343) * (B2 * Q1 - Q1 * B2) -
         (B2 * Q1 - Q1 * B2) * (-6025 * B0 / 4116 + 2875 * B2 / 343)
    Q6 = B3 * (20 * (Q3 + Q4) / 7 + 820 * dt * Q5 / 189) -
         (20 * (Q3 + Q4) / 7 + 820 * dt * Q5 / 189) * B3
    Q7 = (-1 / 42) * (B0 * (B0 * (Q3 - Q4 / 3 + dt * Q5) - (Q3 - Q4 / 3 + dt * Q5) * B0) -
          (B0 * (Q3 - Q4 / 3 + dt * Q5) - (Q3 - Q4 / 3 + dt * Q5) * B0) * B0)

    Ω1 = dt * B0
    Ω2 = (dt^2) * (Q1 + Q2)
    Ω3_4_5_6 = (dt^3) * (Q3 + Q4) + (dt^4) * (Q5 + Q6) + (dt^5) * Q7
    if alg.krylov
        u .= expv(1.0, Ω1 + Ω2 + Ω3_4_5_6, uprev; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!(Ω1 + Ω2 + Ω3_4_5_6, exp_method, exp_cache) * uprev
    end
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusNC6Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusNC6Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t)
    A0 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt / 4)
    A1 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt / 2)
    A2 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + 3 * dt / 4)
    A3 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt)
    A4 = convert(AbstractMatrix, L)
    B0 = (1 / 90) * (7 * (A0 + A4) + 32 * (A1 + A3) + 12 * (A2))
    B1 = (1 / 90) * ((7 / 2) * (A4 - A0) + 8 * (A3 - A1))
    B2 = (1 / 90) * ((7 / 4) * (A0 + A4) + 2 * (A1 + A3))
    Ω1 = dt * B0
    Ω2 = (dt * dt) * (B1 * (3 * B0 / 2 - 6 * B2) - (3 * B0 / 2 - 6 * B2) * B1)
    Ω3_4 = (dt * dt) *
           (B0 * (B0 * (dt * B2 / 2 - Ω2 / 60) - (dt * B2 / 2 - Ω2 / 60) * B0) -
            (B0 * (dt * B2 / 2 - Ω2 / 60) - (dt * B2 / 2 - Ω2 / 60) * B0) * B0) +
           (3 * dt / 5) * (B1 * Ω2 - Ω2 * B1)
    if alg.krylov
        u .= expv(1.0, Ω1 + Ω2 + Ω3_4, uprev; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!(Ω1 + Ω2 + Ω3_4, exp_method, exp_cache) * uprev
    end
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusGL6Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusGL6Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t - dt * (sqrt(3 / 20)) + dt / 2)
    A1 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt / 2)
    A2 = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt * (sqrt(3 / 20)) + dt / 2)
    A3 = convert(AbstractMatrix, L)
    B0 = (1 / 18) * (5 * (A1 + A3) + 8 * A2)
    B1 = (sqrt(15) / 36) * (A3 - A1)
    B2 = (1 / 24) * (A1 + A3)
    Ω1 = dt * B0
    Ω2 = (dt * dt) * (B1 * (3 * B0 / 2 - 6 * B2) - (3 * B0 / 2 - 6 * B2) * B1)
    Ω3_4 = (dt * dt) *
           (B0 * (B0 * (dt * B2 / 2 - Ω2 / 60) - (dt * B2 / 2 - Ω2 / 60) * B0) -
            (B0 * (dt * B2 / 2 - Ω2 / 60) - (dt * B2 / 2 - Ω2 / 60) * B0) * B0) +
           (3 * dt / 5) * (B1 * Ω2 - Ω2 * B1)
    if alg.krylov
        u .= expv(1.0, Ω1 + Ω2 + Ω3_4, uprev; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!(Ω1 + Ω2 + Ω3_4, exp_method, exp_cache) * uprev
    end
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusGauss4Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusGauss4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, uprev, p, t + dt * (1 / 2 + sqrt(3) / 6))
    A = deepcopy(convert(AbstractMatrix, L))
    update_coefficients!(L, uprev, p, t + dt * (1 / 2 - sqrt(3) / 6))
    B = convert(AbstractMatrix, L)
    if alg.krylov
        u .= expv(dt, (A + B) ./ 2 + (dt * sqrt(3)) .* (B * A - A * B) ./ 12, u;
            m = min(alg.m, size(L, 1)), opnorm = integrator.opts.internalopnorm,
            iop = alg.iop)
    else
        u .= exponential!(
            (dt / 2) .* (A + B) +
            ((dt^2) * (sqrt(3) / 12)) .* (B * A - A * B),
            exp_method,
            exp_cache) * uprev
    end
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::LieEulerCache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LieEulerCache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()

    L = integrator.f.f
    update_coefficients!(L, u, p, t)

    if alg.krylov
        u .= expv(dt, L, u; m = min(alg.m, size(L, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u .= exponential!(dt * convert(AbstractMatrix, L), exp_method, exp_cache) * u
    end

    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::MagnusLeapfrogCache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MagnusLeapfrogCache, repeat_step = false,
        alg_extrapolates = true, iter = 1)
    @unpack t, dt, uprev, uprev2, u, p, iter = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack W, k, tmp, exp_cache = cache
    mass_matrix = integrator.f.mass_matrix
    exp_method = ExpMethodGeneric()
    # println("iter   : $iter")
    if iter == 1
        L = integrator.f.f
        update_coefficients!(L, u, p, t + dt / 2)
        if alg.krylov
            u .= expv(dt, L, u; m = min(alg.m, size(L, 1)),
                opnorm = integrator.opts.internalopnorm, iop = alg.iop)
        else
            A = convert(AbstractMatrix, L)
            u .= exponential!(dt * A, exp_method, exp_cache) * u
        end

        integrator.f(integrator.fsallast, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        iter += 1
    else
        L = integrator.f.f
        update_coefficients!(L, u, p, t)
        if alg.krylov
            u .= expv(2 * dt, L, uprev2; m = min(alg.m, size(L, 1)),
                opnorm = integrator.opts.internalopnorm, iop = alg.iop)
        else
            A = convert(AbstractMatrix, L)
            u .= exponential!(2 * dt * A, exp_method, exp_cache) * uprev2
        end
        uprev = u
        integrator.f(integrator.fsallast, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
end

function initialize!(integrator, cache::LinearExponentialConstantCache)
    # Pre-start fsal
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)

    # Initialize interpolation derivatives
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::LinearExponentialConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, f, p = integrator
    alg = unwrap_alg(integrator, nothing)
    A = convert(AbstractMatrix, f.f) # assume f to be an ODEFunction wrapped around a linear operator

    if alg.krylov == :off
        u = exponential!(dt * A, ExpMethodGeneric()) * integrator.u
    elseif alg.krylov == :simple
        u = expv(dt, A, integrator.u; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
    else
        u = expv_timestep(dt, A, integrator.u; m = min(alg.m, size(A, 1)), iop = alg.iop,
            opnorm = integrator.opts.internalopnorm,
            tol = integrator.opts.reltol)
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::LinearExponentialCache)
    # Pre-start fsal
    integrator.fsalfirst = zero(cache.rtmp)
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)

    # Initialize interpolation derivatives
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::LinearExponentialCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tmp, KsCache, exp_cache = cache
    alg = unwrap_alg(integrator, nothing)
    exp_method = ExpMethodGeneric()
    A = convert(AbstractMatrix, f.f) # assume f to be an ODEFunction wrapped around a linear operator

    if alg.krylov == :off
        E = exponential!(dt * A, exp_method, exp_cache)
        mul!(tmp, E, u)
    elseif alg.krylov == :simple
        Ks, expv_cache = KsCache
        arnoldi!(Ks, A, u; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop)
        expv!(tmp, dt, Ks; cache = expv_cache)
    else
        expv_timestep!(tmp, dt, A, u; adaptive = true, caches = KsCache,
            m = min(alg.m, size(A, 1)), iop = alg.iop,
            opnorm = integrator.opts.internalopnorm,
            tol = integrator.opts.reltol)
    end

    # Update integrator state
    u .= tmp
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

cay!(tmp, A) = mul!(tmp, inv(I - 1 / 2 * A), (I + 1 / 2 * A))
cay!(tmp, A::AbstractSciMLOperator) = @error "cannot multiply two SciMLOperators with mul!"
cay(A) = inv(I - 1 / 2 * A) * (I + 1 / 2 * A)

function initialize!(integrator, cache::CayleyEulerConstantCache)
    # Pre-start fsal
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)

    # Initialize interpolation derivatives
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::CayleyEulerConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator

    if f isa SplitFunction
        A = f.f1.f
    else  # f isa ODEFunction
        A = f.f
    end

    L = update_coefficients(A, uprev, p, t)
    L = convert(AbstractMatrix, L)

    V = cay(L * dt)
    u = (V * uprev) * transpose(V)

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::CayleyEulerCache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::CayleyEulerCache, repeat_step = false)
    @unpack t, dt, uprev, u, p, f = integrator
    @unpack k, V, tmp = cache
    mass_matrix = integrator.f.mass_matrix

    if f isa SplitFunction
        L = f.f1.f
    else  # f isa ODEFunction
        L = f.f
    end

    update_coefficients!(L, uprev, p, t)
    L = convert(AbstractMatrix, L)

    cay!(V, L * dt)
    mul!(tmp, uprev, transpose(V))
    mul!(u, V, tmp)

    # Update integrator state
    integrator.f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end
