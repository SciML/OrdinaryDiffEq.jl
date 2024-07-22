using OrdinaryDiffEq
using LinearAlgebra
using NLsolve
using Test

p_inv = [500.0
         0.084
         4.69
         2.0
         400.0
         20.0
         0.2
         1000.0
         0.59
         736.0
         0.0
         0.0
         0.2
         1.27
         14.3
         0.0
         50.0
         0.2
         0.08
         0.003
         0.074
         0.2
         0.01 # Matt
         1.0  # V_ref
         1.0  # ω ref
         0.7  # Pref
         0    # Q_ref
         0.5  #Xtrans
         0.0  # Rtrans
         1.01 #Vm
         0.0] # Vθ

function vsm(dx, x, p, t)
    #PARAMETERS
    ω_lp = p[1]
    kp_pll = p[2]
    ki_pll = p[3]
    Ta = p[4]
    kd = p[5]
    kω = p[6]
    kq = p[7]
    ωf = p[8]
    kpv = p[9]
    kiv = p[10]
    kffv = p[11]
    rv = p[12]
    lv = p[13]
    kpc = p[14]
    kic = p[15]
    kffi = p[16]
    ωad = p[17]
    kad = p[18]
    lf = p[19]
    rf = p[20]
    cf = p[21]
    lg = p[22]
    rg = p[23]
    Vref = p[24]
    ωref = p[25]
    Pref = p[26]
    Qref = p[27]
    Xtrans = p[28]
    Rtrans = p[29]
    Vm = p[30]
    Vθ = p[31]

    #STATE INDEX AND STATES
    i__vi_filter, vi_filter = 1, x[1]
    i__γd_ic, γd_ic = 2, x[2]
    i__vq_pll, vq_pll = 3, x[3]
    i__γq_ic, γq_ic = 4, x[4]
    i__ir_filter, ir_filter = 5, x[5]
    i__ξd_ic, ξd_ic = 6, x[6]
    i__ϕd_ic, ϕd_ic = 7, x[7]
    i__ε_pll, ε_pll = 8, x[8]
    i__ir_cnv, ir_cnv = 9, x[9]
    i__vr_filter, vr_filter = 10, x[10]
    i__ω_oc, ω_oc = 11, x[11]
    i__ξq_ic, ξq_ic = 12, x[12]
    i__vd_pll, vd_pll = 13, x[13]
    i__q_oc, q_oc = 14, x[14]
    i__ϕq_ic, ϕq_ic = 15, x[15]
    i__θ_pll, θ_pll = 16, x[16]
    i__θ_oc, θ_oc = 17, x[17]
    i__ii_cnv, ii_cnv = 18, x[18]
    i__ii_filter, ii_filter = 19, x[19]
    i__ir_source, ir_source = 20, x[20]
    i__ii_source, ii_source = 21, x[21]

    ω_base = 60.0 * 2 * pi
    ω_sys = 1.0

    #PLL
    δω_pll = kp_pll * atan(vq_pll / vd_pll) + ki_pll * ε_pll
    ω_pll = δω_pll + ω_sys
    vd_filt_pll = sin(θ_pll + pi / 2) * vr_filter - cos(θ_pll + pi / 2) * vi_filter
    vq_filt_pll = cos(θ_pll + pi / 2) * vr_filter + sin(θ_pll + pi / 2) * vi_filter

    dx[i__vd_pll] = ω_lp * (vd_filt_pll - vd_pll)                                 #docs:(1a)
    dx[i__vq_pll] = ω_lp * (vq_filt_pll - vq_pll)                                 #docs:(1b)
    dx[i__ε_pll] = atan(vq_pll / vd_pll)                                         #docs:(1c)
    dx[i__θ_pll] = ω_base * δω_pll                                             #docs:(1d)

    #OUTER LOOP CONTROL
    pe = vr_filter * ir_filter + vi_filter * ii_filter                              #docs:(2d)
    qe = vi_filter * ir_filter - vr_filter * ii_filter                              #docs:(2e)
    v_ref_olc = Vref + kq * (Qref - q_oc)

    dx[i__ω_oc] = (Pref - pe - kd * (ω_oc - ω_pll) - kω * (ω_oc - ωref)) / Ta        #docs:(2a)
    dx[i__θ_oc] = ω_base * (ω_oc - ω_sys)                                            #docs:(2b)
    dx[i__q_oc] = ωf * (qe - q_oc)                                               #docs:(2c)

    #INNER LOOP CONTROL
    #reference transformations
    vd_filt_olc = sin(θ_oc + pi / 2) * vr_filter - cos(θ_oc + pi / 2) * vi_filter
    vq_filt_olc = cos(θ_oc + pi / 2) * vr_filter + sin(θ_oc + pi / 2) * vi_filter
    id_filt_olc = sin(θ_oc + pi / 2) * ir_filter - cos(θ_oc + pi / 2) * ii_filter
    iq_filt_olc = cos(θ_oc + pi / 2) * ir_filter + sin(θ_oc + pi / 2) * ii_filter
    id_cnv_olc = sin(θ_oc + pi / 2) * ir_cnv - cos(θ_oc + pi / 2) * ii_cnv
    iq_cnv_olc = cos(θ_oc + pi / 2) * ir_cnv + sin(θ_oc + pi / 2) * ii_cnv

    #Voltage control equations
    Vd_filter_ref = v_ref_olc - rv * id_filt_olc + ω_oc * lv * iq_filt_olc      #docs:(3g)
    Vq_filter_ref = -rv * iq_filt_olc - ω_oc * lv * id_filt_olc                 #docs:(3h)
    dx[i__ξd_ic] = Vd_filter_ref - vd_filt_olc                                 #docs:(3a)
    dx[i__ξq_ic] = Vq_filter_ref - vq_filt_olc                                 #docs:(3b)

    #current control equations
    Id_cnv_ref = (kpv * (Vd_filter_ref - vd_filt_olc) + kiv * ξd_ic -                           #docs:(3i)
                  cf * ω_oc * vq_filt_olc + kffi * id_filt_olc)
    Iq_cnv_ref = (kpv * (Vq_filter_ref - vq_filt_olc) +
                  kiv * ξq_ic +                           #docs:(3j)
                  cf * ω_oc * vd_filt_olc +
                  kffi * iq_filt_olc)
    dx[i__γd_ic] = Id_cnv_ref - id_cnv_olc                                      #docs:(3c)
    dx[i__γq_ic] = Iq_cnv_ref - iq_cnv_olc                                      #docs:(3d)

    #active damping equations
    Vd_cnv_ref = (kpc * (Id_cnv_ref - id_cnv_olc) + kic * γd_ic - lf * ω_oc * iq_cnv_olc +          #docs:(3k)
                  kffv * vd_filt_olc - kad * (vd_filt_olc - ϕd_ic))
    Vq_cnv_ref = (kpc * (Iq_cnv_ref - iq_cnv_olc) +
                  kic * γq_ic +
                  lf * ω_oc * id_cnv_olc +          #docs:(3l)
                  kffv * vq_filt_olc - kad * (vq_filt_olc - ϕq_ic))
    dx[i__ϕd_ic] = ωad * (vd_filt_olc - ϕd_ic)                                 #docs:(3e)
    dx[i__ϕq_ic] = ωad * (vq_filt_olc - ϕq_ic)                                  #docs:(3f)

    #LCL FILTER
    #reference transformations
    Vr_cnv = sin(θ_oc + pi / 2) * Vd_cnv_ref + cos(θ_oc + pi / 2) * Vq_cnv_ref
    Vi_cnv = -cos(θ_oc + pi / 2) * Vd_cnv_ref + sin(θ_oc + pi / 2) * Vq_cnv_ref

    Vr_pcc = Vm * cos(Vθ)
    Vi_pcc = Vm * sin(Vθ)

    dx[i__ir_cnv] = (ω_base / lf) *                                                #docs:(5a)
                    (Vr_cnv - vr_filter - rf * ir_cnv + ω_sys * lf * ii_cnv)
    dx[i__ii_cnv] = (ω_base / lf) *                                                #docs:(5b)
                    (Vi_cnv - vi_filter - rf * ii_cnv - ω_sys * lf * ir_cnv)
    dx[i__vr_filter] = (ω_base / cf) *                                             #docs:(5c)
                       (ir_cnv - ir_filter + ω_sys * cf * vi_filter)
    dx[i__vi_filter] = (ω_base / cf) *                                             #docs:(5d)
                       (ii_cnv - ii_filter - ω_sys * cf * vr_filter)
    dx[i__ir_filter] = (ω_base / lg) *                                             #docs:(5e)
                       (vr_filter - Vr_pcc - rg * ir_filter + ω_sys * lg * ii_filter)
    dx[i__ii_filter] = +(ω_base / lg) *                                            #docs:(5f)
                       (vi_filter - Vi_pcc - rg * ii_filter - ω_sys * lg * ir_filter)

    # Network interface algebraic equations 0 = I - YV
    line_currents = ((Vr_cnv + Vi_cnv * 1im) - (Vr_pcc + Vi_pcc * 1im)) /
                    (Rtrans + Xtrans * 1im)
    dx[i__ii_source] = ii_source - imag(line_currents)
    dx[i__ir_source] = ir_source - real(line_currents)
    return
end

u0 = [0.0,
    0.01,
    0.01,
    0.01,
    0.01,
    0.01,
    0.01,
    0.01,
    1.0,
    1.0,
    0.0,
    0.01,
    0.01,
    0.01,
    0.01,
    0.01,
    0.01,
    0.01,
    0.0,
    1.0,
    0.0
]
init_f! = (dx, x) -> vsm(dx, x, p_inv, 0)
res = nlsolve(init_f!, u0)

M = diagm(ones(21))
#Last 2 equations are algebraic
M[20, 20] = M[21, 21] = 0.0

f = ODEFunction(vsm, mass_matrix = M)

condition(u, t, integrator) = t == 1.0
affect!(integrator) = integrator.p[28] += 0.2
cb = DiscreteCallback(condition, affect!)

prob = ODEProblem(f, deepcopy(res.zero), (0, 20.0), deepcopy(p_inv))
refsol = solve(prob, Rodas4(), saveat = 0.1, callback = cb, tstops = [1.0], reltol = 1e-12,
    abstol = 1e-17)

for solver in (Rodas4, Rodas4P, Rodas5, Rodas5P, FBDF, QNDF, Rosenbrock23)
    @show solver
    prob = ODEProblem(f, deepcopy(res.zero), (0, 20.0), deepcopy(p_inv))
    sol = solve(
        prob, solver(), saveat = 0.1, callback = cb, tstops = [1.0], reltol = 1e-14,
        abstol = 1e-14)
    @test sol.retcode == ReturnCode.Success
    @test sol.t[end] == 20.0
    @test maximum(sol - refsol) < 1e-11
end

function hardstop!(du, u, p, t)
    pm, pg = p
    y, f_wall, dy = u
    du[1] = dy
    du[2] = ifelse(y <= 0, y, f_wall)
    du[3] = (-ifelse(t < 2, -pg * pm, pg * pm) - f_wall) / (-pm)
end

hardstop!(u, p, t) = (du = similar(u); hardstop!(du, u, p, t); du)

fun = ODEFunction(hardstop!, mass_matrix = Diagonal([1, 0, 1]))
prob1 = ODEProblem(fun, [5, 0, 0.0], (0, 4.0), [100, 10.0])
prob2 = ODEProblem(fun, [5, 0, 0.0], (0, 4.0), [100, 10.0])
for prob in [prob1, prob2]
    @test solve(prob, ImplicitEuler(), dt = 1 / 2^10, adaptive = false).retcode ==
          ReturnCode.ConvergenceFailure
end

condition2 = (u, t, integrator) -> t == 2
affect2! = integrator -> integrator.u[1] = 1e-6
cb = DiscreteCallback(condition2, affect2!)

@isdefined(N_FAILS) || const N_FAILS = Ref(0)
function choice_function(integrator)
    integrator.do_error_check = false
    if integrator.force_stepfail
        N_FAILS[] += 1
    else
        N_FAILS[] = 0
    end

    (N_FAILS[] > 3) + 1
end
simple_implicit_euler = ImplicitEuler(nlsolve = NLNewton(check_div = false,
    always_new = true))
alg_switch = CompositeAlgorithm((ImplicitEuler(), simple_implicit_euler), choice_function)

for prob in [prob1, prob2], alg in [simple_implicit_euler, alg_switch]
    sol = solve(prob, alg, callback = cb, dt = 1 / 2^10, adaptive = false)
    @test sol.retcode == ReturnCode.Success
    @test sol(0, idxs = 1) == 5
    @test sol(2, idxs = 1) == 0
    @test sol(4, idxs = 1) > 10
end