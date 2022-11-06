using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, DiffEqDevTools,
      SparseArrays
gr()
default(titlefont = 8, legendfontsize = 4, guidefont = 5, tickfont = 5, markersize = 2.5)#,  markerstrokewidth = 0) #thickness_scaling=0.8)

setups = [
    Dict(:alg => Rodas4()),
    Dict(:alg => Rodas4P()),
    Dict(:alg => Rodas4P2()),
    Dict(:alg => Rodas5()),
    Dict(:alg => Rodas5P()),
]
setups_1 = [
    Dict(:alg => Rodas4(; autodiff = false)),
    Dict(:alg => Rodas4P(; autodiff = false)),
    Dict(:alg => Rodas4P2(; autodiff = false)),
    Dict(:alg => Rodas5(; autodiff = false)),
    Dict(:alg => Rodas5P(; autodiff = false)),
]

#---- Problem definitions --------------------------------------------------------------------

title_1 = "Parabol"
nx = 250;
dx = 1.0 / nx;
n = nx - 1;
x = zeros(n);
M_1 = sparse(1.0I, n, n);
J = spzeros(n, n);
#M_1 = Matrix(1.0I, n, n); J = zeros(n,n)
numrun_1 = 10;
for i in 1:n
    x[i] = i * dx
end
b1 = zeros(n);
b1[1] = 3 / (4 * dx^2);
b1[n] = 3 / (4 * dx^2);
for i in 1:n
    J[i, i] = -2 / (dx^2)
    if i > 1
        J[i, i - 1] = 1 / (dx^2)
    end
    if i < n
        J[i, i + 1] = 1 / (dx^2)
    end
end
p_1 = x, b1, J
function ode_1(du, u, p, t)
    x, b1, J = p
    du[:] = J * u .+ b1 / (1 + t) .- (x .+ 1 / 2) .* (3 / 2 .- x) ./ (t + 1) .^ 2 .+
            2 / (1 + t)
    nothing
end
function jac_1(J, u, p, t)
    #  x, b1, J = p  
    J[:] = p[3]
    nothing
end
function sol_1(t, p)
    x, b1, J = p
    u = (x .+ 1 / 2) .* (3 / 2 .- x) / (1 + t)
end
y0_1 = sol_1(0.0, p_1);
tspan_1 = [0.0, 0.1];

title_2 = "Hyperbol"
nx = 250;
dx = 1.0 / nx;
n = nx;
x = zeros(n);
M_2 = sparse(1.0I, n, n);
J = spzeros(n, n);
#M_2 = Matrix(1.0I, n, n); J = zeros(n,n)
numrun_2 = 10;
for i in 1:n
    x[i] = i * dx
end
b1 = zeros(n);
b1[1] = -1 / dx;
for i in 1:n
    J[i, i] = 1 / dx
    if i > 1
        J[i, i - 1] = -1 / dx
    end
end
p_2 = x, b1, J
function ode_2(du, u, p, t)
    x, b1, J = p
    du[:] = -J * u .- b1 / (1 + t) .+ (t .- x) ./ (t .+ 1) .^ 2
    nothing
end
function jac_2(J, u, p, t)
    J[:] = -p[3]
    nothing
end
function sol_2(t, p)
    x, b1, J = p
    u = (1 .+ x) ./ (1 .+ t)
end
y0_2 = sol_2(0.0, p_2);
tspan_2 = [0.0, 1];

title_3 = "Pendulum Index-1";
title_4 = "Pendulum Index-2";
M_3 = Matrix(1.0I, 5, 5);
M_3[5, 5] = 0.0;
numrun_3 = 20;
numrun_4 = 5;
p_3 = 1;
p_4 = 2;
y0_3 = [2.0, 0.0, 0.0, 0.0, 0.0];
tspan_3 = [0.0, 10.0];
function ode_3(du, u, p, t)
    du[1] = u[2]
    du[2] = u[1] * u[5]
    du[3] = u[4]
    du[4] = -9.81 + u[3] * u[5]
    if p == 1
        du[5] = u[2]^2 + u[4]^2 + u[5] * (u[1]^2 + u[3]^2) - 9.81 * u[3] #-- Index-1
    else
        du[5] = u[1] * u[2] + u[3] * u[4] #-- Index-2
    end
end
function jac_3(J, y, p, t)
    if p == 1
        J[:] = [0 1 0 0 0; y[5] 0 0 0 y[1]; 0 0 0 1 0; 0 0 y[5] 0 y[3];
                2*y[5]*y[1] 2*y[2] 2 * y[5] * y[3]-9.81 2*y[4] y[1]^2+y[3]^2] #-- Index-1
    else
        J[:] = [0 1 0 0 0; y[5] 0 0 0 y[1]; 0 0 0 1 0; 0 0 y[5] 0 y[3];
                y[2] y[1] y[4] y[3] 0] #-- Index-2
    end
    nothing
end

#-- transistor amplifier from IVP testset
title_5 = "Transistor amplifier";
numrun_5 = 10;
ub = 6.0;
uf = 0.026;
alpha = 0.99;
beta = 1.0e-6;
r0 = 1000;
r1 = 9000;
r2 = 9000;
r3 = 9000;
r4 = 9000;
r5 = 9000;
r6 = 9000;
r7 = 9000;
r8 = 9000;
r9 = 9000;
p_5 = ub, uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9
tspan_5 = [0, 0.2];
y0_5 = [
    0.0,
    ub / (r2 / r1 + 1.0),
    ub / (r2 / r1 + 1.0),
    ub,
    ub / (r6 / r5 + 1.0),
    ub / (r6 / r5 + 1.0),
    ub,
    0.0,
]
c1 = 1.0e-6;
c2 = 2.0e-6;
c3 = 3.0e-6;
c4 = 4.0e-6;
c5 = 5.0e-6;
M_5 = zeros(8, 8);
M_5[1, 1] = -c1;
M_5[1, 2] = c1;
M_5[2, 1] = c1;
M_5[2, 2] = -c1;
M_5[3, 3] = -c2;
M_5[4, 4] = -c3;
M_5[4, 5] = c3;
M_5[5, 4] = c3;
M_5[5, 5] = -c3;
M_5[6, 6] = -c4;
M_5[7, 7] = -c5;
M_5[7, 8] = c5;
M_5[8, 7] = c5;
M_5[8, 8] = -c5;

function ode_5(du, u, p, t)
    ub, uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9 = p
    uet = 0.1 * sin(200 * pi * t)
    arg = min((u[2] - u[3]) / uf, 100.0)
    fac1 = beta * (exp(arg) - 1.0)
    arg = min((u[5] - u[6]) / uf, 100.0)
    fac2 = beta * (exp(arg) - 1.0)
    #  fac1  = beta*(exp((u[2]-u[3])/uf)-1.0)
    #  fac2  = beta*(exp((u[5]-u[6])/uf)-1.0)
    du[1] = (u[1] - uet) / r0
    du[2] = u[2] / r1 + (u[2] - ub) / r2 + (1.0 - alpha) * fac1
    du[3] = u[3] / r3 - fac1
    du[4] = (u[4] - ub) / r4 + alpha * fac1
    du[5] = u[5] / r5 + (u[5] - ub) / r6 + (1.0 - alpha) * fac2
    du[6] = u[6] / r7 - fac2
    du[7] = (u[7] - ub) / r8 + alpha * fac2
    du[8] = u[8] / r9
    nothing
end
function jac_5(J, u, p, t)
    ub, uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9 = p
    arg = (u[2] - u[3]) / uf
    fac1 = beta * exp(arg) / uf
    arg = (u[5] - u[6]) / uf
    fac2 = beta * exp(arg) / uf
    J .= 0.0
    J[1, 1] = 1.0 / r0
    J[2, 2] = 1.0 / r1 + 1.0 / r2 + (1.0 - alpha) * fac1
    J[2, 3] = -(1.0 - alpha) * fac1
    J[3, 2] = -fac1
    J[3, 3] = 1.0 / r3 + fac1
    J[4, 2] = alpha * fac1
    J[4, 3] = -alpha * fac1
    J[4, 4] = 1.0 / r4
    J[5, 5] = 1.0 / r5 + 1.0 / r6 + (1.0 - alpha) * fac2
    J[5, 6] = -(1.0 - alpha) * fac2
    J[6, 5] = -fac2
    J[6, 6] = 1.0 / r7 + fac2
    J[7, 5] = alpha * fac2
    J[7, 6] = -alpha * fac2
    J[7, 7] = 1.0 / r8
    J[8, 8] = 1.0 / r9
    nothing
end

title_6 = "PV-Network";
numrun_6 = 20;
c1 = -3.1037;
c2 = 1.0015;
c3 = 0.0032;
c4 = 1.3984e-09;
c5 = 0.4303;
c6 = 1.5 * 0.9562;
R0 = 0.2;
R1 = 0.5;
C = 4000.0;
q_max = 36000.0;
p_6 = c1, c2, c3, c4, c5, c6, R0, R1, C, q_max
y0_6 = [0.0, 11.856598910310167, 0.0, 2.9409008015416687, -2.940900801550821, 0.0, 9000.0]
M_6 = zeros(7, 7);
M_6[6, 6] = 1.0;
M_6[7, 7] = 1.0;
tspan_6 = [0.0, 36000.0];

function ode_6(dy, y, p, t)
    c1, c2, c3, c4, c5, c6, R0, R1, C, q_max = p
    U0 = y[1]
    U1 = y[2]
    iV = y[3]
    iPV = y[4]
    iB = y[5]
    uB = y[6]
    qB = y[7]
    dy[1] = U0
    dy[2] = iB + iPV - iV
    dy[3] = verbraucher(t) - iV * (U1 - U0)
    dy[4] = c1 + c2 * iPV + c3 * (U1 - U0) + c4 * (exp(c5 * iPV + c6 * (U1 - U0)) - 1)
    dy[5] = U1 - U0 - (ruhespannung(qB / q_max) - uB - R0 * iB)
    dy[6] = iB / C - uB / (R1 * C)
    dy[7] = -iB
    nothing
end
function verbraucher(t)
    ts = range(3600.0, 36000.0, step = 3600.0)
    return 50.0 * einaus(t, ts, 60.0)
    #  return 50
end
function ruhespannung(soc)
    pp = [6.8072, -10.5555, 6.2199, 10.2668] #-- Polynomkoeffizienten
    return ((pp[1] * soc + pp[2]) * soc + pp[3]) * soc + pp[4] #-- Polynomauswertung 
end
function fstep(t, t0, dauer)
    #-- atanh(0.999) = 3.8002, atanh(0.99) = 2.647
    s = 2.647 / dauer
    return (tanh(s * (t - t0)) + 1.0) / 2.0
end
function einaus(t, ts, schaltzeit)
    y = 0.0
    for i in 1:2:length(ts)
        y = y + fstep(t, ts[i], schaltzeit)
    end
    for i in 2:2:length(ts)
        y = y - fstep(t, ts[i], schaltzeit)
    end
    return y
end

title_7 = "Water tube problem";
numrun_7 = 10;
nu = 1.31e-6;
g = 9.8;
rho = 1.0e3;
rcrit = 2.3e3;
len = 1.0e3;
k = 2.0e-4;
d = 1.0;
b = 2.0e2;
a = pi * d^2 / 4.0;
mu = nu * rho;
c = b / (rho * g);
v = rho * len / a;
nnodes = 13;
nedges = 18;
ANN = zeros(nnodes, nnodes)
ANN[1, 2] = 1;
ANN[2, 3] = 1;
ANN[2, 6] = 1;
ANN[3, 4] = 1;
ANN[3, 5] = 1;
ANN[4, 5] = 1;
ANN[5, 10] = 1;
ANN[6, 5] = 1;
ANN[7, 4] = 1;
ANN[7, 8] = 1;
ANN[8, 5] = 1;
ANN[8, 10] = 1;
ANN[9, 8] = 1;
ANN[11, 9] = 1;
ANN[11, 12] = 1;
ANN[12, 7] = 1;
ANN[12, 8] = 1;
ANN[13, 11] = 1;
ANE = zeros(nnodes, nedges);
iann = [1, 2, 2, 3, 3, 4, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 12, 13]
jann = [2, 3, 6, 4, 5, 5, 10, 5, 4, 8, 5, 10, 8, 9, 12, 7, 8, 11]
for i in 1:nnodes
    for j in 1:nedges
        if iann[j] == i
            ANE[i, j] = -1
        end
        if jann[j] == i
            ANE[i, j] = 1
        end
    end
end
dein = zeros(nnodes);
rghres = zeros(nedges);
fdba = zeros(nedges);
p_7 = nu, rho, rcrit, len, k, d, a, mu, c, v, nnodes, nedges, iann, jann, ANE, dein, rghres,
      fdba
tspan_7 = [0.0, 17 * 3600.0];
neq = 49;
y0_7 = zeros(neq);
y0_7[19:36] .= 0.47519404529185289807e-1;
y0_7[37:49] .= 109800.0;
M_7 = zeros(neq, neq);
for i in 1:18
    M_7[i, i] = 1
end
i5 = 2 * nedges + 5;
i8 = 2 * nedges + 8;
M_7[i5, i5] = 1;
M_7[i8, i8] = 1;

function ode_7(du, u, p, t)
    nu, rho, rcrit, len, k, d, a, mu, c, v, nnodes, nedges, iann, jann, ANE, dein, rghres, fdba = p
    phi = u[1:nedges]
    lambda = u[(nedges + 1):(2 * nedges)]
    pp = u[(2 * nedges + 1):end]
    that = t / 3600.0
    dein[1] = -sin(exp(-that) - 1.0) * exp(-that) / (3600.0 * 200)
    dein[13] = -sin(exp(-that) - 1.0) * exp(-that) / (3600.0 * 80)
    dein[10] = -(12 * that^3 - 3 * 92 * that^2 + 2 * 720 * that) / (1.0e6 * 3600)
    for n in 1:nedges
        pL = pp[iann[n]]
        pR = pp[jann[n]]
        rtla = sqrt(lambda[n])
        r = abs(phi[n] * d / (nu * a))
        rghres[n] = 1.0 / rtla - 1.74 +
                    2.0 * log10(2.0 * k / d + 18.7 / (max(r, rcrit) * rtla))
        fdba1 = pL - pR - lambda[n] * rho * len * phi[n]^2 / (a^2 * d)
        fdba0 = pL - pR - 32.0 * mu * len * phi[n] / (a * d^2)
        fac0 = 0.9
        fac1 = 1.1
        alpha = (r / rcrit - fac0) / (fac1 - fac0) #-- continous transition
        alpha = max(alpha, 0.0)
        alpha = min(alpha, 1.0)
        fdba[n] = ((1.0 - alpha) * fdba0 + alpha * fdba1) / v
        #    if r > rcrit
        #      fdba[n] = fdba1/v
        #    else
        #      fdba[n] = fdba0/v
        #    end
    end
    #-- 
    netflow = dein .+ ANE * fdba
    AP = (ANE * phi) / c  #-- Storage nodes
    netflow[5] = AP[5]
    netflow[8] = AP[8]
    #--
    du[:] = [fdba; rghres; netflow]
    nothing
end

title_8 = "Pollution";
numrun_8 = 20;
k1 = 0.35e0;
k2 = 0.266e2;
k3 = 0.123e5;
k4 = 0.86e-3;
k5 = 0.82e-3;
k6 = 0.15e5;
k7 = 0.13e-3;
k8 = 0.24e5;
k9 = 0.165e5;
k10 = 0.9e4;
k11 = 0.22e-1;
k12 = 0.12e5;
k13 = 0.188e1;
k14 = 0.163e5;
k15 = 0.48e7;
k16 = 0.35e-3;
k17 = 0.175e-1;
k18 = 0.1e9;
k19 = 0.444e12;
k20 = 0.124e4;
k21 = 0.21e1;
k22 = 0.578e1;
k23 = 0.474e-1;
k24 = 0.178e4;
k25 = 0.312e1;

p_8 = k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19,
      k20, k21, k22, k23, k24, k25

y0_8 = zeros(20)
y0_8[2] = 0.2;
y0_8[4] = 0.04;
y0_8[7] = 0.1;
y0_8[8] = 0.3;
y0_8[9] = 0.01;
y0_8[17] = 0.007;
tspan_8 = [0.0, 60.0];
M_8 = Matrix(1.0I, 20, 20);

function ode_8(dy, y, p, t)
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25 = p
    r1 = k1 * y[1]
    r2 = k2 * y[2] * y[4]
    r3 = k3 * y[5] * y[2]
    r4 = k4 * y[7]
    r5 = k5 * y[7]
    r6 = k6 * y[7] * y[6]
    r7 = k7 * y[9]
    r8 = k8 * y[9] * y[6]
    r9 = k9 * y[11] * y[2]
    r10 = k10 * y[11] * y[1]
    r11 = k11 * y[13]
    r12 = k12 * y[10] * y[2]
    r13 = k13 * y[14]
    r14 = k14 * y[1] * y[6]
    r15 = k15 * y[3]
    r16 = k16 * y[4]
    r17 = k17 * y[4]
    r18 = k18 * y[16]
    r19 = k19 * y[16]
    r20 = k20 * y[17] * y[6]
    r21 = k21 * y[19]
    r22 = k22 * y[19]
    r23 = k23 * y[1] * y[4]
    r24 = k24 * y[19] * y[1]
    r25 = k25 * y[20]

    dy[1] = -r1 - r10 - r14 - r23 - r24 + r2 + r3 + r9 + r11 + r12 + r22 + r25
    dy[2] = -r2 - r3 - r9 - r12 + r1 + r21
    dy[3] = -r15 + r1 + r17 + r19 + r22
    dy[4] = -r2 - r16 - r17 - r23 + r15
    dy[5] = -r3 + r4 + r4 + r6 + r7 + r13 + r20
    dy[6] = -r6 - r8 - r14 - r20 + r3 + r18 + r18
    dy[7] = -r4 - r5 - r6 + r13
    dy[8] = r4 + r5 + r6 + r7
    dy[9] = -r7 - r8
    dy[10] = -r12 + r7 + r9
    dy[11] = -r9 - r10 + r8 + r11
    dy[12] = r9
    dy[13] = -r11 + r10
    dy[14] = -r13 + r12
    dy[15] = r14
    dy[16] = -r18 - r19 + r16
    dy[17] = -r20
    dy[18] = r20
    dy[19] = -r21 - r22 - r24 + r23 + r25
    dy[20] = -r25 + r24
end

function jac_8(J, y, p, t)
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25 = p
    J .= 0.0
    J[1, 1] = -k1 - k10 * y[11] - k14 * y[6] - k23 * y[4] - k24 * y[19]
    J[1, 11] = -k10 * y[1] + k9 * y[2]
    J[1, 6] = -k14 * y[1]
    J[1, 4] = -k23 * y[1] + k2 * y[2]
    J[1, 19] = -k24 * y[1] + k22
    J[1, 2] = k2 * y[4] + k9 * y[11] + k3 * y[5] + k12 * y[10]
    J[1, 13] = k11
    J[1, 20] = k25
    J[1, 5] = k3 * y[2]
    J[1, 10] = k12 * y[2]

    J[2, 4] = -k2 * y[2]
    J[2, 5] = -k3 * y[2]
    J[2, 11] = -k9 * y[2]
    J[2, 10] = -k12 * y[2]
    J[2, 19] = k21
    J[2, 1] = k1
    J[2, 2] = -k2 * y[4] - k3 * y[5] - k9 * y[11] - k12 * y[10]

    J[3, 1] = k1
    J[3, 4] = k17
    J[3, 16] = k19
    J[3, 19] = k22
    J[3, 3] = -k15

    J[4, 4] = -k2 * y[2] - k16 - k17 - k23 * y[1]
    J[4, 2] = -k2 * y[4]
    J[4, 1] = -k23 * y[4]
    J[4, 3] = k15

    J[5, 5] = -k3 * y[2]
    J[5, 2] = -k3 * y[5]
    J[5, 7] = 2k4 + k6 * y[6]
    J[5, 6] = k6 * y[7] + k20 * y[17]
    J[5, 9] = k7
    J[5, 14] = k13
    J[5, 17] = k20 * y[6]

    J[6, 6] = -k6 * y[7] - k8 * y[9] - k14 * y[1] - k20 * y[17]
    J[6, 7] = -k6 * y[6]
    J[6, 9] = -k8 * y[6]
    J[6, 1] = -k14 * y[6]
    J[6, 17] = -k20 * y[6]
    J[6, 2] = k3 * y[5]
    J[6, 5] = k3 * y[2]
    J[6, 16] = 2k18

    J[7, 7] = -k4 - k5 - k6 * y[6]
    J[7, 6] = -k6 * y[7]
    J[7, 14] = k13

    J[8, 7] = k4 + k5 + k6 * y[6]
    J[8, 6] = k6 * y[7]
    J[8, 9] = k7

    J[9, 9] = -k7 - k8 * y[6]
    J[9, 6] = -k8 * y[9]

    J[10, 10] = -k12 * y[2]
    J[10, 2] = -k12 * y[10] + k9 * y[11]
    J[10, 9] = k7
    J[10, 11] = k9 * y[2]

    J[11, 11] = -k9 * y[2] - k10 * y[1]
    J[11, 2] = -k9 * y[11]
    J[11, 1] = -k10 * y[11]
    J[11, 9] = k8 * y[6]
    J[11, 6] = k8 * y[9]
    J[11, 13] = k11

    J[12, 11] = k9 * y[2]
    J[12, 2] = k9 * y[11]

    J[13, 13] = -k11
    J[13, 11] = k10 * y[1]
    J[13, 1] = k10 * y[11]

    J[14, 14] = -k13
    J[14, 10] = k12 * y[2]
    J[14, 2] = k12 * y[10]

    J[15, 1] = k14 * y[6]
    J[15, 6] = k14 * y[1]

    J[16, 16] = -k18 - k19
    J[16, 4] = k16

    J[17, 17] = -k20 * y[6]
    J[17, 6] = -k20 * y[17]

    J[18, 17] = k20 * y[6]
    J[18, 6] = k20 * y[17]

    J[19, 19] = -k21 - k22 - k24 * y[1]
    J[19, 1] = -k24 * y[19] + k23 * y[4]
    J[19, 4] = k23 * y[1]
    J[19, 20] = k25

    J[20, 20] = -k25
    J[20, 1] = k24 * y[19]
    J[20, 19] = k24 * y[1]

    return
end

#--- collect problems -----
title_s = [title_1, title_2, title_3, title_4, title_5, title_6, title_7, title_8]
ode_s = [ode_1, ode_2, ode_3, ode_3, ode_5, ode_6, ode_7, ode_8]
y0_s = [y0_1, y0_2, y0_3, y0_3, y0_5, y0_6, y0_7, y0_8]
tspan_s = [tspan_1, tspan_2, tspan_3, tspan_3, tspan_5, tspan_6, tspan_7, tspan_8]
M_s = [M_1, M_2, M_3, M_3, M_5, M_6, M_7, M_8]
Jac_s = [jac_1, jac_2, jac_3, jac_3, jac_5, [], [], jac_8]
param_s = [p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8]
numrun_s = [numrun_1, numrun_2, numrun_3, numrun_4, numrun_5, numrun_6, numrun_7, numrun_8]

n_probs = length(ode_s);
Pl = [];

#--- Benchmarks for Rodas-type methods - --
for ip in 1:n_probs
    println(title_s[ip])
    P_1 = plot()
    P_2 = plot()
    if Jac_s[ip] == []
        f = ODEFunction(ode_s[ip], mass_matrix = M_s[ip])
    else
        f = ODEFunction(ode_s[ip], mass_matrix = M_s[ip], jac = Jac_s[ip])
    end
    tspan = tspan_s[ip]
    prob = ODEProblem(f, y0_s[ip], tspan, param_s[ip])
    if ip == 4 #-- take test_sol from ip=3
        abstols = 1.0 ./ 10.0 .^ (3:6)
        reltols = 1.0 ./ 10.0 .^ (3:6)
    else
        if ip == 7
            sol = solve(prob, Rodas5P(; autodiff = false), abstol = 1.0e-14,
                        reltol = 1.0e-14, maxiters = Int(1e8))
        else
            sol = solve(prob, Rodas5P(), abstol = 4.0e-14, reltol = 4.0e-14,
                        maxiters = Int(1e8))
        end
        display(plot(sol))
        global test_sol = TestSolution(sol)
        abstols = 1.0 ./ 10.0 .^ (3:9)
        reltols = 1.0 ./ 10.0 .^ (3:9)
        if (ip == 6)
            abstols = 1.0 ./ 10.0 .^ (7:12)
            reltols = 1.0 ./ 10.0 .^ (7:12)
        end
    end
    println("Test_sol ready")
    if ip == 7
        wp = WorkPrecisionSet(prob, abstols, reltols, setups_1; save_everystep = true,
                              error_estimate = :l2, appxsol = test_sol, maxiters = Int(1e8),
                              numruns = numrun_s[ip])
    else
        wp = WorkPrecisionSet(prob, abstols, reltols, setups; save_everystep = true,
                              error_estimate = :l2, appxsol = test_sol, maxiters = Int(1e8),
                              numruns = numrun_s[ip])
    end
    display(plot(wp, title = title_s[ip]))
    P_1 = plot(wp)
    if ip == 7
        wp = WorkPrecisionSet(prob, abstols, reltols, setups_1; save_everystep = true,
                              error_estimate = :L2, dense_errors = true, appxsol = test_sol,
                              maxiters = Int(1e8), numruns = numrun_s[ip])
    else
        wp = WorkPrecisionSet(prob, abstols, reltols, setups; save_everystep = true,
                              error_estimate = :L2, dense_errors = true, appxsol = test_sol,
                              maxiters = Int(1e8), numruns = numrun_s[ip])
    end
    P_2 = plot(wp)
    display(plot(wp, title = title_s[ip]))
    push!(Pl, plot(P_1, P_2, layout = (1, 2), title = title_s[ip]))
end
#plot(Pl[1])
p = plot(Pl[1], Pl[2], layout = (2, 1))
savefig(p, "rodas5p_bench_4.pdf")
p = plot(Pl[3], Pl[4], layout = (2, 1))
savefig(p, "rodas5p_bench_5.pdf")
p = plot(Pl[5], Pl[6], layout = (2, 1))
savefig(p, "rodas5p_bench_6.pdf")
p = plot(Pl[7], Pl[8], layout = (2, 1))
savefig(p, "rodas5p_bench_7.pdf")
