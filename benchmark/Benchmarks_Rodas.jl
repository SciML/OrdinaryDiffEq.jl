using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra
gr()
default(titlefont = 8, legendfontsize = 4, guidefont = 5, tickfont = 5, markersize = 2.5)

methods = [Rodas4(), Rodas4P(), Rodas4P2(), Rodas5(), Rodas5P()];
s_methods = ["Rodas4", "Rodas4P", "Rodas4P2", "Rodas5", "Rodas5P"]
n_method = length(methods)

#---- Problem definitions --------------------------------------------------------------------
title_1 = "Polynomial of order 6"
y0_1 = [0.0; 0.0];
tspan_1 = [0.0, 1.0];
M_1 = zeros(2, 2);
M_1[1, 1] = 1.0;
p_1 = [];
n_1 = [1, 2, 4, 8, 16, 32]; #- number of time steps
function ode_1(dy, y, p, t)
    dy[1] = 6 * t^5
    dy[2] = y[1] - y[2]
    nothing
end
function sol_1(t, p)
    [t .^ 6, t .^ 6]
end

title_2 = "Prothero Robinson"
M_2 = zeros(1, 1);
M_2[1, 1] = 1.0;
y0_2 = [0.0];
tspan_2 = [0.0, 2.0];
p_2 = [];
n_2 = [2, 4, 8, 16, 32, 64, 128, 256];
function ode_2(dy, y, p, t)
    lam = 1.0e5
    g = 10.0 - (10 + t) * exp(-t)
    dg = exp(-t) * (9 + t)
    dy[1] = -lam * (y[1] - g) + dg
    nothing
end
function sol_2(t, p)
    [10.0 - (10 + t) * exp(-t)]
end

title_3 = "Index 2"
M_3 = zeros(2, 2);
M_3[1, 1] = 1.0;
y0_3 = [-1.0; 1.0];
tspan_3 = [1.0, 2.0];
p_3 = [];
n_3 = [32, 64, 128, 256, 512, 1024];
function ode_3(du, u, p, t)
    du[1] = u[2]
    du[2] = u[1]^2 - 1 / t^2
end
function sol_3(t, p)
    [-1 / t, 1 / t .^ 2]
end

title_4 = "Parabol"
nx = 1000;
dx = 1.0 / nx;
n = nx - 1;
x = zeros(n);
M_4 = Matrix(1.0I, n, n);
for i in 1:n
    x[i] = i * dx
end
b1 = zeros(n);
b1[1] = 3 / (4 * dx^2);
b1[n] = 3 / (4 * dx^2);
J = zeros(n, n)
for i in 1:n
    J[i, i] = -2 / (dx^2)
    if i > 1
        J[i, i - 1] = 1 / (dx^2)
    end
    if i < n
        J[i, i + 1] = 1 / (dx^2)
    end
end
p_4 = x, b1, J
function ode_4(du, u, p, t)
    x, b1, J = p
    g = -(x .+ 1 / 2) .* (3 / 2 .- x) ./ (t + 1) .^ 2 .+ 2 / (1 + t)
    rb = b1 / (1 + t)
    du[:] = J * u .+ rb .+ g
end
function sol_4(t, p)
    x, b1, J = p
    u = (x .+ 1 / 2) .* (3 / 2 .- x) / (1 + t)
end
y0_4 = sol_4(0.0, p_4);
tspan_4 = [0.0, 0.1];
n_4 = [1, 2, 4, 8];

title_5 = "Hyperbol"
nx = 1000;
dx = 1.0 / nx;
n = nx;
x = zeros(n);
M_5 = Matrix(1.0I, n, n);
for i in 1:n
    x[i] = i * dx
end
b1 = zeros(n);
b1[1] = -1 / dx;
J = zeros(n, n);
for i in 1:n
    J[i, i] = 1 / dx
    if i > 1
        J[i, i - 1] = -1 / dx
    end
end
p_5 = x, b1, J
function ode_5(du, u, p, t)
    x, b1, J = p
    g = (t .- x) ./ (t .+ 1) .^ 2
    rb = b1 / (1 + t)
    du[:] = -J * u .- rb .+ g
end
function sol_5(t, p)
    x, b1, J = p
    u = (1 .+ x) ./ (1 .+ t)
end
y0_5 = sol_5(0.0, p_5);
tspan_5 = [0.0, 1];
n_5 = [4, 8, 16, 32, 64, 128];

#--- problems with constant stepsize -----
title_s = [title_1, title_2, title_3, title_4, title_5]
ode_s = [ode_1, ode_2, ode_3, ode_4, ode_5]
sol_s = [sol_1, sol_2, sol_3, sol_4, sol_5]
y0_s = [y0_1, y0_2, y0_3, y0_4, y0_5]
tspan_s = [tspan_1, tspan_2, tspan_3, tspan_4, tspan_5]
M_s = [M_1, M_2, M_3, M_4, M_5]
param_s = [p_1, p_2, p_3, p_4, p_5]
n_s = [n_1, n_2, n_3, n_4, n_5]

n_probs = length(ode_s);
Pl = [];

#--- Benchmarks for Rodas-type methods - --
for ip in 1:n_probs
    P_1 = plot()
    P_2 = plot()
    f = ODEFunction(ode_s[ip], mass_matrix = M_s[ip])
    tspan = tspan_s[ip]
    prob = ODEProblem(f, y0_s[ip], tspan, param_s[ip])
    n = n_s[ip]
    max_order = length(n)
    h = (tspan[2] - tspan[1]) ./ n
    err_end = zeros(max_order)
    err_dense = zeros(max_order)
    t_dense = range(tspan[1], tspan[2], length = 100)
    for m in 1:n_method
        for i in 1:max_order
            global sol = solve(prob, methods[m], dt = h[i], adaptive = false)
            #       println(sol.destats);
            err_end[i] = max(1.0e-16,
                             maximum(abs.(sol(tspan[2]) .- sol_s[ip](tspan[2], param_s[ip]))))
            #       err_end[i] = 1.0e-16
            #       for j = 1:length(sol.t)
            #         err_end[i] = max(err_end[i], maximum(abs.(sol.u[j] .- sol_s[ip](sol.t[j],param_s[ip]))))
            #       end
            err_dense[i] = 1.0e-16
            for j in 1:length(t_dense)
                err_dense[i] = max(err_dense[i],
                                   maximum(abs.(sol(t_dense[j]) .-
                                                sol_s[ip](t_dense[j], param_s[ip]))))
            end
        end
        P_1 = plot!(P_1, h, err_end, xaxis = :log, yaxis = :log, xlabel = "stepsize",
                    ylabel = L"$L_{\infty}$ error at $t_{\textrm{end}}$",
                    label = s_methods[m], marker = ([:circle :d]), legend = :bottomright)
        P_2 = plot!(P_2, h, err_dense, xaxis = :log, yaxis = :log, xlabel = "stepsize",
                    ylabel = L"dense $L_{\infty}$ error", label = s_methods[m],
                    marker = ([:circle :d]), legend = :bottomright)
    end
    push!(Pl, plot(P_1, P_2, layout = (1, 2), title = title_s[ip]))
end
p = plot(Pl[1], Pl[2], layout = (2, 1))
savefig(p, "rodas5p_bench_1.pdf")
p = plot(Pl[3], Pl[4], layout = (2, 1))
savefig(p, "rodas5p_bench_2.pdf")
p = plot(Pl[5], layout = (1, 1))
savefig(p, "rodas5p_bench_3.pdf")
