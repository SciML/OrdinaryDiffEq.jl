using OrdinaryDiffEq, Test

function fcn(du, u, p, t)
    du[1] = u[1]^2 - 1 / t^4
end

function f_jac(J, y, P, t)
    #-- numerical Jac by FD
    y_thres = P
    del = sqrt(eps(1.0))
    n = length(y)
    f0 = similar(y)
    f1 = similar(y)
    fcn(f0, y, P, t)
    for i = 1:n
        del_1 = del * max(abs(y[i]), y_thres)
        y1 = copy(y)
        y1[i] = y1[i] + del_1
        fcn(f1, y1, P, t)
        J[:, i] = (f1 - f0) / del_1
    end
end

tspan = (1.0, 1.0e3);
M = zeros(1, 1);

println("--- AD ---")

f = ODEFunction(fcn, mass_matrix = M)
problem = ODEProblem(f, [-1.0], tspan);
problem2 = ODEProblem(fcn, [-1.0], tspan);

sol = solve(problem, Rodas4P2(), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 6100

sol = solve(problem2, TRBDF2(), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 4200

println("--- FD central ---")

sol = solve(problem, Rodas4P2(autodiff = false, diff_type = Val{:central}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 6100
@test sol.t[end] == 1000

sol = solve(problem2, TRBDF2(autodiff = false, diff_type = Val{:central}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 4200
@test sol.t[end] == 1000

println("--- FD forward ---")

sol = solve(problem, Rodas4P2(autodiff = false, diff_type = Val{:forward}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 6100
@test sol.t[end] == 1000

sol = solve(problem2, TRBDF2(autodiff = false, diff_type = Val{:forward}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 4200
@test sol.t[end] == 1000

println("--- FD forward, y_thres = 1 ---")

y_thres = 1.0;
f = ODEFunction(fcn, mass_matrix = M, jac = f_jac)
problem = ODEProblem(f, [-1.0], tspan, y_thres);
sol = solve(problem, Rodas4P2(autodiff = false, diff_type = Val{:forward}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 6100
@test sol.t[end] == 1000

println("--- FD forward, y_thres = 1.0e-5 ---")

y_thres = 1.0e-5;
f = ODEFunction(fcn, mass_matrix = M, jac = f_jac)
problem = ODEProblem(f, [-1.0], tspan, y_thres);
sol = solve(problem, Rodas4P2(autodiff = false, diff_type = Val{:forward}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 6100
@test sol.t[end] == 1000

println("--- FD forward, y_thres = sqrt(eps) ---")

y_thres = sqrt(eps(1.0));
f = ODEFunction(fcn, mass_matrix = M, jac = f_jac)
problem = ODEProblem(f, [-1.0], tspan, y_thres);
sol = solve(problem, Rodas4P2(autodiff = false, diff_type = Val{:forward}), maxiters = Int(1e7), reltol = 1.0e-12, abstol = 1.0e-12);
@test sol.destats.naccept < 6100
@test sol.t[end] == 1000
