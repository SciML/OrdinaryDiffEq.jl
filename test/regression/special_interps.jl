using OrdinaryDiffEq, Test, LinearAlgebra

function f(du, u, p, t)
    du[1] = dx = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = dy = -p[3] * u[2] + p[4] * u[1] * u[2]
end
function foop(u, p, t)
    dx = p[1] * u[1] - p[2] * u[1] * u[2]
    dy = -p[3] * u[2] + p[4] * u[1] * u[2]
    [dx, dy]
end

p = [1.5, 1.0, 3.0, 1.0]
prob = ODEProblem(f, [1.0; 1.0], (0.0, 10.0), p)
proboop = ODEProblem(foop, [1.0; 1.0], (0.0, 10.0), p)

SPECIAL_INTERPS = [Tsit5(), DP5(), SSPRK22(), OwrenZen3(), OwrenZen4(), OwrenZen5(),
    BS5(), Vern6(), Vern7(), Vern8(), Vern9(), DP8(), Rosenbrock23(),
    Rodas4(), Rodas5()]

y1 = zeros(2);
y2 = zeros(2);
for alg in SPECIAL_INTERPS
    @show alg
    sol = solve(prob, alg, dt = 0.0033, abstol = 1e-14, reltol = 1e-14)
    soloop = solve(proboop, alg, adaptive = false, tstops = sol.t, abstol = 1e-14,
                   reltol = 1e-14)
    @test maximum(norm(soloop(t) - sol(t)) for t in 0:0.001:10) < 1e-10
    @test maximum(norm(soloop(y1, t) - sol(y2, t)) for t in 0:0.001:10) < 1e-10
end
