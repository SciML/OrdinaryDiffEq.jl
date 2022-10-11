using OrdinaryDiffEq, SparseArrays, LinearAlgebra, Test

function explicit_fun(du, _, _, _)
    du .= 0
    nothing
end

#Capillary leveling in an axisymmetric polar domain
@views function capillary_leveling(du, u, dr, nodes, boundaries)
    #Boundary conditions
    h = similar(u, length(u) + 4)
    h[3:(end - 2)] .= u
    h[2:-1:1] .= u[1:2]
    h[(end - 1):end] .= u[end]
    tmp = similar(h)

    #r*dh/dr
    @. tmp[1:(end - 1)] = (h[2:end] - h[1:(end - 1)]) / dr * boundaries
    #Capillary pressure p=-1/r*d/dr(r*dh/dr)
    @. tmp[2:(end - 1)] = -(tmp[2:(end - 1)] - tmp[1:(end - 2)]) / dr / nodes[2:(end - 1)]
    #Disjoining pressure A/h^3
    @. tmp[2:(end - 1)] -= 1E-3 / h[2:(end - 1)]^3

    #r*h^3*dp/dr
    @. tmp[2:(end - 2)] = (tmp[3:(end - 1)] - tmp[2:(end - 2)]) / dr *
                          boundaries[2:(end - 1)] * (h[2:(end - 2)]^3 + h[3:(end - 1)]^3) /
                          2
    #1/3*1/r*d/dr(r*h^3*dp/dr)
    @. du = (tmp[3:(end - 2)] - tmp[2:(end - 3)]) / dr / nodes[3:(end - 2)] / 3
    nothing
end

#Initial droplet shape
function droplet(r, h_p)
    if r <= 1
        h_p + (2 - h_p) * (1 - r^2)
    else
        h_p
    end
end

tspan = [0, 1E-2]
N = 1000
dr = 4 / N
r = @. dr / 2 + dr * (0:(N - 1))
h0 = map(x -> droplet(x, 1E-4), r)
r = [r[2]; r[1]; r; r[end]; r[end]]
boundaries = @. (r[1:(end - 1)] + r[2:end]) / 2

jac_proto = spdiagm(0 => ones(N), -1 => ones(N - 1), 1 => ones(N - 1), -2 => ones(N - 2),
                    2 => ones(N - 2))

fun = ODEFunction((dh, h, p, t) -> capillary_leveling(dh, h, dr, r, boundaries),
                  jac_prototype = jac_proto)
prob = ODEProblem(fun, h0, tspan)

#Should be functionally equivalent to above, explicit_fun does nothing
sfun = SplitFunction(fun, explicit_fun, jac_prototype = jac_proto)
sprob = ODEProblem(sfun, h0, tspan)

#CFNLIRK3 has same erroneous FSAL logic as KenCarp solvers
#Can't efficiently test with stiff problem (to cause dtmin issue) because it requires constant dt
@test_broken solve(sprob, CFNLIRK3(), reltol = 1E-8, dt = 1E-5).retcode == :Success

for Alg in (KenCarp3, KenCarp4, KenCarp5, KenCarp47, KenCarp58)
    print(Alg)
    sol = solve(prob, Alg(), reltol = 1E-8)
    @test sol.retcode == :Success

    split_sol = solve(sprob, Alg(), reltol = 1E-8)
    @test split_sol.retcode == :Success

    L2 = norm(sol[end] .- split_sol[end]) / sqrt(N)
    println(" ", L2)
    @test L2 < 1E-6
end
