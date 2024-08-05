using OrdinaryDiffEq
prob = ODEProblem((u, p, t) -> -u, BigFloat(1.0), (0.0, 1.0))
solve(prob, Tsit5())
solve(prob, KenCarp4())

# Function initial condition
f(u, p, t) = u
u0(p, t0) = ones(2)
prob = ODEProblem(f, u0, (0.0, 1.0))
sol = solve(prob, Tsit5())

# Test array partition outside of symplectic

u0 = fill(0.0, 2)
v0 = ones(2)

function f_ap(du, u, p, t)
    du.x[1] .= -2u.x[2]
    du.x[2] .= u.x[1]
end

u = ArrayPartition((u0, v0))

prob = ODEProblem(f_ap, u, (0.0, 5.0))
sol = solve(prob, Euler(), dt = 1 / 100)
