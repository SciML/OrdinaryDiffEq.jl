using OrdinaryDiffEqTaylorSeries
using Test

# Define a simple linear DAE system
function linear_dae!(res, du, u, p, t)
    # x' = y
    # 0 = x - t
    res[1] = du[1] - u[2]     # Differential equation
    res[2] = u[1] - t         # Algebraic constraint
    return nothing
end

# Initial conditions (consistent with constraints)
u0 = [0.0, 1.0]  # [x(0), y(0)]
du0 = [1.0, 0.0] # [x'(0), y'(0)]
tspan = (0.0, 5.0)
p = nothing

# Create problem and solve
prob = DAEProblem(linear_dae!, du0, u0, tspan, p)
sol = solve(prob, DAETS(), dt=0.1)

# Check solution against analytical solution: x(t) = t, y(t) = 1
for i in 1:length(sol.t)
    t = sol.t[i]
    @test isapprox(sol[1,i], t, rtol=1e-2)  # x(t) = t
    @test isapprox(sol[2,i], 1.0, rtol=1e-2) # y(t) = 1
end

println("Test completed successfully!") 