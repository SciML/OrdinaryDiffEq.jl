using OrdinaryDiffEq, OrdinaryDiffEqTaylorSeries, TaylorIntegration, StaticArrays,
      DiffEqDevTools, BenchmarkTools, Plots

@taylorize function pcr3bp!(dq, q, param, t)
    local μ = param[1]
    local onemμ = 1 - μ
    x1 = q[1] - μ
    x1sq = x1^2
    y = q[2]
    ysq = y^2
    r1_1p5 = (x1sq + ysq)^1.5
    x2 = q[1] + onemμ
    x2sq = x2^2
    r2_1p5 = (x2sq + ysq)^1.5
    dq[1] = q[3] + q[2]
    dq[2] = q[4] - q[1]
    dq[3] = (-((onemμ * x1) / r1_1p5) - ((μ * x2) / r2_1p5)) + q[4]
    dq[4] = (-((onemμ * y) / r1_1p5) - ((μ * y) / r2_1p5)) - q[3]
    return nothing
end

const μ = 0.01
V(x, y) = -(1 - μ) / sqrt((x - μ)^2 + y^2) - μ / sqrt((x + 1 - μ)^2 + y^2)
H(x, y, px, py) = (px^2 + py^2) / 2 - (x * py - y * px) + V(x, y)
H(x) = H(x...)
J0 = -1.58
function py!(q0, J0)
    @assert iszero(q0[2]) && iszero(q0[3]) # q0[2] and q0[3] have to be equal to zero
    q0[4] = q0[1] + sqrt(q0[1]^2 - 2(V(q0[1], q0[2]) - J0))
    nothing
end
q0 = [-0.8, 0.0, 0.0, 0.0]
py!(q0, J0)
tspan = (0.0, 2000.0)
p = SA[μ]
prob = ODEProblem{true, SciMLBase.FullSpecialize}(pcr3bp!, q0, tspan, p)

ref = solve(prob, Vern9(), abstol = 1e-10, reltol = 1e-10)

@benchmark solve($prob, $(TaylorMethod(25)), abstol = 1e-15)
@benchmark solve($prob, $(ExplicitTaylor(order = Val(25))), abstol = 1e-15, reltol = 1e-15)
@benchmark solve($prob, $(Vern9()), abstol = 1e-15, reltol = 1e-15)
