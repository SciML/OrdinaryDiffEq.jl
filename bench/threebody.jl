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
tspan = (0.0, 100.0)
p = SA[μ]
prob = ODEProblem{true, SciMLBase.FullSpecialize}(pcr3bp!, q0, tspan, p)

ref = solve(prob, TaylorMethod(32), abstol = 1e-20)

@benchmark solve($prob, $(TaylorMethod(25)), abstol = 1e-15)
@benchmark solve($prob, $(ExplicitTaylor(order = Val(11))), abstol = 1e-10, reltol = 1e-10)
@benchmark solve($prob, $(Vern9()), abstol = 1e-10, reltol = 1e-10)

setups = [Dict(:alg => Vern7())
          Dict(:alg => Vern9())
          Dict(:alg => TaylorMethod(12))
          Dict(:alg => TaylorMethod(16))
          Dict(:alg => TaylorMethod(20))
          Dict(:alg => TaylorMethod(24))]
abstols = 10.0 .^ (-19:-15)
reltols = 10.0 .^ (-15:-11)
names = ["Vern7", "Vern9", "TI12", "TI16",
    "TI20", "TI24"]
wp = WorkPrecisionSet([prob], abstols, reltols, setups; names = names, appxsol = [ref],
    save_everystep = false, numruns = 100)
plot(wp)
