using OrdinaryDiffEq

function SinCosDiffEqToSolve!(ddu, du, u, p, t)
    ddu[1] = -u[1]
end

prob = SecondOrderODEProblem(SinCosDiffEqToSolve!, [2.0], [3.0], (0.0, 1.0))
sol = solve(prob, Tsit5())
sol = solve(prob, Vern9())
sol = solve(prob, ROCK4())
sol = solve(prob, Rodas4())
sol = solve(prob, AutoVern9(Rodas4()))
sol = solve(prob, TRBDF2())
