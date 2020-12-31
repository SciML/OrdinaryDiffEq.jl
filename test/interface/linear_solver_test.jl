using Test, OrdinaryDiffEq

prob = ODEProblem(function (du,u,p,t)
                    du[1] = 10.0(u[2]-u[1])
                    du[2] = u[1]*(28.0-u[3]) - u[2]
                    du[3] = u[1]*u[2] - (8/3)*u[3]
                  end,
                  [1.0;0.0;0.0], (0.0, 100.0))

function buildW(...)::NTuple{3}
    nothing #= left preconditioner =#, nothing #= right preconditioner =#
end

nlsolve = NLNewton(buildW, )
sol = solve(prob, TRBDF2(nlsolve=nlsolve))
