using OrdinaryDiffEq, Test
f(du, u, p, t) = du .= u
prob = ODEProblem(f, [1.0], (0.0, 1.0))

i = init(prob, Tsit5())
resize!(i, 5)
@test length(i.cache.u) == 5
@test length(i.cache.uprev) == 5
@test length(i.cache.k1) == 5
@test length(i.cache.k2) == 5
@test length(i.cache.k3) == 5
@test length(i.cache.k4) == 5
@test length(i.cache.k5) == 5
@test length(i.cache.k6) == 5
@test length(i.cache.k7) == 5
solve!(i)

i = init(prob, ImplicitEuler())
resize!(i, 5)
@test length(i.cache.atmp) == 5
@test length(i.cache.uprev) == 5
# nlsolver fields
@test length(i.cache.nlsolver.z) == 5
@test length(i.cache.nlsolver.ztmp) == 5
@test length(i.cache.nlsolver.tmp) == 5
# nlsolver cache fields
@test length(i.cache.nlsolver.cache.ustep) == 5
@test length(i.cache.nlsolver.cache.k) == 5
@test length(i.cache.nlsolver.cache.atmp) == 5
@test length(i.cache.nlsolver.cache.dz) == 5
@test size(i.cache.nlsolver.cache.J) == (5, 5)
@test size(i.cache.nlsolver.cache.W) == (5, 5)
@test length(i.cache.nlsolver.cache.du1) == 5
@test length(i.cache.nlsolver.cache.jac_config.fx) == 5
@test length(i.cache.nlsolver.cache.jac_config.dx) == 5
@test length(i.cache.nlsolver.cache.jac_config.t) == 5
@test length(i.cache.nlsolver.cache.jac_config.p) == 5
@test length(i.cache.nlsolver.cache.weight) == 5
solve!(i)

i = init(prob, ImplicitEuler(; autodiff = false))
resize!(i, 5)
@test length(i.cache.atmp) == 5
@test length(i.cache.uprev) == 5
# nlsolver fields
@test length(i.cache.nlsolver.z) == 5
@test length(i.cache.nlsolver.ztmp) == 5
@test length(i.cache.nlsolver.tmp) == 5
# nlsolver cache fields
@test length(i.cache.nlsolver.cache.ustep) == 5
@test length(i.cache.nlsolver.cache.k) == 5
@test length(i.cache.nlsolver.cache.atmp) == 5
@test length(i.cache.nlsolver.cache.dz) == 5
@test size(i.cache.nlsolver.cache.J) == (5, 5)
@test size(i.cache.nlsolver.cache.W) == (5, 5)
@test length(i.cache.nlsolver.cache.du1) == 5
@test length(i.cache.nlsolver.cache.jac_config.x1) == 5
@test length(i.cache.nlsolver.cache.jac_config.fx) == 5
@test length(i.cache.nlsolver.cache.jac_config.fx1) == 5
@test length(i.cache.nlsolver.cache.weight) == 5
solve!(i)

i = init(prob, Rosenbrock23())
resize!(i, 5)
@test length(i.cache.u) == 5
@test length(i.cache.uprev) == 5
@test length(i.cache.k₁) == 5
@test length(i.cache.k₂) == 5
@test length(i.cache.k₃) == 5
@test length(i.cache.du1) == 5
@test length(i.cache.du2) == 5
@test length(i.cache.f₁) == 5
@test length(i.cache.fsalfirst) == 5
@test length(i.cache.fsallast) == 5
@test length(i.cache.dT) == 5
@test length(i.cache.tmp) == 5
@test size(i.cache.J) == (5, 5)
@test size(i.cache.W) == (5, 5)
@test length(i.cache.linsolve_tmp) == 5
@test length(i.cache.jac_config.fx) == 5
@test length(i.cache.jac_config.dx) == 5
@test length(i.cache.jac_config.t) == 5
@test length(i.cache.jac_config.p) == 5
solve!(i)

i = init(prob, Rosenbrock23(; autodiff = false))
resize!(i, 5)
@test length(i.cache.u) == 5
@test length(i.cache.uprev) == 5
@test length(i.cache.k₁) == 5
@test length(i.cache.k₂) == 5
@test length(i.cache.k₃) == 5
@test length(i.cache.du1) == 5
@test length(i.cache.du2) == 5
@test length(i.cache.f₁) == 5
@test length(i.cache.fsalfirst) == 5
@test length(i.cache.fsallast) == 5
@test length(i.cache.dT) == 5
@test length(i.cache.tmp) == 5
@test size(i.cache.J) == (5, 5)
@test size(i.cache.W) == (5, 5)
@test length(i.cache.linsolve_tmp) == 5
@test length(i.cache.jac_config.x1) == 5
@test length(i.cache.jac_config.fx) == 5
@test length(i.cache.jac_config.fx1) == 5
solve!(i)

function f(du, u, p, t)
    du[1] = 2.0 * u[1] - 1.2 * u[1] * u[2]
    du[2] = -3 * u[2] + u[1] * u[2]
    for i in 3:length(u)
        du[i] = 0.0
    end
end
function f_jac(J, u, p, t)
    J[1, 1] = 2.0 - 1.2 * u[2]
    J[1, 2] = -1.2 * u[1]
    J[2, 1] = 1 * u[2]
    J[2, 2] = -3 + u[1]
    for i in 3:length(u)
        for j in 3:length(u)
            if i == j
                J[i, j] = 1.0
            else
                J[i, j] = 0.0
            end
        end
    end
    nothing
end
ff = ODEFunction(f; jac = f_jac, jac_prototype = [1.0 1.0; 1.0 1.0])

cb = DiscreteCallback((u, t, integ) -> true, integ -> @views(integ.u[3:5]) .= 0)
prob = ODEProblem(ff, [1.0, 1.0], (0.0, 1.0))
i = init(prob, ImplicitEuler(), callback = cb)
resize!(i, 5)
solve!(i)

function dsdt(ds, s, _, t)
    # state looks like x1,v1, x2,v2, x3,v3,...
    ds[1:2:end] .= s[2:2:end] # velocity changes position
    ds[2:2:end] .= -1.0 # (constant downward acceleration)
end

function splitCheck(s, t, intgr)
    #If any of the position coordinates are negative, we need to bounce (and resize).
    if any(s[1:2:end] .< 0.0)
        return true
    else
        return false
    end
end

function splitMod!(intgr)
    s = intgr.u

    # flip the velocity sign
    for i in 1:2:length(s)
        if s[i] < 0.0
            s[i] = 0.0
            s[i + 1] = -s[i + 1]
        end
    end

    # Add a particle to the system.
    # comment out these lines and it will work with Rosenbrock32.
    resize!(intgr, length(s) + 2) # (resizes s -> intgr.u)
    s[end - 1] = rand() # new position
    s[end] = rand() # new velocity
end

function runSim(method)
    s0 = rand(2)
    tspan = (0.0, 20.0)
    prob = ODEProblem(dsdt, s0, tspan)

    # callback to bounce / split system.
    cb = DiscreteCallback(splitCheck, splitMod!)

    solve(prob, method, callback = cb, dtmax = 0.01)
    # setting dtmax here so the discrete callback doesn't miss the zero-crossing too badly.
    # ...no real reason not to use a continuous callback here, I just chose not to.
end

runSim(BS3())

runSim(Rosenbrock23())
runSim(Rosenbrock23(autodiff = false))

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1990
@testset "resize! with SplitODEProblem" begin
    f!(du, u, p, t) = du .= u
    ode = SplitODEProblem(f!, f!, [1.0], (0.0, 1.0))
    integrator = init(ode, Tsit5())
    @test_nowarn step!(integrator)
    @test_nowarn resize!(integrator, 2)
    @test_nowarn step!(integrator)
end
