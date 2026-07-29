function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

function lorenz_oop(u, p, t)
    return [10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end

# Parameterized variant for the AutoDePSpecialize workloads: the opaque-p path
# needs an isbits non-NullParameters `p` to trigger, and the precompiled
# solve is shared across *all* isbits parameter types since the wrapped
# signature carries `OpaqueParams` instead of `typeof(p)`.
function lorenz_p(du, u, p, t)
    du[1] = p.σ * (u[2] - u[1])
    du[2] = u[1] * (p.ρ - u[3]) - u[2]
    return du[3] = u[1] * u[2] - p.β * u[3]
end

const lorenz_p_params = (σ = 10.0, ρ = 28.0, β = 8 / 3)

# Non-isbits variant for the AutoDePSpecialize workloads: a non-`isbits` `p`
# (here a `Vector`) is packed by reference into an `OpaqueRef`, and — like the
# isbits `OpaqueParams` case — the precompiled solve is shared across *all*
# non-isbits parameter types since the wrapped signature carries `OpaqueRef`
# instead of `typeof(p)`. Without this the `OpaqueRef` path compiles on first
# use and non-isbits `p` sees no TTFX benefit.
function lorenz_pref(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    return du[3] = u[1] * u[2] - p[3] * u[3]
end

const lorenz_pref_params = [10.0, 28.0, 8 / 3]

PrecompileTools.@compile_workload begin
    ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0))
    ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])
    ODEProblem{true, SciMLBase.AutoSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0)
    )
    ODEProblem{true, SciMLBase.AutoSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0), Float64[]
    )
    ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0)
    )
    ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0), Float64[]
    )
    ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0))
    ODEProblem{true, SciMLBase.NoSpecialize}(
        lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
        Float64[]
    )
    ODEProblem{true, AutoDePSpecialize}(
        lorenz_p, [1.0; 0.0; 0.0],
        (0.0, 1.0), lorenz_p_params
    )
    ODEProblem{true, AutoDePSpecialize}(
        lorenz_pref, [1.0; 0.0; 0.0],
        (0.0, 1.0), lorenz_pref_params
    )

    lorenz([1.0; 0.0; 0.0], [1.0; 0.0; 0.0], SciMLBase.NullParameters(), 0.0)
    lorenz([1.0; 0.0; 0.0], [1.0; 0.0; 0.0], Float64[], 0.0)
    lorenz_p([1.0; 0.0; 0.0], [1.0; 0.0; 0.0], lorenz_p_params, 0.0)
    lorenz_pref([1.0; 0.0; 0.0], [1.0; 0.0; 0.0], lorenz_pref_params, 0.0)
    lorenz_oop([1.0; 0.0; 0.0], SciMLBase.NullParameters(), 0.0)
    lorenz_oop([1.0; 0.0; 0.0], Float64[], 0.0)
end
