# Skip Enzyme tests on Julia 1.12+ prerelease versions
if !isempty(VERSION.prerelease)
    @warn "Skipping Enzyme tests on Julia prerelease version $(VERSION)"
    exit(0)
end

using Enzyme, OrdinaryDiffEqTsit5, StaticArrays, DiffEqBase, ForwardDiff, Test

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

const _saveat = SA[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

function f_dt(y::Array{Float64}, u0::Array{Float64})
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    y .= sol[1, :]
    return nothing
end;

function f_dt(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    return sol[1, :]
end;

u0 = [1.0; 0.0; 0.0]
fdj = ForwardDiff.jacobian(f_dt, u0)

ezj = stack(
    map(1:3) do i
        d_u0 = zeros(3)
        dy = zeros(13)
        y = zeros(13)
        d_u0[i] = 1.0
        Enzyme.autodiff(Forward, f_dt, Duplicated(y, dy), Duplicated(u0, d_u0))
        dy
    end
)

@test ezj ≈ fdj

function f_dt2(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), dt = 0.1, saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    return sum(sol[1, :])
end

fdg = ForwardDiff.gradient(f_dt2, u0)
d_u0 = zeros(3)
Enzyme.autodiff(Reverse, f_dt2, Active, Duplicated(u0, d_u0));

@test d_u0 ≈ fdg
