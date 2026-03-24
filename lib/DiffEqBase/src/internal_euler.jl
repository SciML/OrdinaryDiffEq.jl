# This is a simple example showing how forward and backward Euler
# could be wrapped

module InternalEuler

using DiffEqBase, LinearAlgebra

# make a algorithm type
abstract type EulerAlgs <: DiffEqBase.AbstractODEAlgorithm end

struct FwdEulerAlg <: EulerAlgs end

struct BwdEulerAlg <: EulerAlgs end

function DiffEqBase.solve(
        prob::DiffEqBase.AbstractODEProblem{uType, tType, isinplace},
        Alg::FwdEulerAlg;
        dt = (prob.tspan[2] - prob.tspan[1]) / 100,
        tstops = tType[],
        kwargs...
    ) where {uType, tType, isinplace}
    u0 = prob.u0
    f = prob.f
    tspan = prob.tspan
    p = prob.p

    if isempty(tstops)
        tstops = tspan[1]:dt:tspan[2]
    end
    @assert tstops[1] == tspan[1]

    nt = length(tstops)
    out = Vector{uType}(undef, nt)
    out[1] = copy(u0)
    tmp = copy(u0)
    for i in 2:nt
        t = tstops[i]
        dt = t - tstops[i - 1]
        if isinplace
            f(tmp, out[i - 1], p, t)
            out[i] = out[i - 1] + dt * tmp
        else
            out[i] = out[i - 1] + dt * f(out[i - 1], p, t)
        end
    end
    # make solution type
    return DiffEqBase.build_solution(prob, Alg, tstops, out)
end

function DiffEqBase.solve(
        prob::DiffEqBase.AbstractODEProblem{uType, tType, isinplace},
        Alg::BwdEulerAlg;
        dt = (prob.tspan[2] - prob.tspan[1]) / 100,
        tstops = tType[],
        tol = 1.0e-5,
        maxiter = 100,
        kwargs...
    ) where {uType, tType, isinplace}
    u0 = prob.u0
    f = prob.f
    tspan = prob.tspan
    p = prob.p
    # TODO: fix numparameters as it picks up the Jacobian
    #    @assert !isinplace "Only out of place functions supported"
    @assert DiffEqBase.has_jac(f) "Provide Jacobian as f(::Val{:jac}, ...)"

    if isempty(tstops)
        tstops = tspan[1]:dt:tspan[2]
    end
    @assert tstops[1] == tspan[1]

    nt = length(tstops)
    out = Vector{uType}(undef, nt)
    out[1] = copy(u0)
    for i in 2:nt
        t = tstops[i]
        dt = t - tstops[i - 1]
        out[i] = newton(t, dt, out[i - 1], p, f, f.jac, tol, maxiter)
    end
    # make solution type
    return DiffEqBase.build_solution(prob, Alg, tstops, out)
end

function newton(t, dt, u_last, p, f, jac, tol, maxiter)
    res = (u) -> u - u_last - dt * f(u, p, t)
    u = u_last + dt * f(u_last, p, t) # forward Euler step as first guess
    for i in 1:maxiter
        du = -(I - dt * jac(u, p, t)) \ res(u)
        u += du
        norm(du, Inf) < tol && return u
    end
    error("Newton not converged")
end

end
