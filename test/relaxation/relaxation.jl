using Roots
using LinearAlgebra
using UnPack

struct Relaxation{OPT, INV, FSAL}
    opt::OPT
    invariant::INV
    has_fsal::Bool
    update_fsal::FSAL

    Relaxation(opt, inv, fsal = nothing) = new{typeof(opt), typeof(inv), typeof(fsal)}(opt, inv, fsal isa Nothing ? false : true, fsal)
end

function (r::Relaxation)(integrator)

    @unpack t, dt, uprev, u = integrator

    # We fix here the bounds of interval where we are going to look for the relaxation
    #(gamma_min, gamma_max) = apriori_bounds_dt(integrator) ./ dt
    
    S_u = u - uprev

    # Find relaxation paramter gamma
    gamma_min = 0.5
    gamma_max = 1.5

    if (r.invariant(gamma_min*S_u .+ uprev) .- r.invariant(uprev)) * (r.invariant(gamma_max*S_u .+ uprev) .- r.invariant(uprev)) ≤ 0
        gamma = find_zero(gamma -> r.invariant(gamma*S_u .+ uprev) .- r.invariant(uprev),
                            (gamma_min, gamma_max),
                            r.opt())

        change_dt!(integrator, gamma*dt)
        change_u!(integrator, uprev + gamma*S_u)

        if r.has_fsal
            integrator.fsallast = r.update_fsal(gamma, integrator.fsalfirst, integrator.fsallast)
        end
    end
end


fsal_r(gamma, fsalfirst, fsallast) = fsalfirst + gamma * (fsallast - fsalfirst)
r_fsal(gamma, fsalfirst, fsallast) = fsalfirst + 1/gamma * (fsallast - fsalfirst)
fsal_r_int(integrator) = integrator.fsalllast = fsal_r(gamma, fsalfirst, fsallast)

