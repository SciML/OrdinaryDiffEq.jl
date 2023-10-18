function initialize!(integrator, cache::NewmarkCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    duprev, uprev = integrator.uprev.x
    ku = integrator.f(duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.fsalfirst = ku
end

@muladd function perform_step!(integrator, cache::NewmarkCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    @unpack β, γ, nlsolver = cache

    mass_matrix = f.mass_matrix

    nlsolver.z = uprev

    # _tmp = mass_matrix * @.. broadcast=false (α₁ * uprev+α₂ * uprev2)
    # nlsolver.tmp = @.. broadcast=false _tmp/(dt * β₀)

    # u = nlsolve!(nlsolver, integrator, cache, repeat_step)
    # nlsolvefail(nlsolver) && return

    f(integrator.fsallast, u, p, t + dt)
    integrator.stats.nf += 1
end
