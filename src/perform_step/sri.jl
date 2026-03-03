#=
@muladd function perform_step!(integrator,cache::SRICache)
  (; c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,error_terms) = cache.tab
  (; H0,H1,A0temp,A1temp,B0temp,B1temp,A0temp2,A1temp2,B0temp2,B1temp2,atemp,btemp,E₁,E₂,E₁temp,ftemp,gtemp,chi1,chi2,chi3,tmp) = cache
  (; t,dt,uprev,u,W,p,f) = integrator
  @.. chi1 = .5*(W.dW.^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
  @.. chi2 = .5*(W.dW + W.dZ/sqrt3) #I_(1,0)/h
  @.. chi3 = 1/6 * (W.dW.^3 - 3*W.dW*dt)/dt #I_(1,1,1)/h
  for i=1:stages
    fill!(H0[i],zero(eltype(integrator.u)))
    fill!(H1[i],zero(eltype(integrator.u)))
  end
  for i = 1:stages
    fill!(A0temp,zero(eltype(integrator.u)))
    fill!(B0temp,zero(eltype(integrator.u)))
    fill!(A1temp,zero(eltype(integrator.u)))
    fill!(B1temp,zero(eltype(integrator.u)))
    for j = 1:i-1
      integrator.f((t + c₀[j]*dt),H0[j],ftemp)
      integrator.f.g((t + c₁[j]*dt),H1[j],gtemp)
      @.. A0temp = A0temp + A₀[j,i] * ftemp
      @.. B0temp = B0temp + B₀[j,i] * gtemp
      @.. A1temp = A1temp + A₁[j,i] * ftemp
      @.. B1temp = B1temp + B₁[j,i] * gtemp
    end
    @.. H0[i] = uprev + A0temp*dt + B0temp*chi2
    @.. H1[i] = uprev + A1temp*dt + B1temp*integrator.sqdt
  end
  fill!(atemp,zero(eltype(integrator.u)))
  fill!(btemp,zero(eltype(integrator.u)))
  fill!(E₂,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))
  for i = 1:stages
    integrator.f((t+c₀[i]*dt),H0[i],ftemp)
    integrator.f.g((t+c₁[i]*dt),H1[i],gtemp)
    @.. atemp = atemp + α[i] * ftemp
    @.. btemp = btemp + (β₁[i] * W.dW + β₂[i] * chi1) * gtemp
    @.. E₂    = E₂    + (β₃[i] * chi2 + β₄[i] * chi3) * gtemp
    if i <= error_terms
      @.. E₁temp += ftemp
    end
  end

  @.. u = uprev + (dt*atemp + btemp) + E₂

  if integrator.opts.adaptive
    @.. E₁ = dt * E₁temp

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
  integrator.u = u
end
=#

@muladd function perform_step!(integrator, cache::SRICache)
    (; c₀, c₁, A₀, A₁, B₀, B₁, α, β₁, β₂, β₃, β₄, stages, error_terms) = cache.tab
    (;
        H0, H1, A0temp, A1temp, B0temp, B1temp, A0temp2, A1temp2, B0temp2, B1temp2,
        atemp, btemp, E₁, E₂, E₁temp, ftemp, gtemp, chi1, chi2, chi3, tmp,
    ) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    sqrt3 = sqrt(3one(eltype(W.dW)))
    if W.dW isa Union{SArray, Number}
        chi1 = (W.dW .^ 2 - abs(dt)) / 2integrator.sqdt #I_(1,1)/sqrt(h)
        chi2 = (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
        chi3 = (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    else
        @.. chi1 = (W.dW .^ 2 - dt) / 2integrator.sqdt #I_(1,1)/sqrt(h)
        @.. chi2 = (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
        @.. chi3 = (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    end

    for i in 1:stages
        fill!(H0[i], zero(eltype(integrator.u)))
        fill!(H1[i], zero(eltype(integrator.u)))
    end
    for i in 1:stages
        fill!(A0temp, zero(eltype(integrator.u)))
        fill!(B0temp, zero(eltype(integrator.u)))
        fill!(A1temp, zero(eltype(integrator.u)))
        fill!(B1temp, zero(eltype(integrator.u)))
        for j in 1:(i - 1)
            integrator.f(ftemp, H0[j], p, t + c₀[j] * dt)
            integrator.f.g(gtemp, H1[j], p, t + c₁[j] * dt)
            @.. A0temp = A0temp + A₀[j, i] * ftemp
            @.. B0temp = B0temp + B₀[j, i] * gtemp
            @.. A1temp = A1temp + A₁[j, i] * ftemp
            @.. B1temp = B1temp + B₁[j, i] * gtemp
        end
        @.. H0[i] = uprev + A0temp * dt + B0temp * chi2
        @.. H1[i] = uprev + A1temp * dt + B1temp * integrator.sqdt
    end
    fill!(atemp, zero(eltype(integrator.u)))
    fill!(btemp, zero(eltype(integrator.u)))
    fill!(E₂, zero(eltype(integrator.u)))
    fill!(E₁temp, zero(eltype(integrator.u)))
    for i in 1:stages
        integrator.f(ftemp, H0[i], p, t + c₀[i] * dt)
        integrator.f.g(gtemp, H1[i], p, t + c₁[i] * dt)
        @.. atemp = atemp + α[i] * ftemp
        @.. btemp = btemp + (β₁[i] * W.dW + β₂[i] * chi1) * gtemp
        @.. E₂ = E₂ + (β₃[i] * chi2 + β₄[i] * chi3) * gtemp
        if i <= error_terms
            @.. E₁temp += ftemp
        end
    end

    @.. u = uprev + (dt * atemp + btemp) + E₂

    if integrator.opts.adaptive
        @.. E₁ = dt * E₁temp

        calculate_residuals!(
            tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
end

#=
@muladd function perform_step!(integrator,cache::SRIW1Cache)
  (; t,dt,uprev,u,W,p,f) = integrator
  (; chi1,chi2,chi3,fH01o4,g₁o2,H0,H11,H12,H13,g₂o3,Fg₂o3,g₃o3,Tg₃o3,mg₁,E₁,E₂,fH01,fH02,g₁,g₂,g₃,g₄,tmp) = cache
  @.. chi1 = (W.dW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
  @.. chi2 = (W.dW + W.dZ/sqrt3)/2 #I_(1,0)/h
  @.. chi3 = (W.dW.^3 - 3W.dW*dt)/6dt #I_(1,1,1)/h
  integrator.f(fH01,uprev,p,t)
  @.. fH01 = dt*fH01
  integrator.f.g(t,uprev,g₁)
  dto4 = dt/4
  @.. fH01o4 = fH01/4
  @.. g₁o2 = g₁/2
  @.. H0 =  uprev + 3 * (fH01o4 + chi2 * g₁o2)
  @.. H11 = uprev + fH01o4 + integrator.sqdt * g₁o2
  @.. H12 = uprev + fH01 - integrator.sqdt * g₁
  integrator.f.g(t+dto4,H11,g₂)
  integrator.f.g(t+dt,H12,g₃)
  @.. H13 = uprev + fH01o4 + integrator.sqdt * (-5 * g₁ + 3 * g₂ + g₃ / 2)

  integrator.f.g(t+dto4,H13,g₄)
  integrator.f(fH02,H0,p,t+3dto4)

  @.. fH02 = fH02*dt
  @.. g₂o3 = g₂/3
  @.. Fg₂o3 = 4g₂o3
  @.. g₃o3 = g₃/3
  @.. Tg₃o3 = 2g₃o3
  @.. mg₁ = -g₁
  @.. E₂ = chi2 * (2 * g₁ - Fg₂o3 - Tg₃o3) + chi3 * (2 * mg₁ + 5 * g₂o3 - Tg₃o3 + g₄)

  @.. u = uprev +  (fH01 + 2 * fH02) / 3 + W.dW * (mg₁ + Fg₂o3 + Tg₃o3) + chi1 * (mg₁ + Fg₂o3 - g₃o3) + E₂

  if integrator.opts.adaptive
    @.. E₁ = fH01 + fH02

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
  integrator.u = u
end
=#

@muladd function perform_step!(integrator, cache::SRIW1Cache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (;
        chi1, chi2, chi3, fH01o4, g₁o2, H0, H11, H12, H13, g₂o3, Fg₂o3, g₃o3,
        Tg₃o3, mg₁, E₁, E₂, fH01, fH02, g₁, g₂, g₃, g₄, tmp,
    ) = cache

    sqrt3 = sqrt(3one(eltype(W.dW)))
    if W.dW isa Union{SArray, Number}
        chi1 = (W.dW .^ 2 - abs(dt)) / 2integrator.sqdt #I_(1,1)/sqrt(h)
        chi2 = (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
        chi3 = (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    else
        @.. chi1 = (W.dW .^ 2 - dt) / 2integrator.sqdt #I_(1,1)/sqrt(h)
        @.. chi2 = (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
        @.. chi3 = (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    end

    integrator.f(fH01, uprev, p, t)
    @.. fH01 = dt * fH01
    integrator.f.g(g₁, uprev, p, t)
    dto4 = dt / 4

    @.. fH01o4 = fH01 / 4
    @.. g₁o2 = g₁ / 2
    @.. H0 = uprev + 3 * (fH01o4 + chi2 * g₁o2)
    @.. H11 = uprev + fH01o4 + integrator.sqdt * g₁o2
    @.. H12 = uprev + fH01 - integrator.sqdt * g₁

    integrator.f.g(g₂, H11, p, t + dto4)
    integrator.f.g(g₃, H12, p, t + dt)

    @.. H13 = uprev + fH01o4 + integrator.sqdt * (-5g₁ + 3g₂ + g₃ / 2)

    integrator.f.g(g₄, H13, p, t + dto4)
    integrator.f(fH02, H0, p, t + 3dto4)

    @.. fH02 = fH02 * dt
    @.. g₂o3 = g₂ / 3
    @.. Fg₂o3 = 4g₂o3
    @.. g₃o3 = g₃ / 3
    @.. Tg₃o3 = 2g₃o3
    @.. mg₁ = -g₁

    @.. E₂ = chi2 * (2 * g₁ - Fg₂o3 - Tg₃o3) + chi3 * (2 * mg₁ + 5 * g₂o3 - Tg₃o3 + g₄)

    @.. u = uprev + (fH01 + 2 * fH02) / 3 + W.dW * (mg₁ + Fg₂o3 + Tg₃o3) + chi1 * (mg₁ + Fg₂o3 - g₃o3) + E₂

    if integrator.opts.adaptive
        @.. E₁ = fH01 + fH02

        calculate_residuals!(
            tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
end

@muladd function perform_step!(integrator, cache::SRIW1ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    sqrt3 = sqrt(3one(eltype(W.dW)))
    chi1 = @.. (W.dW .^ 2 - abs(dt)) / 2integrator.sqdt #I_(1,1)/sqrt(h)
    chi2 = @.. (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
    chi3 = @.. (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    fH01 = dt * integrator.f(uprev, p, t)

    g₁ = integrator.f.g(uprev, p, t)
    fH01o4 = fH01 / 4
    dto4 = dt / 4
    g₁o2 = g₁ / 2
    H0 = @.. uprev + 3 * (fH01o4 + chi2 * g₁o2)
    H11 = @.. uprev + fH01o4 + integrator.sqdt * g₁o2
    H12 = @.. uprev + fH01 - integrator.sqdt * g₁
    g₂ = integrator.f.g(H11, p, t + dto4)
    g₃ = integrator.f.g(H12, p, t + dt)
    H13 = @.. uprev + fH01o4 + integrator.sqdt * (-5 * g₁ + 3 * g₂ + g₃ / 2)

    g₄ = integrator.f.g(H13, p, t + dto4)
    fH02 = dt * integrator.f(H0, p, t + 3dto4)

    g₂o3 = g₂ / 3
    Fg₂o3 = 4g₂o3
    g₃o3 = g₃ / 3
    Tg₃o3 = 2g₃o3
    mg₁ = -g₁
    E₂ = @.. chi2 * (2 * g₁ - Fg₂o3 - Tg₃o3) + chi3 * (2 * mg₁ + 5 * g₂o3 - Tg₃o3 + g₄)

    u = uprev + (fH01 + 2fH02) / 3 + W.dW .* (mg₁ + Fg₂o3 + Tg₃o3) + chi1 .* (mg₁ + Fg₂o3 - g₃o3) + E₂
    if integrator.opts.adaptive
        E₁ = fH01 .+ fH02

        resids = calculate_residuals(
            E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(resids, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SRIConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (; c₀, c₁, A₀, A₁, B₀, B₁, α, β₁, β₂, β₃, β₄, stages, H0, H1, error_terms) = cache
    sqrt3 = sqrt(3one(eltype(W.dW)))
    chi1 = 0.5 * (W.dW .^ 2 - abs(dt)) / integrator.sqdt #I_(1,1)/sqrt(h)
    chi2 = 0.5 * (W.dW + W.dZ / sqrt3) #I_(1,0)/h
    chi3 = 1 / 6 * (W.dW .^ 3 - 3 * W.dW * dt) / dt #I_(1,1,1)/h

    fill!(H0, zero(u))
    fill!(H1, zero(u))
    @inbounds for i in 1:stages
        A0temp = zero(u)
        B0temp = zero(u)
        A1temp = zero(u)
        B1temp = zero(u)
        for j in 1:(i - 1)
            ftmp = integrator.f(H0[j], p, t + c₀[j] * dt)
            A0temp = @.. A0temp + A₀[j, i] * ftmp
            A1temp = @.. A1temp + A₁[j, i] * ftmp

            gtmp = integrator.f.g(H1[j], p, t + c₁[j] * dt)
            B0temp = @.. B0temp + B₀[j, i] * gtmp
            B1temp = @.. B1temp + B₁[j, i] * gtmp
        end
        H0[i] = @.. uprev + dt * A0temp + chi2 * B0temp
        H1[i] = @.. uprev + dt * A1temp + integrator.sqdt * B1temp
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂ = zero(u)
    E₁temp = zero(u)
    @inbounds for i in 1:stages
        ftmp = integrator.f(H0[i], p, t + c₀[i] * dt)
        atemp = @.. atemp + α[i] * ftmp

        gtmp = integrator.f.g(H1[i], p, t + c₁[i] * dt)
        btemp = @.. btemp + (β₁[i] * W.dW + β₂[i] * chi1) * gtmp
        E₂ = @.. E₂ + (β₃[i] * chi2 + β₄[i] * chi3) * gtmp
        if i <= error_terms #1 or 2
            E₁temp = E₁temp .+ ftmp
        end
    end

    u = (uprev + dt * atemp) + btemp + E₂

    if integrator.opts.adaptive
        E₁ = dt .* E₁temp

        resids = calculate_residuals(
            E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(resids, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::FourStageSRIConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (;
        a021,
        a031, a032, a041, a042, a043, a121, a131, a132, a141, a142, a143, b021, b031, b032,
        b041, b042, b043, b121, b131, b132, b141, b142, b143, c02, c03, c04, c11, c12,
        c13, c14, α1, α2, α3, α4, beta11, beta12, beta13, beta14, beta21, beta22, beta23,
        beta24, beta31, beta32, beta33, beta34, beta41, beta42, beta43, beta44,
    ) = cache

    sqrt3 = sqrt(3one(eltype(W.dW)))
    chi1 = (W.dW .^ 2 .- abs(dt)) / (2integrator.sqdt) #I_(1,1)/sqrt(h)
    chi2 = (W.dW .+ W.dZ ./ sqrt3) ./ 2 #I_(1,0)/h
    chi3 = (W.dW .^ 3 .- 3 * W.dW * dt) / (6dt) #I_(1,1,1)/h

    k1 = integrator.f(uprev, p, t)
    g1 = integrator.f.g(uprev, p, t + c11 * dt)

    H01 = uprev + dt * a021 * k1 + b021 * chi2 .* g1
    H11 = uprev + dt * a121 * k1 + integrator.sqdt * b121 * g1

    k2 = integrator.f(H01, p, t + c02 * dt)
    g2 = integrator.f.g(H11, p, t + c12 * dt)

    H02 = uprev + dt * (a031 * k1 + a032 * k2) + chi2 .* (b031 * g1 + b032 * g2)
    H12 = uprev + dt * (a131 * k1 + a132 * k2) + integrator.sqdt * (b131 * g1 + b132 * g2)

    k3 = integrator.f(H02, p, t + c03 * dt)
    g3 = integrator.f.g(H12, p, t + c13 * dt)

    H03 = uprev + dt * (a041 * k1 + a042 * k2 + a043 * k3) + chi2 .* (b041 * g1 + b042 * g2 + b043 * g3)
    H13 = uprev + dt * (a141 * k1 + a142 * k2 + a143 * k3) + integrator.sqdt * (b141 * g1 + b142 * g2 + b143 * g3)

    k4 = integrator.f(H03, p, t + c04 * dt)
    g4 = integrator.f.g(H13, p, t + c14 * dt)

    E₂ = chi2 .* (beta31 * g1 + beta32 * g2 + beta33 * g3 + beta34 * g4) + chi3 .* (beta41 * g1 + beta42 * g2 + beta43 * g3 + beta44 * g4)

    u = uprev + dt * (α1 * k1 + α2 * k2 + α3 * k3 + α4 * k4) + E₂ + W.dW .* (beta11 * g1 + beta12 * g2 + beta13 * g3 + beta14 * g4) + chi1 .* (beta21 * g1 + beta22 * g2 + beta23 * g3 + beta24 * g4)

    E₁ = dt * (k1 + k2 + k3 + k4)

    if integrator.alg isa StochasticCompositeAlgorithm && integrator.alg.algs[1] isa SOSRI2
        ϱu = integrator.opts.internalnorm(k4 - k3, t)
        ϱd = integrator.opts.internalnorm(H03 - H02, t)
        integrator.eigen_est = ϱu / ϱd
    end

    if integrator.opts.adaptive
        E₁ = @.. dt * (k1 + k2 + k3 + k4)

        resids = calculate_residuals(
            E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(resids, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::FourStageSRICache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (;
        chi1, chi2, chi3, tab, g1, g2, g3, g4, k1, k2, k3, k4, E₁, E₂, tmp, H02,
        H03,
    ) = cache
    (;
        a021,
        a031, a032, a041, a042, a043, a121, a131, a132, a141, a142, a143, b021, b031, b032,
        b041, b042, b043, b121, b131, b132, b141, b142, b143, c02, c03, c04, c11, c12,
        c13, c14, α1, α2, α3, α4, beta11, beta12, beta13, beta14, beta21, beta22, beta23,
        beta24, beta31, beta32, beta33, beta34, beta41, beta42, beta43, beta44,
    ) = cache.tab

    sqdt = integrator.sqdt
    sqrt3 = sqrt(3one(eltype(W.dW)))
    if W.dW isa Union{SArray, Number}
        chi1 = (W.dW .^ 2 - abs(dt)) / 2integrator.sqdt #I_(1,1)/sqrt(h)
        chi2 = (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
        chi3 = (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    else
        @.. chi1 = (W.dW .^ 2 - dt) / 2integrator.sqdt #I_(1,1)/sqrt(h)
        @.. chi2 = (W.dW + W.dZ / sqrt3) / 2 #I_(1,0)/h
        @.. chi3 = (W.dW .^ 3 - 3W.dW * dt) / 6dt #I_(1,1,1)/h
    end

    integrator.f(k1, uprev, p, t)
    integrator.f.g(g1, uprev, p, t + c11 * dt)

    @.. tmp = uprev + dt * a021 * k1 + chi2 * b021 * g1
    integrator.f(k2, tmp, p, t + c02 * dt)

    @.. tmp = uprev + dt * a121 * k1 + sqdt * b121 * g1
    integrator.f.g(g2, tmp, p, t + c12 * dt)

    @.. H02 = uprev + dt * (a031 * k1 + a032 * k2) + chi2 * (b031 * g1 + b032 * g2)
    integrator.f(k3, H02, p, t + c03 * dt)
    @.. tmp = uprev + dt * (a131 * k1 + a132 * k2) + sqdt * (b131 * g1 + b132 * g2)
    integrator.f.g(g3, tmp, p, t + c13 * dt)

    @.. H03 = uprev + dt * (a041 * k1 + a042 * k2 + a043 * k3) + chi2 * (b041 * g1 + b042 * g2 + b043 * g3)
    integrator.f(k4, H03, p, t + c04 * dt)

    @.. tmp = uprev + dt * (a141 * k1 + a142 * k2 + a143 * k3) + sqdt * (b141 * g1 + b142 * g2 + b143 * g3)
    integrator.f.g(g4, tmp, p, t + c14 * dt)

    if integrator.alg isa StochasticCompositeAlgorithm && integrator.alg.algs[1] isa SOSRI2
        @.. tmp = k4 - k3
        ϱu = integrator.opts.internalnorm(tmp, t)
        @.. tmp = H03 - H02
        ϱd = integrator.opts.internalnorm(tmp, t)
        integrator.eigen_est = ϱu / ϱd
    end

    @.. E₂ = chi2 * (beta31 * g1 + beta32 * g2 + beta33 * g3 + beta34 * g4) + chi3 * (beta41 * g1 + beta42 * g2 + beta43 * g3 + beta44 * g4)
    @.. u = uprev + dt * (α1 * k1 + α2 * k2 + α3 * k3 + α4 * k4) + E₂ + W.dW * (beta11 * g1 + beta12 * g2 + beta13 * g3 + beta14 * g4) + chi1 * (beta21 * g1 + beta22 * g2 + beta23 * g3 + beta24 * g4)

    if integrator.opts.adaptive
        @.. E₁ = dt * (k1 + k2 + k3 + k4)

        calculate_residuals!(
            tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
end
