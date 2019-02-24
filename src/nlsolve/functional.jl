"""
    nlsolve!(nlsolve, nlcache::Union{NLFunctionalCache,NLAndersonCache,NLFunctionalConstantCache,NLAnNLAndersonConstantCache}, integrator)

Perform functional iteration that is used by implicit methods.

It solves

```math
G(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅h)
z = G(z)
```

by iterating

```math
zᵏ⁺¹ = G(zᵏ),
```

where `dt` is the step size and `γ` is a constant.

It returns the tuple `(z, η, iter, fail_convergence)`, where `z` is the solution, `η` is
used to measure the iteration error (see [^HW96]), `iter` is the number of iterations, and
`fail_convergence` reports whether the algorithm succeeded.

[^HW96]: Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
"""
@muladd function nlsolve!(nlsolver::NLSolver, nlcache::Union{NLFunctionalConstantCache,NLAndersonConstantCache}, integrator)
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,tmp,κtol,c,γ,max_iter = nlsolver

  if nlcache isa NLAndersonConstantCache
    @unpack Δz₊s,Q,R,γs,aa_start,droptol = nlcache
  end

  # precalculations
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  tstep = t + c*dt
  η = max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
  if nlcache isa NLAndersonConstantCache
    history = 0
    max_history = length(Δz₊s)
  end

  # fixed point iteration
  local ndz
  if nlcache isa NLAndersonConstantCache
    local Δdz, dzold, z₊old # cache variables
  end
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1

    # evaluate function
    u = @. tmp + γ*z
    if mass_matrix == I
      z₊ = dt .* f(u, p, tstep)
      dz = z₊ .- z
    else
      mz = _reshape(mass_matrix * _vec(z), axes(z))
      dz = dt .* f(u, p, tstep) .- mz
      z₊ = z .+ dz
    end

    # compute norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = integrator.opts.internalnorm(dz, tstep)

    # check divergence (after first iteration)
    if iter > 1
      θ = ndz / ndzprev
      if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
        # fixed-point iteration diverges
        break
      end
    end

    # update iterate
    z = z₊

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κtol && (iter > 1 || !iszero(integrator.success_iter))
      # fixed-point iteration converges
      fail_convergence = false
      break
    end

    # perform Anderson acceleration
    if nlcache isa NLAndersonConstantCache && iter < max_iter
      if iter == aa_start
        # update cached values for next step of Anderson acceleration
        dzold = dz
        z₊old = z₊
      elseif iter > aa_start
        # increase size of history
        history += 1

        # remove oldest history if maximum size is exceeded
        if history > max_history
          # circularly shift differences of z₊
          for i in 1:(max_history-1)
            Δz₊s[i] = Δz₊s[i + 1]
          end

          # delete left-most column of QR decomposition
          qrdelete!(Q, R, max_history)

          # update size of history
          history = max_history
        end

        # update history of differences of z₊
        Δz₊s[history] = @. z₊ - z₊old

        # replace/add difference of residuals as right-most column to QR decomposition
        qradd!(Q, R, vec(dz .- dzold), history)

        # update cached values
        dzold = dz
        z₊old = z₊

        # define current Q and R matrices
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

        # check condition (TODO: incremental estimation)
        if droptol !== nothing
          while cond(R) > droptol && history > 1
            qrdelete!(Q, R, history)
            history -= 1
            Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
          end
        end

        # solve least squares problem
        γscur = view(γs, 1:history)
        ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
        integrator.destats.nsolve += 1

        # update next iterate
        for i in 1:history
          z = @. z - γs[i] * Δz₊s[i]
        end

        # update norm of residuals
        ndz = integrator.opts.internalnorm(z .- z₊ .+ dz, tstep)
      end
    end
  end

  integrator.force_stepfail = fail_convergence
  z, η, iter, fail_convergence
end

@muladd function nlsolve!(nlsolver::NLSolver, nlcache::Union{NLFunctionalCache,NLAndersonCache}, integrator)
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,dz,tmp,ztmp,k,κtol,c,γ,max_iter = nlsolver

  if nlcache isa NLFunctionalCache
    @unpack z₊ = nlcache
  else
    @unpack z₊,dzold,z₊old,Δz₊s,Q,R,γs,aa_start,droptol = nlcache
  end

  # precalculations
  vecztmp = vec(ztmp); vecz = vec(z); vecz₊ = vec(z₊)
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  tstep = t + c*dt
  η = max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
  if nlcache isa NLAndersonCache
    history = 0
    max_history = length(Δz₊s)
  end

  # fixed-point iteration without Newton
  local ndz
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1

    # evaluate function
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. z₊ = dt*k
      @. dz = z₊ - z
    else
      mul!(vecztmp, mass_matrix, vecz)
      @. dz = dt*k - ztmp
      @. z₊ = z + dz
    end

    # compute norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = integrator.opts.internalnorm(dz, tstep)

    # check divergence (after first iteration)
    if iter > 1
      θ = ndz / ndzprev
      if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
        # fixed-point iteration diverges
        break
      end
    end

    # update iterate
    @. z = z₊

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κtol && (iter > 1 || !iszero(integrator.success_iter))
      # fixed-point iteration converges
      fail_convergence = false
      break
    end

    # perform Anderson acceleration
    if nlcache isa NLAndersonCache && iter < max_iter
      if iter == aa_start
        # update cached values for next step of Anderson acceleration
        @. dzold = dz
        @. z₊old = z₊
      elseif iter > aa_start
        # increase size of history
        history += 1

        # remove oldest history if maximum size is exceeded
        if history > max_history
          # circularly shift differences of z₊
          ptr = Δgs[1]
          for i in 1:(max_history-1)
            Δz₊s[i] = Δz₊s[i + 1]
          end
          Δz₊s[max_history] = ptr

          # delete left-most column of QR decomposition
          qrdelete!(Q, R, max_history)

          # update size of history
          history = max_history
        end

        # update history of differences of z₊
        @. Δz₊s[history] = z₊ - z₊old

        # replace/add difference of residuals as right-most column to QR decomposition
        @. dzold = dz - dzold
        qradd!(Q, R, vec(dzold), history)

        # update cached values
        @. dzold = dz
        @. z₊old = z₊

        # define current Q and R matrices
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

        # check condition (TODO: incremental estimation)
        if droptol !== nothing
          while cond(R) > droptol && history > 1
            qrdelete!(Q, R, history)
            history -= 1
            Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
          end
        end

        # solve least squares problem
        γscur = view(γs, 1:history)
        ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
        integrator.destats.nsolve += 1

        # update next iterate
        for i in 1:history
          @. z = z - γs[i] * Δz₊s[i]
        end

        # update norm of residuals
        @. dz = z - z₊ + dz
        ndz = integrator.opts.internalnorm(dz, tstep)
      end
    end
  end

  integrator.force_stepfail = fail_convergence
  z, η, iter, fail_convergence
end
