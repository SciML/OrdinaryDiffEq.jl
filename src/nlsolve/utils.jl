## Anderson acceleration

"""
    anderson(z, cache, integrator)

Return the next iterate of the fixed-point iteration `z = g(z)` by performing Anderson
acceleration based on the current iterate `z` and the settings and history in the `cache`.
"""
@muladd function anderson(z, cache, integrator)
  @unpack dz,Δgzs,Q,R,γs,history,droptol = cache

  # increase size of history
  history += 1
  
  # remove oldest history if maximum size is exceeded
  max_history = length(Δgzs)
  if history > max_history
    # circularly shift differences of G(z)
    for i in 1:(max_history-1)
      Δgzs[i] = Δgzs[i + 1]
    end
  
    # delete left-most column of QR decomposition
    qrdelete!(Q, R, max_history)
  
    # update size of history
    history = max_history
  end
  
  # update history of differences of G(z)
  Δgzs[history] = @.. z - cache.gzprev
  
  # replace/add difference of residuals as right-most column to QR decomposition
  qradd!(Q, R, _vec(dz .- cache.dzprev), history)
  
  # update cached values
  cache.dzprev = dz
  cache.gzprev = z
  
  # define current Q and R matrices
  Qcur = view(Q, :, 1:history)
  Rcur = UpperTriangular(view(R, 1:history, 1:history))
  
  # check condition (TODO: incremental estimation)
  if droptol !== nothing
    while cond(R) > droptol && history > 1
      qrdelete!(Q, R, history)
      history -= 1
      Qcur = view(Q, :, 1:history)
      Rcur = UpperTriangular(view(R, 1:history, 1:history))
    end
  end
  
  # solve least squares problem
  γscur = view(γs, 1:history)
  ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end
  
  # update iterate
  for i in 1:history
    z = @.. z - γs[i] * Δgzs[i]
  end
  
  # update cached values
  cache.history = history

  z
end

"""
    anderson!(z, cache, integrator)

Update the current iterate `z` of the fixed-point iteration `z = g(z)` in-place
by performing Anderson acceleration based on the settings and history in the `cache`.
"""
@muladd function anderson!(z, cache, integrator)
  @unpack gzprev,dz,dzprev,Δgzs,Q,R,γs,history,droptol = cache

  # increase size of history
  history += 1

  # remove oldest history if maximum size is exceeded
  max_history = length(Δgzs)
  if history > max_history
    # circularly shift differences of z
    ptr = Δgzs[1]
    for i in 1:(max_history-1)
      Δgzs[i] = Δgzs[i + 1]
    end
    Δgzs[max_history] = ptr

    # delete left-most column of QR decomposition
    qrdelete!(Q, R, max_history)

    # update size of history
    history = max_history
  end

  # update history of differences of g(z)
  @.. Δgzs[history] = z - gzprev

  # replace/add difference of residuals as right-most column to QR decomposition
  @.. dzprev = dz - dzprev
  qradd!(Q, R, vec(dzprev), history)

  # update cached values
  @.. dzprev = dz
  @.. gzprev = z

  # define current Q and R matrices
  Qcur = view(Q, :, 1:history)
  Rcur = UpperTriangular(view(R, 1:history, 1:history))

  # check condition (TODO: incremental estimation)
  if droptol !== nothing
    while cond(R) > droptol && history > 1
      qrdelete!(Q, R, history)
      history -= 1
      Qcur = view(Q, :, 1:history)
      Rcur = UpperTriangular(view(R, 1:history, 1:history))
    end
  end

  # solve least squares problem
  γscur = view(γs, 1:history)
  ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  # update next iterate
  for i in 1:history
    @.. z = z - γs[i] * Δgzs[i]
  end

  # update cached values
  cache.history = history

  nothing
end
