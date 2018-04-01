function calc_tderivative!(integrator, cache, dtd1, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du2,fsalfirst,dT,tf,linsolve_tmp = cache

    # Time derivative
    if !repeat_step # skip calculation if step is repeated
      if has_tgrad(f)
        f(Val{:tgrad}, dT, uprev, p, t)
      else
        tf.uprev = uprev
        derivative!(dT, tf, t, du2, integrator)
      end
    end

    f(fsalfirst, uprev, p, t)
    @. linsolve_tmp = fsalfirst + dtd1*dT
  end
end

function calc_J(integrator, cache, is_compos)
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du1,uf,J,jac_config = cache
    if has_jac(f)
      f(Val{:jac}, J, uprev, p, t)
    else
      uf.t = t
      jacobian!(J, uf, uprev, du1, integrator, jac_config)
      if is_compos
        integrator.eigen_est = norm(J, Inf)
      end
    end
end

function calc_W!(integrator, cache, dtgamma, repeat_step, check_newW=false, W_transform=false)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p,alg = integrator
    @unpack J,W,jac_config = cache
    mass_matrix = integrator.sol.prob.mass_matrix
    is_compos = typeof(integrator.alg) <: CompositeAlgorithm

    # calculate W
    if has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && W_transform ? f(Val{:invW_t}, W, uprev, p, dtgamma, t) :
                                    f(Val{:invW}, W, uprev, p, dtgamma, t) # W == inverse W
      is_compos && calc_J(integrator, cache, true)

    else
      # skip calculation of J if step is repeated
      if repeat_step || (alg_can_repeat_jac(alg) &&
                         (!integrator.last_stepfail && cache.newton_iters == 1 &&
                          cache.ηold < integrator.alg.new_jac_conv_bound))
        new_jac = false
      else
        new_jac = true
        calc_J(integrator, cache, is_compos)
      end
      # skip calculation of W if step is repeated
      if !repeat_step && (!alg_can_repeat_jac(alg) ||
                          (integrator.iter < 1 || new_jac ||
                           abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t))))
        new_W = true
        if W_transform
          for j in 1:length(u), i in 1:length(u)
              W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
          end
        else
          for j in 1:length(u), i in 1:length(u)
              W[i,j] = mass_matrix[i,j] - dtgamma*J[i,j]
          end
        end
      else
        new_W = false
      end
    end
    return check_newW ? new_W : nothing
  end
end

function calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, W_transform)
  calc_tderivative!(integrator, cache, dtd1, repeat_step)
  calc_W!(integrator, cache, dtgamma, repeat_step, false, W_transform)
end
