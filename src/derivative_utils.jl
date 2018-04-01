function calc_tderivative!(integrator, cache, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du2,fsalfirst,dT,tf,linsolve_tmp = cache
    d1 = cache.tab.d1
    dtd1 = dt*d1

    # Time derivative
    if !repeat_step # skip calculation if step is repeated
      if has_tgrad(f)
        f(Val{:tgrad}, dT, uprev, p, t)
      else
        tf.uprev = uprev
        derivative!(dT, tf, t, du2, integrator)
      end
    end

    f( fsalfirst,  uprev, p, t)

    @. linsolve_tmp = fsalfirst + dtd1*dT

  end
end

function calc_jacobian!(integrator, cache, repeat_step, W_transform=false)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du1,du2,fsalfirst,dT,J,W,uf,tf,linsolve_tmp,linsolve_tmp_vec,jac_config,vectmp = cache
    d1 = cache.tab.d1
    dtgamma = dt*d1
    if !W_transform
      gamma = cache.tab.gamma
      dtgamma = dt*gamma
    end
    mass_matrix = integrator.sol.prob.mass_matrix

    # Jacobian
    if has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && f(Val{:invW}, W, uprev, p, dtgamma, t) # W == inverse W
      if typeof(integrator.alg) <: CompositeAlgorithm
        if has_jac(f)
          f(Val{:jac}, J, uprev, p, t)
        else
          uf.t = t
          jacobian!(J, uf, uprev, du1, integrator, jac_config)
        end
        integrator.eigen_est = norm(J, Inf)
      end

      A_mul_B!(vectmp, W, linsolve_tmp_vec)
    else
      if !repeat_step # skip calculation of J and W if step is repeated
        if has_jac(f)
          f(Val{:jac}, J, uprev, p, t)
        else
          uf.t = t
          jacobian!(J, uf, uprev, du1, integrator, jac_config)
        end
        if typeof(integrator.alg) <: CompositeAlgorithm
          integrator.eigen_est = norm(J, Inf)
        end
        if !W_transform
          for j in 1:length(u), i in 1:length(u)
              W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
          end
        else
          for j in 1:length(u), i in 1:length(u)
              W[i,j] = mass_matrix[i,j] - dtgamma*J[i,j]
          end
        end
      end

      # use existing factorization of W if step is repeated
      cache.linsolve(vectmp, W, linsolve_tmp_vec, !repeat_step)
    end
  end
end

function calc_differentiation!(integrator, cache, repeat_step, W_transform=false)
  calc_tderivative!(integrator, cache, repeat_step)
  calc_jacobian!(integrator, cache, repeat_step, W_transform)
end
