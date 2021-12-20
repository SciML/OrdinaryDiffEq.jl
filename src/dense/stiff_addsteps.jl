function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin
    @unpack tf,uf,d = cache
    γ = dt*d
    tf.u = uprev
    if cache.autodiff
      dT = ForwardDiff.derivative(tf, t)
    else
      dT = FiniteDiff.finite_difference_derivative(tf, t, dir = sign(dt))
    end

    mass_matrix = f.mass_matrix
    if typeof(uprev) <: AbstractArray
      J = ForwardDiff.jacobian(uf, uprev)
      W = mass_matrix - γ*J
    else
      J = ForwardDiff.derivative(uf, uprev)
      W = 1 - γ*J
    end
    k₁ = W\(f(uprev,p,t) + dt*d*dT)
    f₁ = f(uprev+dt*k₁/2,p,t+dt/2)
    k₂ = W\(f₁-k₁) + k₁
    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin
    @unpack k₁,k₂,k₃,du1,du2,f₁,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,weight = cache
    @unpack c₃₂,d = cache.tab
    uidx = eachindex(uprev)

    # Assignments
    sizeu  = size(u)
    mass_matrix = f.mass_matrix
    γ = dt*d
    dto2 = dt/2
    dto6 = dt/6

    @.. linsolve_tmp = @muladd fsalfirst + γ*dT

    ### Jacobian does not need to be re-evaluated after an event
    ### Since it's unchanged
    jacobian2W!(W, mass_matrix, γ, J, false)

    linsolve = cache.linsolve
    linsolve = LinearSolve.set_A(linsolve,W)
    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. veck₁ = -vecu

    @.. tmp = uprev + dto2*k₁
    f(f₁,tmp,p,t+dto2)

    if mass_matrix == I
      tmp .= k₁
    else
      mul!(tmp,mass_matrix,k₁)
    end

    @.. linsolve_tmp = f₁ - tmp

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck2 = _vec(k₂)

    @.. veck2 = -vecu

    @.. k₂ += k₁

    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache::Rosenbrock23Cache{<:Array},always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin
    @unpack k₁,k₂,k₃,du1,du2,f₁,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,weight = cache
    @unpack c₃₂,d = cache.tab
    uidx = eachindex(uprev)

    # Assignments
    sizeu  = size(u)
    mass_matrix = f.mass_matrix
    γ = dt*d
    dto2 = dt/2
    dto6 = dt/6

    @inbounds @simd ivdep for i in eachindex(u)
      linsolve_tmp[i] = @muladd fsalfirst[i] + γ*dT[i]
    end

    ### Jacobian does not need to be re-evaluated after an event
    ### Since it's unchanged
    jacobian2W!(W, mass_matrix, γ, J, false)

    linsolve = cache.linsolve
    linsolve = LinearSolve.set_A(linsolve,W)
    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @inbounds @simd ivdep for i in eachindex(u)
      veck₁[i] = -vecu[i]
      tmp[i] = uprev[i] + dto2*k₁[i]
    end
    f(f₁,tmp,p,t+dto2)

    if mass_matrix == I
      copyto!(tmp,k₁)
    else
      mul!(tmp,mass_matrix,k₁)
    end

    @inbounds @simd ivdep for i in eachindex(u)
      linsolve_tmp[i] = f₁[i] - tmp[i]
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck2 = _vec(k₂)

    @inbounds @simd ivdep for i in eachindex(u)
      veck2[i] = -vecu[i]
      k₂[i] += k₁[i]
    end

    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache::Rodas4ConstantCache,always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin
    @unpack tf,uf = cache
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

    # Precalculations
    dtC21 = C21/dt
    dtC31 = C31/dt
    dtC32 = C32/dt
    dtC41 = C41/dt
    dtC42 = C42/dt
    dtC43 = C43/dt
    dtC51 = C51/dt
    dtC52 = C52/dt
    dtC53 = C53/dt
    dtC54 = C54/dt
    dtC61 = C61/dt
    dtC62 = C62/dt
    dtC63 = C63/dt
    dtC64 = C64/dt
    dtC65 = C65/dt

    dtd1 = dt*d1
    dtd2 = dt*d2
    dtd3 = dt*d3
    dtd4 = dt*d4
    dtgamma = dt*gamma
    mass_matrix = f.mass_matrix

    # Time derivative
    tf.u = uprev
    if cache.autodiff
      dT = ForwardDiff.derivative(tf, t)
    else
      dT = FiniteDiff.finite_difference_derivative(tf, t, dir = sign(dt))
    end

    # Jacobian
    uf.t = t
    if typeof(uprev) <: AbstractArray
      J = ForwardDiff.jacobian(uf, uprev)
      W = mass_matrix/dtgamma - J
    else
      J = ForwardDiff.derivative(uf, uprev)
      W = 1/dtgamma - J
    end

    du = f(uprev, p, t)

    linsolve_tmp =  du + dtd1*dT

    k1 = W\linsolve_tmp
    u = uprev  + a21*k1
    du = f(u, p, t+c2*dt)

    linsolve_tmp =  du + dtd2*dT + dtC21*k1

    k2 = W\linsolve_tmp
    u = uprev  + a31*k1 + a32*k2
    du = f(u, p, t+c3*dt)

    linsolve_tmp =  du + dtd3*dT + (dtC31*k1 + dtC32*k2)

    k3 = W\linsolve_tmp
    u = uprev  + a41*k1 + a42*k2 + a43*k3
    du = f(u, p, t+c4*dt)

    linsolve_tmp =  du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)

    k4 = W\linsolve_tmp
    u = uprev  + a51*k1 + a52*k2 + a53*k3 + a54*k4
    du = f(u, p, t+dt)

    linsolve_tmp =  du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)

    k5 = W\linsolve_tmp

    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    k₁ =  h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    k₂ =  h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache::Rodas4Cache,always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin

    @unpack du,du1,du2,tmp,k1,k2,k3,k4,k5,k6,dT,J,W,uf,tf,linsolve_tmp,jac_config,fsalfirst,weight = cache
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

    # Assignments
    sizeu  = size(u)
    uidx = eachindex(uprev)
    mass_matrix = f.mass_matrix
    atmp = du # does not work with units - additional unitless array required!

    # Precalculations
    dtC21 = C21/dt
    dtC31 = C31/dt
    dtC32 = C32/dt
    dtC41 = C41/dt
    dtC42 = C42/dt
    dtC43 = C43/dt
    dtC51 = C51/dt
    dtC52 = C52/dt
    dtC53 = C53/dt
    dtC54 = C54/dt
    dtC61 = C61/dt
    dtC62 = C62/dt
    dtC63 = C63/dt
    dtC64 = C64/dt
    dtC65 = C65/dt

    dtd1 = dt*d1
    dtd2 = dt*d2
    dtd3 = dt*d3
    dtd4 = dt*d4
    dtgamma = dt*gamma

    @.. linsolve_tmp = @muladd fsalfirst + dtgamma*dT

    ### Jacobian does not need to be re-evaluated after an event
    ### Since it's unchanged
    jacobian2W!(W, mass_matrix, dtgamma, J, true)

    linsolve = cache.linsolve
    linsolve = LinearSolve.set_A(linsolve,W)
    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @.. veck1 = -vecu
    @.. tmp = uprev + a21*k1
    f( du,  tmp, p, t+c2*dt)

    if mass_matrix == I
      @.. linsolve_tmp = du + dtd2*dT + dtC21*k1
    else
      @.. du1 = dtC21*k1
      mul!(du2,mass_matrix,du1)
      @.. linsolve_tmp = du + dtd2*dT + du2
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck2 = _vec(k2)
    @.. veck2 = -vecu
    @.. tmp = uprev + a31*k1 + a32*k2
    f( du,  tmp, p, t+c3*dt)

    if mass_matrix == I
      @.. linsolve_tmp = du + dtd3*dT + (dtC31*k1 + dtC32*k2)
    else
      @.. du1 = dtC31*k1 + dtC32*k2
      mul!(du2,mass_matrix,du1)
      @.. linsolve_tmp = du + dtd3*dT + du2
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck3 = _vec(k3)
    @.. veck3 = -vecu
    @.. tmp = uprev + a41*k1 + a42*k2 + a43*k3
    f( du,  tmp, p, t+c4*dt)

    if mass_matrix == I
      @.. linsolve_tmp = du + dtd4*dT + (dtC41*k1 + dtC42*k2 + dtC43*k3)
    else
      @.. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
      mul!(du2,mass_matrix,du1)
      @.. linsolve_tmp = du + dtd4*dT + du2
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck4 = _vec(k4)
    @.. veck4 = -vecu
    @.. tmp = uprev + a51*k1 + a52*k2 + a53*k3 + a54*k4
    f( du,  tmp, p, t+dt)

    if mass_matrix == I
      @.. linsolve_tmp = du + (dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3)
    else
      @.. du1 = dtC52*k2 + dtC54*k4 + dtC51*k1 + dtC53*k3
      mul!(du2,mass_matrix,du1)
      @.. linsolve_tmp = du + du2
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck5 = _vec(k5)
    @.. veck5 = -vecu
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab
    @.. k6 = h21*k1 + h22*k2 + h23*k3 + h24*k4 + h25*k5
    copyat_or_push!(k,1,copy(k6))

    @.. k6 = h31*k1 + h32*k2 + h33*k3 + h34*k4 + h35*k5
    copyat_or_push!(k,2,copy(k6))

  end
  nothing
end

function DiffEqBase.addsteps!(k,t,uprev,u,dt,f,p,cache::Rodas4Cache{<:Array},always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<2 || always_calc_begin

    @unpack du,du1,du2,tmp,k1,k2,k3,k4,k5,k6,dT,J,W,uf,tf,linsolve_tmp,jac_config,fsalfirst = cache
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,C21,C31,C32,C41,C42,C43,C51,C52,C53,C54,C61,C62,C63,C64,C65,gamma,c2,c3,c4,d1,d2,d3,d4 = cache.tab

    # Assignments
    sizeu  = size(u)
    uidx = eachindex(uprev)
    mass_matrix = f.mass_matrix
    atmp = du # does not work with units - additional unitless array required!

    # Precalculations
    dtC21 = C21/dt
    dtC31 = C31/dt
    dtC32 = C32/dt
    dtC41 = C41/dt
    dtC42 = C42/dt
    dtC43 = C43/dt
    dtC51 = C51/dt
    dtC52 = C52/dt
    dtC53 = C53/dt
    dtC54 = C54/dt
    dtC61 = C61/dt
    dtC62 = C62/dt
    dtC63 = C63/dt
    dtC64 = C64/dt
    dtC65 = C65/dt

    dtd1 = dt*d1
    dtd2 = dt*d2
    dtd3 = dt*d3
    dtd4 = dt*d4
    dtgamma = dt*gamma

    @inbounds @simd ivdep for i in eachindex(u)
      linsolve_tmp[i] = @muladd fsalfirst[i] + dtgamma*dT[i]
    end

    ### Jacobian does not need to be re-evaluated after an event
    ### Since it's unchanged
    jacobian2W!(W, mass_matrix, dtgamma, J, true)

    linsolve = cache.linsolve
    linsolve = LinearSolve.set_A(linsolve,W)
    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @inbounds @simd ivdep for i in eachindex(u)
      veck1[i] = -vecu[i]
    end
    @inbounds @simd ivdep for i in eachindex(u)
      tmp[i] = uprev[i] + a21*k1[i]
    end
    f( du,  tmp, p, t+c2*dt)

    if mass_matrix == I
      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + dtd2*dT[i] + dtC21*k1[i]
      end
    else
      @inbounds @simd ivdep for i in eachindex(u)
        du1[i] = dtC21*k1[i]
      end
      mul!(du2,mass_matrix,du1)

      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + dtd2*dT[i] + du2[i]
      end
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck2 = _vec(k2)
    @inbounds @simd ivdep for i in eachindex(u)
      veck2[i] = -vecu[i]
    end
    @inbounds @simd ivdep for i in eachindex(u)
      tmp[i] = uprev[i] + a31*k1[i] + a32*k2[i]
    end
    f( du,  tmp, p, t+c3*dt)

    if mass_matrix == I
      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + dtd3*dT[i] + (dtC31*k1[i] + dtC32*k2[i])
      end
    else
      @inbounds @simd ivdep for i in eachindex(u)
        du1[i] = dtC31*k1[i] + dtC32*k2[i]
      end
      mul!(du2,mass_matrix,du1)
      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + dtd3*dT[i] + du2[i]
      end
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck3 = _vec(k3)
    @inbounds @simd ivdep for i in eachindex(u)
      veck3[i] = -vecu[i]
      tmp[i] = uprev[i] + a41*k1[i] + a42*k2[i] + a43*k3[i]
    end
    f( du,  tmp, p, t+c4*dt)

    if mass_matrix == I
      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + dtd4*dT[i] + (dtC41*k1[i] + dtC42*k2[i] + dtC43*k3[i])
      end
    else
      @inbounds @simd ivdep for i in eachindex(u)
        du1[i] = dtC41*k1[i] + dtC42*k2[i] + dtC43*k3[i]
      end
      mul!(du2,mass_matrix,du1)

      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + dtd4*dT[i] + du2[i]
      end
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck4 = _vec(k4)
    @inbounds @simd ivdep for i in eachindex(u)
      veck4[i] = -vecu[i]
      tmp[i] = uprev[i] + a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]
    end
    f( du,  tmp, p, t+dt)

    if mass_matrix == I
      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + (dtC52*k2[i] + dtC54*k4[i] + dtC51*k1[i] + dtC53*k3[i])
      end
    else
      @inbounds @simd ivdep for i in eachindex(u)
        du1[i] = dtC52*k2[i] + dtC54*k4[i] + dtC51*k1[i] + dtC53*k3[i]
      end
      mul!(du2,mass_matrix,du1)
      @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = du[i] + du2[i]
      end
    end

    linsolve = LinearSolve.set_b(linsolve,_vec(linsolve_tmp))
    linres = solve(linsolve)
    vecu = _vec(linres.u)
    veck5 = _vec(k5)
    @.. veck5 = -vecu
    @inbounds @simd ivdep for i in eachindex(u)
      veck5[i] = -vecu[i]
    end
    @unpack h21,h22,h23,h24,h25,h31,h32,h33,h34,h35 = cache.tab

    @inbounds @simd ivdep for i in eachindex(u)
      k6[i] = h21*k1[i] + h22*k2[i] + h23*k3[i] + h24*k4[i] + h25*k5[i]
    end
    copyat_or_push!(k,1,copy(k6))

    @inbounds @simd ivdep for i in eachindex(u)
      k6[i] = h31*k1[i] + h32*k2[i] + h33*k3[i] + h34*k4[i] + h35*k5[i]
    end
    copyat_or_push!(k,2,copy(k6))

  end
  nothing
end
