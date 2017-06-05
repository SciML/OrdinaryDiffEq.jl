function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::Rosenbrock23ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<2 || calcVal
    dT = ForwardDiff.derivative(tf,t)
    J = ForwardDiff.derivative(uf,uprev)
    W = 1-dt*d*J
    k₁ = W\(f(t,uprev) + dt*d*dT)
    f₁ = f(t+dt/2,uprev+dt*k₁/2)
    k₂ = W\(f₁-k₁) + k₁
    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::Rosenbrock32ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<2 || calcVal
    dT = ForwardDiff.derivative(tf,t)
    J = ForwardDiff.derivative(uf,uprev)
    W = 1-dt*d*J
    k₁ = W\(f(t,uprev) .+ dt.*d.*dT)
    f₁ = f(t+dt/2,uprev.+dt.*k₁/2)
    k₂ = W\(f₁-k₁) + k₁
    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::Rosenbrock23Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<2 || calcVal
    @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec = cache
    @unpack c₃₂,d = cache.tab
    uidx = eachindex(uprev)

    γ = dt*d

    @. linsolve_tmp = @muladd fsalfirst + γ*dT

    if has_invW(f)
      f(Val{:invW},t,u,γ,W) # W == inverse W
      A_mul_B!(vectmp,W,linsolve_tmp_vec)
    else
      ### Jacobian does not need to be re-evaluated after an event
      ### Since it's unchanged
      for i in 1:length(u), j in 1:length(u)
        W[i,j] = @muladd I[i,j]-γ*J[i,j]
      end
      cache.linsolve(vectmp,W,linsolve_tmp_vec,true)
    end

    recursivecopy!(k₁,reshape(vectmp,size(u)...))
    @. tmp=uprev+dt*k₁/2
    f(t+dt/2,tmp,f₁)

    @. linsolve_tmp = f₁-k₁

    if has_invW(f)
      A_mul_B!(vectmp2,W,linsolve_tmp_vec)
    else
      cache.linsolve(vectmp2,W,linsolve_tmp_vec)
    end

    @. k₂ = tmp + k₁
    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::Rosenbrock32Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<2 || calcVal
    @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,uf,tf,linsolve_tmp,linsolve_tmp_vec = cache
    @unpack c₃₂,d = cache.tab
    uidx = eachindex(uprev)

    γ = dt*d

    @. linsolve_tmp = @muladd fsalfirst + γ*dT

    if has_invW(f)
      f(Val{:invW},t,u,γ,W) # W == inverse W
      A_mul_B!(vectmp,W,linsolve_tmp_vec)
    else
      ### Jacobian does not need to be re-evaluated after an event
      ### Since it's unchanged
      for i in 1:length(u), j in 1:length(u)
        W[i,j] = @muladd I[i,j]-γ*J[i,j]
      end
      cache.linsolve(vectmp,W,linsolve_tmp_vec,true)
    end

    recursivecopy!(k₁,reshape(vectmp,size(u)...))
    @. tmp=uprev+dt*k₁/2
    f(t+dt/2,tmp,f₁)
    @. linsolve_tmp = f₁-k₁

    if has_invW(f)
      A_mul_B!(vectmp2,W,linsolve_tmp_vec)
    else
      cache.linsolve(vectmp2,W,linsolve_tmp_vec)
    end

    @. k₂ = tmp + k₁
    copyat_or_push!(k,1,k₁)
    copyat_or_push!(k,2,k₂)
  end
  nothing
end
