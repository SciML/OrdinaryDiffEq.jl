@inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::FunctionMapCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  nothing
end

@inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::FunctionMapConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  nothing
end

@inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::Union{SSPRK22ConstantCache,SSPRK33ConstantCache,SSPRK432ConstantCache},always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<1 || calcVal
    copyat_or_push!(k,1,f(uprev,p,t))
  end
  nothing
end

@inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::Union{SSPRK22Cache,SSPRK33Cache,SSPRK432Cache},always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<1 || calcVal
    rtmp = similar(u, eltype(eltype(k)))
    f(rtmp,uprev,p,t)
    copyat_or_push!(k,1,rtmp)
  end
  nothing
end

#=
@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::OwrenZen4Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack k1,k2,k3,k4,k5,k6,tmp = cache
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4 = cache.tab
    a = dt*a21
    @. tmp = uprev+a*k1
    f(k2,tmp,p,t+c1*dt)
    @. tmp = uprev+dt*(a31*k1+a32*k2)
    f(k3,tmp,p,t+c2*dt)
    @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(k4,tmp,p,t+c3*dt)
    @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(k5,tmp,p,t+c4*dt)
    @. u = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
    f(k6,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
  end
  nothing
end

@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::OwrenZen5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,tmp = cache
    @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6 = cache.tab
    a = dt*a21
    @. tmp = uprev+a*k1
    f(k2,tmp,p,t+c1*dt)
    @. tmp = uprev+dt*(a31*k1+a32*k2)
    f(k3,tmp,p,t+c2*dt)
    @. tmp = uprev+dt*(a41*k1+a42*k2+k3)
    f(k4,tmp,p,t+c3*dt)
    @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(k5,tmp,p,t+c4*dt)
    @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(k6,tmp,p,t+dt)
    @. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    f(k7,tmp,p,t+c6*dt)
    @. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
    f(k8,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
  end
  nothing
end

@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::DP5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,c1,c2,c3,c4,c5,c6 = cache.tab
    @unpack d1,d3,d4,d5,d6,d7 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,tmp = cache
    f(k1,uprev,p,t)
    @. tmp = uprev+dt*(a21*k1)
    f(k2,tmp,p,t+c1*dt)
    @. tmp = uprev+dt*(a31*k1+a32*k2)
    f(k3,tmp,p,t+c2*dt)
    @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(k4,tmp,p,t+c3*dt)
    @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(k5,tmp,p,t+c4*dt)
    @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(k6,tmp,p,t+dt)
    @. update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
    @. tmp = uprev+dt*update
    f(k7,tmp,p,t+dt)
    copyat_or_push!(k,1,update)
    @. bspl = k1 - update
    @. dense_tmp3 = update - k7 - bspl
    @. dense_tmp4 = d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,dense_tmp3)
    copyat_or_push!(k,4,dense_tmp4)
  end
  nothing
end

@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::DP5ThreadedCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,c1,c2,c3,c4,c5,c6 = cache.tab
    @unpack d1,d3,d4,d5,d6,d7 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,tmp = cache
    uidx = eachindex(uprev)
    f(k1,uprev,p,t)
    @tight_loop_macros for i in uidx
      tmp = uprev+dt*(a21*k1)
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      tmp = uprev+dt*(a31*k1+a32*k2)
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      tmp =uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    end
    f(k6,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
      tmp = uprev+dt*update
    end
    f(k7,tmp,p,t+dt)
    copyat_or_push!(k,1,update)
    @tight_loop_macros for i in uidx
      bspl = k1 - update
      dense_tmp3 = update - k7 - bspl
      dense_tmp4 = (d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7)
    end
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,dense_tmp3)
    copyat_or_push!(k,4,dense_tmp4)
  end
  nothing
end

@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::Tsit5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<7 || calcVal
    @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,tmp = cache
    @. tmp = uprev+dt*(a21*k1)
    f(k2,tmp,p,t+c1*dt)
    @. tmp = uprev+dt*(a31*k1+a32*k2)
    f(k3,tmp,p,t+c2*dt)
    @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(k4,tmp,p,t+c3*dt)
    @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(k5,tmp,p,t+c4*dt)
    @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(k6,tmp,p,t+dt)
    @. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    f(k7,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
  end
  nothing
end

"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
@muladd @inline function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,p,cache::BS5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k) < 8 || calcVal
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,tmp = cache
    @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87 = cache.tab
    @. tmp = uprev+dt*a21*k1
    f(k2,tmp,p,t+c1*dt)
    @. tmp = uprev+dt*(a31*k1+a32*k2)
    f(k3,tmp,p,t+c2*dt)
    @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(k4,tmp,p,t+c3*dt)
    @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(k5,tmp,p,t+c4*dt)
    @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(k6,tmp,p,t+c5*dt)
    @. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    f(k7,tmp,p,t+dt)
    @. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
    f(k8,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
    copyat_or_push!(k,8,k8)
  end
  if (calcVal2 && length(k)< 11) || calcVal3 # Have not added the extra stages yet
    rtmp = similar(cache.k1)
    @unpack tmp = cache
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache.tab
    @. tmp = uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8])
    f(rtmp,tmp,p,t+c6*dt); copyat_or_push!(k,9,rtmp)
    @. tmp = uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9])
    f(rtmp,tmp,p,t+c7*dt); copyat_or_push!(k,10,rtmp)
    @. tmp = uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10])
    f(rtmp,tmp,p,t+c8*dt); copyat_or_push!(k,11,rtmp)
  end
  nothing
end
=#

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::OwrenZen3ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,c1,c2 = cache
    k1 = f(uprev,p,t)
    a1 = dt*a21
    k2 = f( uprev+a1*k1,p,t+c1*dt)
    tmp = uprev+ dt*(a31*k1 + a32*k2)
    k3 = f(tmp,p,t+c2*dt)
    u = uprev+dt*(a41*k1+a42*k2+a43*k3)
    k4 = f(u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::OwrenZen3Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<4 || calcVal
    @unpack k1,k2,k3,k4,tmp = cache
    @unpack a21,a31,a32,a41,a42,a43,c1,c2 = cache.tab
    # NOTE: k1 does not need to be evaluated since it is aliased with integrator.fsalfirst.
    a1 = dt*a21
    @. tmp = uprev+a1*k1
    f(k2,tmp,p,t+c1*dt)
    @. tmp = uprev+dt*(a31*k1+a32*k2)
    f(k3,tmp,p,t+c2*dt)
    # NOTE: We should not change u here.
    @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(k4,tmp,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::OwrenZen4ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<6 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4 = cache
    k1 = f(uprev,p,t)
    a = dt*a21
    k2 = f( uprev+a*k1,p,t+c1*dt)
    k3 = f( uprev+dt*(a31*k1+a32*k2),p,t+c2*dt)
    k4 = f( uprev+dt*(a41*k1+a42*k2+a43*k3),p,t+c3*dt)
    k5 = f( uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4),p,t+c4*dt)
    u = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
    k6 = f(u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::OwrenZen4Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<6 || calcVal
    uidx = eachindex(uprev)
    @unpack k1,k2,k3,k4,k5,k6,tmp = cache
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4 = cache.tab
    # NOTE: k1 does not need to be evaluated since it is aliased with integrator.fsalfirst.
    a = dt*a21
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    # NOTE: We should not change u here.
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::OwrenZen5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<8 || calcVal
    @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6 = cache
    k1 = f(uprev,p,t)
    a = dt*a21
    k2 = f( uprev+a*k1,p,t+c1*dt)
    k3 = f( uprev+dt*(a31*k1+a32*k2),p,t+c2*dt)
    k4 = f( uprev+dt*(a41*k1+a42*k2+k3),p,t+c3*dt)
    k5 = f( uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4),p,t+c4*dt)
    k6 = f( uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5),p,t+c5*dt)
    k7 = f( uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6),p,t+c6*dt)
    u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
    k8 = f(u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
    copyat_or_push!(k,8,k8)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::OwrenZen5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<8 || calcVal
    uidx = eachindex(uprev)
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,tmp = cache
    @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6 = cache.tab
    # NOTE: k1 does not need to be evaluated since it is aliased with integrator.fsalfirst.
    a = dt*a21
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+c5*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    f(k7,tmp,p,t+c6*dt)
    # NOTE: We should not change u here.
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
    end
    f(k8,tmp,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
    copyat_or_push!(k,8,k8)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::DP5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,c1,c2,c3,c4,c5,c6 = cache
    @unpack d1,d3,d4,d5,d6,d7 = cache
    k1 = f(uprev,p,t)
    k2 = f(uprev+dt*(a21*k1),p,t+c1*dt)
    k3 = f(uprev+dt*(a31*k1+a32*k2),p,t+c2*dt)
    k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3),p,t+c3*dt)
    k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4),p,t+c4*dt)
    k6 = f(uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5),p,t+dt)
    update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
    k7 = f(uprev+dt*update,p,t+dt)
    copyat_or_push!(k,1,update)
    bspl = k1 - update
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,update - k7 - bspl)
    copyat_or_push!(k,4,d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::DP5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,c1,c2,c3,c4,c5,c6 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
    @unpack d1,d3,d4,d5,d6,d7 = cache.tab
    uidx = eachindex(uprev)
    f(k1,uprev,p,t)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a21*k1[i])
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      tmp[i] =uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds update[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
      @inbounds tmp[i] = uprev[i]+dt*update[i]
    end
    f(k7,tmp,p,t+dt)
    copyat_or_push!(k,1,update)
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    @tight_loop_macros for i in uidx
      #integrator.k[4] == k5
      @inbounds k5[i] = d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i]
      #bspl == k3
      @inbounds bspl[i] = k1[i] - update[i]
      # k6 === integrator.k[3] === k2
      @inbounds k6[i] = update[i] - k7[i] - bspl[i]
    end
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,k6)
    copyat_or_push!(k,4,k5)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::DP5ThreadedCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,c1,c2,c3,c4,c5,c6 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
    @unpack d1,d3,d4,d5,d6,d7 = cache.tab
    uidx = eachindex(uprev)
    f(k1,uprev,p,t)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a21*k1[i])
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      tmp[i] =uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds update[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
      @inbounds u[i] = uprev[i]+dt*update[i]
    end
    f(k7,tmp,p,t+dt)
    copyat_or_push!(k,1,update)
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    @tight_loop_macros for i in uidx
      #integrator.k[4] == k5
      @inbounds k5[i] = d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i]
      #bspl == k3
      @inbounds bspl[i] = k1[i] - update[i]
      # k6 === integrator.k[3] === k2
      @inbounds k6[i] = update[i] - k7[i] - bspl[i]
    end
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,k6)
    copyat_or_push!(k,4,k5)
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::Tsit5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<7 || calcVal
    @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76 = cache
    copyat_or_push!(k,1,f(uprev,p,t))
    copyat_or_push!(k,2,f(uprev+dt*(a21*k[1]),p,t+c1*dt))
    copyat_or_push!(k,3,f(uprev+dt*(a31*k[1]+a32*k[2]),p,t+c2*dt))
    copyat_or_push!(k,4,f(uprev+dt*(a41*k[1]+a42*k[2]+a43*k[3]),p,t+c3*dt))
    copyat_or_push!(k,5,f(uprev+dt*(a51*k[1]+a52*k[2]+a53*k[3]+a54*k[4]),p,t+c4*dt))
    copyat_or_push!(k,6,f(uprev+dt*(a61*k[1]+a62*k[2]+a63*k[3]+a64*k[4]+a65*k[5]),p,t+dt))
    utmp = uprev+dt*(a71*k[1]+a72*k[2]+a73*k[3]+a74*k[4]+a75*k[5]+a76*k[6])
    copyat_or_push!(k,7,f(utmp,p,t+dt))
  end
  nothing
end

@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::Tsit5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k)<7 || calcVal
    @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,tmp = cache
    uidx = eachindex(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a21*k1[i])
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    f(k7,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
  end
  nothing
end

"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::BS5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k) < 8 || calcVal
    @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87 = cache
    copyat_or_push!(k,1,f(uprev,p,t))
    copyat_or_push!(k,2,f(uprev+dt*a21*k[1],p,t+c1*dt))
    copyat_or_push!(k,3,f(uprev+dt*(a31*k[1]+a32*k[2]),p,t+c2*dt))
    copyat_or_push!(k,4,f(uprev+dt*(a41*k[1]+a42*k[2]+a43*k[3]),p,t+c3*dt))
    copyat_or_push!(k,5,f(uprev+dt*(a51*k[1]+a52*k[2]+a53*k[3]+a54*k[4]),p,t+c4*dt))
    copyat_or_push!(k,6,f(uprev+dt*(a61*k[1]+a62*k[2]+a63*k[3]+a64*k[4]+a65*k[5]),p,t+c5*dt))
    copyat_or_push!(k,7,f(uprev+dt*(a71*k[1]+a72*k[2]+a73*k[3]+a74*k[4]+a75*k[5]+a76*k[6]),p,t+dt))
    copyat_or_push!(k,8,f(uprev+dt*(a81*k[1]+a83*k[3]+a84*k[4]+a85*k[5]+a86*k[6]+a87*k[7]),p,t+dt))
  end
  if (calcVal2 && length(k)< 11) || calcVal3 # Have not added the extra stages yet
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache
    copyat_or_push!(k,9,f(uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8]),p,t+c6*dt))
    copyat_or_push!(k,10,f(uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9]),p,t+c7*dt))
    copyat_or_push!(k,11,f(uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10]),p,t+c8*dt))
  end
  nothing
end

"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
@muladd @inline function ode_addsteps!(k,t,uprev,u,dt,f,p,cache::BS5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false}) where {calcVal,calcVal2,calcVal3}
  if length(k) < 8 || calcVal
    uidx = eachindex(uprev)
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,tmp = cache
    @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87 = cache.tab
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*a21*k1[i]
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+c5*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    f(k7,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds u[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
    end
    f(k8,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
    copyat_or_push!(k,8,k8)
  end
  if (calcVal2 && length(k)< 11) || calcVal3 # Have not added the extra stages yet
    uidx = eachindex(uprev)
    rtmp = similar(cache.k1)
    @unpack tmp = cache
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache.tab
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a91*k[1][i]+a92*k[2][i]+a93*k[3][i]+a94*k[4][i]+a95*k[5][i]+a96*k[6][i]+a97*k[7][i]+a98*k[8][i])
    end
    f(rtmp,tmp,p,t+c6*dt); copyat_or_push!(k,9,rtmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a101*k[1][i]+a102*k[2][i]+a103*k[3][i]+a104*k[4][i]+a105*k[5][i]+a106*k[6][i]+a107*k[7][i]+a108*k[8][i]+a109*k[9][i])
    end
    f(rtmp,tmp,p,t+c7*dt); copyat_or_push!(k,10,rtmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a111*k[1][i]+a112*k[2][i]+a113*k[3][i]+a114*k[4][i]+a115*k[5][i]+a116*k[6][i]+a117*k[7][i]+a118*k[8][i]+a119*k[9][i]+a1110*k[10][i])
    end
    f(rtmp,tmp,p,t+c8*dt); copyat_or_push!(k,11,rtmp)
  end
  nothing
end
