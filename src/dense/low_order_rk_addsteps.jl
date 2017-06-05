function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::DiscreteCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::DiscreteConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::DP5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache
    @unpack d1,d3,d4,d5,d6,d7 = cache
    k1 = f(t,uprev)
    k2 = f(t+c1*dt,@.(@muladd(uprev+dt*(a21*k1))))
    k3 = f(t+c2*dt,@.(@muladd(uprev+dt*(a31*k1+a32*k2))))
    k4 = f(t+c3*dt,@.(@muladd(uprev+dt*(a41*k1+a42*k2+a43*k3))))
    k5 = f(t+c4*dt,@.(@muladd(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))))
    k6 = f(t+dt,@.(@muladd(uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))))
    update = @. @muladd a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
    k7 = f(t+dt,@.(@muladd(uprev+dt*update)))
    copyat_or_push!(k,1,update)
    bspl = k1 - update
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,@.(update - k7 - bspl))
    copyat_or_push!(k,4,@.(@muladd(d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7)))
  end
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::DP5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache.tab
    @unpack d1,d3,d4,d5,d6,d7 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
    f(t,uprev,k1)
    @. tmp = @muladd uprev+dt*(a21*k1)
    f(t+c1*dt,tmp,k2)
    @. tmp = @muladd uprev+dt*(a31*k1+a32*k2)
    f(t+c2*dt,tmp,k3)
    @. tmp = @muladd uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(t+c3*dt,tmp,k4)
    @. tmp = @muladd uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(t+c4*dt,tmp,k5)
    @. tmp = @muladd uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(t+dt,tmp,k6)
    @. update = @muladd a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
    @. tmp = @muladd uprev+dt*update
    f(t+dt,tmp,k7)
    copyat_or_push!(k,1,update)
    @. bspl = k1 - update
    @. dense_tmp3 = update - k7 - bspl
    @. dense_tmp4 = @muladd (d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7)
    copyat_or_push!(k,2,bspl)
    copyat_or_push!(k,3,dense_tmp3)
    copyat_or_push!(k,4,dense_tmp4)
  end
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::DP5ThreadedCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<4 || calcVal
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache.tab
    @unpack d1,d3,d4,d5,d6,d7 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
    uidx = eachindex(uprev)
    f(t,uprev,k1)
    for i in uidx
      tmp = uprev+dt*(a21*k1)
    end
    f(t+c1*dt,tmp,k2)
    for i in uidx
      tmp = uprev+dt*(a31*k1+a32*k2)
    end
    f(t+c2*dt,tmp,k3)
    for i in uidx
      tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
    end
    f(t+c3*dt,tmp,k4)
    for i in uidx
      tmp =uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    end
    f(t+c4*dt,tmp,k5)
    for i in uidx
      tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    end
    f(t+dt,tmp,k6)
    for i in uidx
      update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
      tmp = uprev+dt*update
    end
    f(t+dt,tmp,k7)
    copyat_or_push!(k,1,update)
    for i in uidx
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

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::Tsit5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<7 || calcVal
    @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = cache
    copyat_or_push!(k,1,f(t,uprev))
    copyat_or_push!(k,2,f(t+c1*dt,@.(@muladd(uprev+dt*(a21*k[1])))))
    copyat_or_push!(k,3,f(t+c2*dt,@.(@muladd(uprev+dt*(a31*k[1]+a32*k[2])))))
    copyat_or_push!(k,4,f(t+c3*dt,@.(@muladd(uprev+dt*(a41*k[1]+a42*k[2]+a43*k[3])))))
    copyat_or_push!(k,5,f(t+c4*dt,@.(@muladd(uprev+dt*(a51*k[1]+a52*k[2]+a53*k[3]+a54*k[4])))))
    copyat_or_push!(k,6,f(t+dt,@.(@muladd(uprev+dt*(a61*k[1]+a62*k[2]+a63*k[3]+a64*k[4]+a65*k[5])))))
    utmp = @. @muladd uprev+dt*(a71*k[1]+a72*k[2]+a73*k[3]+a74*k[4]+a75*k[5]+a76*k[6])
    copyat_or_push!(k,7,f(t+dt,utmp))
  end
  nothing
end

function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::Tsit5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k)<7 || calcVal
    @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache
    @. tmp = @muladd uprev+dt*(a21*k1)
    f(t+c1*dt,tmp,k2)
    @. tmp = @muladd uprev+dt*(a31*k1+a32*k2)
    f(t+c2*dt,tmp,k3)
    @. tmp = @muladd uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(t+c3*dt,tmp,k4)
    @. tmp = @muladd uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(t+c4*dt,tmp,k5)
    @. tmp = @muladd uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(t+dt,tmp,k6)
    @. tmp = @muladd uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    f(t+dt,u,k7)
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
function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::BS5ConstantCache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k) < 8 || calcVal
    @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache
    copyat_or_push!(k,1,f(t,uprev))
    copyat_or_push!(k,2,f(t+c1*dt,@.(@muladd(uprev+dt*a21*k[1]))))
    copyat_or_push!(k,3,f(t+c2*dt,@.(@muladd(uprev+dt*(a31*k[1]+a32*k[2])))))
    copyat_or_push!(k,4,f(t+c3*dt,@.(@muladd(uprev+dt*(a41*k[1]+a42*k[2]+a43*k[3])))))
    copyat_or_push!(k,5,f(t+c4*dt,@.(@muladd(uprev+dt*(a51*k[1]+a52*k[2]+a53*k[3]+a54*k[4])))))
    copyat_or_push!(k,6,f(t+c5*dt,@.(@muladd(uprev+dt*(a61*k[1]+a62*k[2]+a63*k[3]+a64*k[4]+a65*k[5])))))
    copyat_or_push!(k,7,f(t+dt,@.(@muladd(uprev+dt*(a71*k[1]+a72*k[2]+a73*k[3]+a74*k[4]+a75*k[5]+a76*k[6])))))
    copyat_or_push!(k,8,f(t+dt,@.(@muladd(uprev+dt*(a81*k[1]+a83*k[3]+a84*k[4]+a85*k[5]+a86*k[6]+a87*k[7])))))
  end
  if (calcVal2 && length(k)< 11) || calcVal3 # Have not added the extra stages yet
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache
    copyat_or_push!(k,9,f(t+c6*dt, @.(@muladd(uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8])))))
    copyat_or_push!(k,10,f(t+c7*dt,@.(@muladd(uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9])))))
    copyat_or_push!(k,11,f(t+c8*dt,@.(@muladd(uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10])))))
  end
  nothing
end

"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
function ode_addsteps!{calcVal,calcVal2,calcVal3}(k,t,uprev,u,dt,f,cache::BS5Cache,always_calc_begin::Type{Val{calcVal}} = Val{false},allow_calc_end::Type{Val{calcVal2}} = Val{true},force_calc_end::Type{Val{calcVal3}} = Val{false})
  if length(k) < 8 || calcVal
    @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde = cache
    @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab
    @. tmp = @muladd uprev+dt*a21*k1
    f(t+c1*dt,tmp,k2)
    @. tmp = @muladd uprev+dt*(a31*k1+a32*k2)
    f(t+c2*dt,tmp,k3)
    @. tmp = @muladd uprev+dt*(a41*k1+a42*k2+a43*k3)
    f(t+c3*dt,tmp,k4)
    @. tmp = @muladd uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    f(t+c4*dt,tmp,k5)
    @. tmp = @muladd uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    f(t+c5*dt,tmp,k6)
    @. tmp = @muladd uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    f(t+dt,tmp,k7)
    @. u = @muladd uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
    f(t+dt,u,k8)
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
    @. tmp = @muladd uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8])
    f(t+c6*dt,tmp,rtmp); copyat_or_push!(k,9,rtmp)
    @. tmp = @muladd uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9])
    f(t+c7*dt,tmp,rtmp); copyat_or_push!(k,10,rtmp)
    @. tmp = @muladd uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10])
    f(t+c8*dt,tmp,rtmp); copyat_or_push!(k,11,rtmp,Val{false})
  end
  nothing
end
