@cache mutable struct AitkenNevilleCache{uType,rateType,arrayType,dtType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  utilde::uType
  atmp::uNoUnitsType
  fsalfirst::rateType
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int
  u_tmps::Array{uType,1}
  k_tmps::Array{rateType,1}
end

@cache mutable struct AitkenNevilleConstantCache{dtType,arrayType} <: OrdinaryDiffEqConstantCache
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int
end

function alg_cache(alg::AitkenNeville,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  utilde = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  cur_order = max(alg.init_order, alg.min_order)
  dtpropose = zero(dt)
  T = Array{typeof(u),2}(undef, alg.max_order, alg.max_order)
  # Array of arrays of length equal to number of threads to store intermediate
  # values of u and k. [Thread Safety]
  u_tmps = Array{typeof(u),1}(undef, Threads.nthreads())
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  # Initialize each element of u_tmps and k_tmps to different instance of
  # zeros array similar to u and k respectively
  for i=1:Threads.nthreads()
      u_tmps[i] = zero(u)
      k_tmps[i] = zero(rate_prototype)
  end
  # Initialize lower triangle of T to different instance of zeros array similar to u
  for i=1:alg.max_order
    for j=1:i
      T[i,j] = zero(u)
    end
  end
  work = zero(dt)
  A = one(Int)
  atmp = similar(u,uEltypeNoUnits)
  step_no = zero(Int)
  AitkenNevilleCache(u,uprev,tmp,k,utilde,atmp,fsalfirst,dtpropose,T,cur_order,work,A,step_no,u_tmps,k_tmps)
end

function alg_cache(alg::AitkenNeville,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dtpropose = zero(dt)
  cur_order = max(alg.init_order, alg.min_order)
  T = Array{typeof(u),2}(undef, alg.max_order, alg.max_order)
  @.. T = u
  work = zero(dt)
  A = one(Int)
  step_no = zero(Int)
  AitkenNevilleConstantCache(dtpropose,T,cur_order,work,A,step_no)
end

@cache mutable struct ImplicitEulerExtrapolationCache{uType,rateType,arrayType,dtType,JType,WType,F,JCType,GCType,uNoUnitsType,TFType,UFType} <: OrdinaryDiffEqMutableCache
  uprev::uType
  u_tmps::Array{uType,1}
  utilde::uType
  tmp::uType
  atmp::uNoUnitsType
  k_tmps::Array{rateType,1}
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int
  du1::rateType
  du2::rateType
  J::JType
  W::WType
  tf::TFType
  uf::UFType
  linsolve_tmps::Array{rateType,1}
  linsolve::Array{F,1}
  jac_config::JCType
  grad_config::GCType
end

@cache mutable struct ImplicitEulerExtrapolationConstantCache{dtType,arrayType,TF,UF} <: OrdinaryDiffEqConstantCache
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int

  tf::TF
  uf::UF
end

function alg_cache(alg::ImplicitEulerExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dtpropose = zero(dt)
  cur_order = max(alg.init_order, alg.min_order)
  T = Array{typeof(u),2}(undef, alg.max_order, alg.max_order)
  @.. T = u
  work = zero(dt)
  A = one(Int)
  step_no = zero(Int)
  tf = TimeDerivativeWrapper(f,u,p)
  uf = UDerivativeWrapper(f,t,p)
  ImplicitEulerExtrapolationConstantCache(dtpropose,T,cur_order,work,A,step_no,tf,uf)
end

function alg_cache(alg::ImplicitEulerExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  u_tmp = similar(u)
  u_tmps = Array{typeof(u_tmp),1}(undef, Threads.nthreads())

  u_tmps[1] = u_tmp
  for i=2:Threads.nthreads()
    u_tmps[i] = zero(u_tmp)
  end

  utilde = similar(u)
  tmp = similar(u)
  k_tmp = zero(rate_prototype)
  k_tmps = Array{typeof(k_tmp),1}(undef, Threads.nthreads())

  k_tmps[1] = k_tmp
  for i=2:Threads.nthreads()
    k_tmps[i] = zero(rate_prototype)
  end

  cur_order = max(alg.init_order, alg.min_order)
  dtpropose = zero(dt)
  T = Array{typeof(u),2}(undef, alg.max_order, alg.max_order)
  # Initialize lower triangle of T to different instance of zeros array similar to u
  for i=1:alg.max_order
    for j=1:i
      T[i,j] = zero(u)
    end
  end
  work = zero(dt)
  A = one(Int)
  atmp = similar(u,uEltypeNoUnits)
  step_no = zero(Int)

  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)

  if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
    W_el = WOperator(f, dt, true)
    J = nothing # is J = W.J better?
  else
    J = false .* vec(rate_prototype) .* vec(rate_prototype)' # uEltype?
    W_el = similar(J)
  end

  W = Array{typeof(W_el),1}(undef, Threads.nthreads())
  W[1] = W_el
  for i=2:Threads.nthreads()
    if W_el isa WOperator
      W[i] = WOperator(f, dt, true)
    else
      W[i] = zero(W_el)
    end
  end

  tf = TimeGradientWrapper(f,uprev,p)
  uf = UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve_tmps = Array{typeof(linsolve_tmp),1}(undef, Threads.nthreads())

  for i=1:Threads.nthreads()
    linsolve_tmps[i] = zero(rate_prototype)
  end

  linsolve_el = alg.linsolve(Val{:init},uf,u)
  linsolve = Array{typeof(linsolve_el),1}(undef, Threads.nthreads())
  linsolve[1] = linsolve_el
  for i=2:Threads.nthreads()
    linsolve[i] = alg.linsolve(Val{:init},uf,u)
  end
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,du1,du2)


  ImplicitEulerExtrapolationCache(uprev,u_tmps,utilde,tmp,atmp,k_tmps,dtpropose,T,cur_order,work,A,step_no,
    du1,du2,J,W,tf,uf,linsolve_tmps,linsolve,jac_config,grad_config)
end



struct extrapolation_coefficients{T1,T2,T3}
  # This structure is used by the caches of the algorithms
  # ExtrapolationMidpointDeuflhard() and  ExtrapolationMidpointHairerWanner().
  # It contains the constant coefficients used to extrapolate the internal discretisations
  # in their perfom_step! function and some additional constant data.

  subdividing_sequence::T1  # subdividing_sequence[n] is used for the (n -1)th internal discretisation

  # Weights and Scaling factors for extrapolation operators
  extrapolation_weights::T2
  extrapolation_scalars::T3

  # Weights and scaling factors for internal extrapolation operators (used for error estimate)
  extrapolation_weights_2::T2
  extrapolation_scalars_2::T3
end

function create_extrapolation_coefficients(T, alg::Union{ExtrapolationMidpointDeuflhard,
                                                         ExtrapolationMidpointHairerWanner,
                                                         ImplicitDeuflhardExtrapolation,
                                                         ImplicitHairerWannerExtrapolation})
  # Compute and return extrapolation_coefficients

  @unpack n_min, n_init, n_max, sequence = alg

  # Initialize subdividing_sequence:
  if sequence == :harmonic
      subdividing_sequence = BigInt.(1:(n_max + 1))
  elseif sequence == :romberg
      subdividing_sequence = BigInt(2).^(0:n_max)
  else # sequence == :bulirsch
      subdividing_sequence = [n==0 ? BigInt(1) : (isodd(n) ? BigInt(2)^((n + 1) ÷ 2) : 3 * BigInt(2)^(n ÷ 2 - 1)) for n = 0:n_max]
  end


  # Compute nodes corresponding to subdividing_sequence
  nodes = BigInt(1) .// subdividing_sequence .^ 2

  # Compute barycentric weights for internal extrapolation operators
  extrapolation_weights_2 = zeros(Rational{BigInt}, n_max, n_max)
  extrapolation_weights_2[1,:] = ones(Rational{BigInt}, 1, n_max)
  for n = 2:n_max
      distance = nodes[2:n] .- nodes[n+1]
      extrapolation_weights_2[1:(n-1), n] = extrapolation_weights_2[1:n-1, n-1] .// distance
      extrapolation_weights_2[n, n] = 1 // prod(-distance)
  end

  # Compute barycentric weights for extrapolation operators
  extrapolation_weights = zeros(Rational{BigInt}, n_max+1, n_max+1)
  for n = 1:n_max
      extrapolation_weights[n+1, (n+1) : (n_max+1)] = extrapolation_weights_2[n, n:n_max] // (nodes[n+1] - nodes[1])
      extrapolation_weights[1, n] = 1 // prod(nodes[1] .- nodes[2:n])
  end
  extrapolation_weights[1, n_max+1] = 1 // prod(nodes[1] .- nodes[2:n_max+1])

  # Rescale barycentric weights to obtain weights of 1. Barycentric Formula
  for m = 1:(n_max+1)
      extrapolation_weights[1:m, m] = - extrapolation_weights[1:m, m] .// nodes[1:m]
      if 2 <= m
          extrapolation_weights_2[1:m-1, m-1] = -extrapolation_weights_2[1:m-1, m-1] .// nodes[2:m]
      end
  end

  # Compute scaling factors for internal extrapolation operators
  extrapolation_scalars_2 = ones(Rational{BigInt}, n_max)
  extrapolation_scalars_2[1] = -nodes[2]
  for n = 1:(n_max-1)
      extrapolation_scalars_2[n+1] = -extrapolation_scalars_2[n] * nodes[n+2]
  end

  # Compute scaling factors for extrapolation operators
  extrapolation_scalars = -nodes[1] * [BigInt(1); extrapolation_scalars_2]

  # Initialize and return extrapolation_coefficients
  extrapolation_coefficients(Int.(subdividing_sequence),
      T.(extrapolation_weights), T.(extrapolation_scalars),
      T.(extrapolation_weights_2), T.(extrapolation_scalars_2))
end

function create_extrapolation_coefficients(T::Type{<:CompiledFloats}, alg::Union{ExtrapolationMidpointDeuflhard,
                                                                                 ExtrapolationMidpointHairerWanner,
                                                                                 ImplicitDeuflhardExtrapolation,
                                                                                 ImplicitHairerWannerExtrapolation})
  # Compute and return extrapolation_coefficients

  @unpack n_min, n_init, n_max, sequence = alg

  n_max > 15 && error("n_max > 15 not allowed for Float32 or Float64 with this algorithm. That's a bad idea.")

  # Initialize subdividing_sequence:
  if sequence == :harmonic
    subdividing_sequence = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    extrapolation_weights = T[-1.0 -1.3333333333333333 -1.5 -1.6 -1.6666666666666667 -1.7142857142857142 -1.75 -1.7777777777777777 -1.8 -1.8181818181818181 -1.8333333333333333 -1.8461538461538463 -1.8571428571428572 -1.8666666666666667 -1.875 -1.8823529411764706; 0.0 5.333333333333333 38.4 204.8 975.2380952380952 4388.571428571428 19114.666666666668 81555.91111111111 343170.32727272727 1.4298763636363635e6 5.915044102564103e6 2.433618145054945e7 9.970459794285715e7 4.0712710826666665e8 1.6579836988235295e9 6.737203601568627e9; 0.0 0.0 -72.9 -1499.6571428571428 -21088.928571428572 -253067.14285714287 -2.79006525e6 -2.9219592436363637e7 -2.9584837341818184e8 -2.925972923916084e9 -2.8449861733434067e10 -2.7311867264096704e11 -2.596334381793193e12 -2.449162486354648e13 -2.2960898309574828e14 -2.141777720860745e15; 0.0 0.0 0.0 1872.4571428571428 83220.31746031746 2.3967451428571427e6 5.6940854303030305e7 1.214738225131313e9 2.422001138107972e10 4.613335501158042e11 8.50611193356378e12 1.5311001480414803e14 2.7059443139242895e15 4.7143563158147624e16 8.120422362168969e17 1.3858854164768373e19; 0.0 0.0 0.0 0.0 -77504.96031746031 -6.341314935064935e6 -3.236712831439394e8 -1.3278821872571873e10 -4.80171683784965e11 -1.6005722792832168e13 -5.043469942533053e14 -1.525755612867142e16 -4.4766093502525523e17 -1.2827711003647664e19 -3.607793719775906e20 -9.995618963881296e21; 0.0 0.0 0.0 0.0 0.0 4.711650077922078e6 6.393346721118882e8 5.260811016234965e10 3.4090055385202573e12 1.9175656154176447e14 9.826959789128542e15 4.7169406987817e17 2.15773437679608e19 9.515608601670712e20 4.078117972144591e22 1.708360692331116e24; 0.0 0.0 0.0 0.0 0.0 0.0 -3.952348909376457e8 -8.263044119869713e10 -1.0248756909925902e13 -9.846844874242534e14 -8.108603230469998e16 -6.022558357283822e18 -4.1560671463889437e20 -2.715297202307443e22 -1.7009177076954297e24 -1.0307396968759164e26; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.3741255122091064e10 1.3338509797230594e13 2.371290630618772e15 3.221627130440662e17 3.711314454267642e19 3.8230073464151254e21 3.6330154661690406e23 3.249405137443117e25 2.772825717284793e27; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -6.174193142616171e12 -2.6321560239574205e15 -6.449440297701669e17 -1.1940678036887662e20 -1.8574538823517637e22 -2.5642554640188345e24 -3.245385821648838e26 -3.845504022726303e28; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.082508822446903e15 6.237312738860727e17 2.041302350899874e20 4.99971155510259e22 1.0207744425001122e25 1.837393996500202e27 3.015210660923408e29; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.30788366240374e17 -1.748372388422729e20 -7.448430618928414e22 -2.355293074113417e25 -6.165659032955555e27 -1.4147218829987503e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.87960511179021e19 5.723442800021062e22 3.1065086459191243e25 1.2426034583676497e28 4.089940525827236e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.7640007344435631e22 -2.1641022343595773e25 -1.4694640618129094e28 -7.307458985088932e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.155890364150897e24 9.361198795139813e27 7.828458512415587e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.472332709198601e27 -4.593753679027078e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1322357960170301e30]
    extrapolation_weights_2 = T[-4.0 -28.8 -153.6 -731.4285714285714 -3291.4285714285716 -14336.0 -61166.933333333334 -257377.74545454545 -1.0724072727272727e6 -4.436283076923077e6 -1.8252136087912086e7 -7.477844845714286e7 -3.053453312e8 -1.2434877741176472e9 -5.052902701176471e9; 0.0 64.8 1333.0285714285715 18745.714285714286 224948.57142857142 2.480058e6 2.5972971054545455e7 2.6297633192727274e8 2.6008648212587414e9 2.5288765985274727e10 2.4277215345863736e11 2.3078527838161714e12 2.1770333212041316e13 2.0409687386288734e14 1.9038024185428845e15; 0.0 0.0 -1755.4285714285713 -78019.04761904762 -2.2469485714285714e6 -5.338205090909091e7 -1.138817086060606e9 -2.2706260669762238e10 -4.325002032335664e11 -7.974479937716044e12 -1.4354063887888878e14 -2.5368227943040215e15 -4.41970904607634e16 -7.612895964533408e17 -1.299267577947035e19; 0.0 0.0 0.0 74404.76190476191 6.087662337662337e6 3.107244318181818e8 1.2747668997668997e10 4.609648164335664e11 1.536549388111888e13 4.8417311448317306e14 1.4647253883524564e16 4.29754497624245e17 1.2314602563501758e19 3.4634819709848696e20 9.595794205326044e21; 0.0 0.0 0.0 0.0 -4.580770909090909e6 -6.215753756643356e8 -5.114677376895105e10 -3.314310940228028e12 -1.8642999038782656e14 -9.553988683874972e15 -4.585914568259986e17 -2.0977973107739664e19 -9.251286140513193e20 -3.964836917362797e22 -1.6609062286552515e24; 0.0 0.0 0.0 0.0 0.0 3.8716887275524473e8 8.094410566402983e10 1.0039598605641701e13 9.64588885640085e14 7.943121531888978e16 5.899649003053539e18 4.071249449523863e20 2.659882973688924e22 1.6662051014159312e24 1.0097041928580407e26; 0.0 0.0 0.0 0.0 0.0 0.0 -4.3057798010808395e10 -1.3130095581648865e13 -2.3342392145153535e15 -3.1712892065275264e17 -3.65332516591971e19 -3.763272856627389e21 -3.5762495995101494e23 -3.198633182170568e25 -2.729500315452218e27; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.097968535917206e12 2.59966027057523e15 6.369817577976957e17 1.1793262258654482e20 1.8345223529400134e22 2.532597989154405e24 3.205319330023544e26 3.7980286644210397e28; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0716837342224339e15 -6.174939611472119e17 -2.0208893273908753e20 -4.949714439551565e22 -1.010566698075111e25 -1.8190200565352e27 -2.9850585543141743e29; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.2888102437061885e17 1.733923029840723e20 7.38687334108603e22 2.3358278420959505e25 6.114703173179063e27 1.4030299666103307e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -5.838774520736112e19 -5.6836966694653606e22 -3.0849356692113527e25 -1.233974267684541e28 -4.061538161064547e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.753562860275258e22 2.1512968956947278e25 1.4607690081927147e28 7.2642195828103e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -6.124482760252167e24 -9.313437576797262e27 -7.788517397556324e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.461344563824385e27 4.57333699600918e30; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.1278129999388386e30]
    extrapolation_scalars = T[-1.0, 0.25, -0.027777777777777776, 0.001736111111111111, -6.944444444444444e-5, 1.9290123456790124e-6, -3.936759889140842e-8, 6.151187326782565e-10, -7.594058428126624e-12, 7.594058428126623e-14, -6.276081345559193e-16, 4.358389823304995e-18, -2.5789288895295828e-20, 1.3157800456783586e-22, -5.8479113141260385e-25, 2.2843403570804838e-27]
    extrapolation_scalars_2 = T[-0.25, 0.027777777777777776, -0.001736111111111111, 6.944444444444444e-5, -1.9290123456790124e-6, 3.936759889140842e-8, -6.151187326782565e-10, 7.594058428126624e-12, -7.594058428126623e-14, 6.276081345559193e-16, -4.358389823304995e-18, 2.5789288895295828e-20, -1.3157800456783586e-22, 5.8479113141260385e-25, -2.2843403570804838e-27]
  elseif sequence == :romberg
    subdividing_sequence = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]
    extrapolation_weights = T[-1.0 -1.3333333333333333 -1.4222222222222223 -1.4447971781305116 -1.4504630494172979 -1.451880901860521 -1.452235451531305 -1.4523240943593299 -1.4523462554044868 -1.4523517956869105 -1.4523531807588372 -1.4523535270269017 -1.4523536135939228 -1.4523536352356785 -1.4523536406461173 -1.452353641998727; 0.0 5.333333333333333 28.444444444444443 121.36296296296297 493.15743680188126 1980.3655501377505 7929.205565360925 31724.5675171852 126906.01579724405 507631.8090356717 2.0305349820189304e6 8.122147673959353e6 3.248859844172289e7 1.2995440151277749e8 5.19817613796996e8 2.07927046293387e9; 0.0 0.0 -91.02222222222223 -1941.8074074074075 -33140.17975308642 -538659.4296374682 -8.652349112921841e6 -1.3857291091506496e8 -2.2177080072599993e9 -3.54854939788296e10 -5.677765672441477e11 -9.084459730370057e12 -1.4535149430390788e14 -2.325624463334606e15 -3.720999363124215e16 -5.953599069714284e17; 0.0 0.0 0.0 5917.889241622575 504993.2152851264 3.4474203496797964e7 2.241370436871182e9 1.4401024799097037e11 9.225665310201598e12 5.905867660750885e14 3.77998601491761e16 2.4192279640364677e18 1.5483118033241417e20 9.909204991428803e21 6.341892706539483e23 4.05881157410929e25; 0.0 0.0 0.0 0.0 -1.5209207425057925e6 -5.1914094677531046e8 -1.4176008786611145e11 -3.6866623485688414e13 -9.474866810815984e15 -2.4279369357326935e18 -6.217036386624773e20 -1.591658462098873e23 -4.074707838080507e25 -1.0431291857706623e28 -2.670413262277438e30 -6.836259581321538e32; 0.0 0.0 0.0 0.0 0.0 1.5589452477944808e9 2.1284799116553977e12 2.3248676581708025e15 2.4184528070774876e18 2.486207422190278e21 2.5483650380553205e24 2.6101630458060033e27 2.6729701040532997e30 2.7371631523457504e33 2.8028657600863993e36 2.870137275504677e39; 0.0 0.0 0.0 0.0 0.0 0.0 -6.386999060908798e12 -3.4881530871309916e16 -1.5239973381214444e20 -6.341377114357268e23 -2.607614058456583e27 -1.069122783562139e31 -4.380196305334128e34 -1.7942379182565492e38 -7.349310654759861e41 -3.010289127531343e45; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0465098000284594e17 2.2861355418221703e21 3.995311436502874e25 6.649821721972122e29 1.0937793665786103e34 1.7937998718897874e38 2.93967940527153e42 4.816664723480877e46 7.891743901406595e50; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -6.858511278043386e21 -5.993058601571351e26 -4.189451610800854e31 -2.7891799442838834e36 -1.835085270122301e41 -1.2038170852496654e46 -7.891262227584487e50 -5.171933283225826e55; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.7979244390088466e27 6.284201388527134e32 1.7571900680469942e38 4.679485289631606e43 1.2315095760466199e49 3.231484209329825e54 8.473210620642257e59; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.8852622144842938e33 -2.635787615753444e39 -2.948078543974702e45 -3.140352413792322e51 -3.305811498811932e57 -3.46978305819599e63; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 7.907364732522996e39 4.422118870277351e46 1.9784224923818425e53 8.429821331796452e59 3.549588914822444e66; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.3266357401568573e47 -2.9676339154575293e54 -5.310787755579382e61 -9.051452272305827e68; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 8.902901879036164e54 7.966181752074431e62 5.702415016525277e70; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.389854534525231e63 -8.553622556652643e71; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.5660867693856473e72]
    extrapolation_weights_2 = T[-4.0 -21.333333333333332 -91.02222222222223 -369.86807760141096 -1485.274162603313 -5946.904174020694 -23793.4256378889 -95179.51184793304 -380723.8567767538 -1.522901236514198e6 -6.091610755469514e6 -2.4366448831292167e7 -9.746580113458312e7 -3.89863210347747e8 -1.5594528472004025e9; 0.0 85.33333333333333 1820.4444444444443 31068.91851851852 504993.2152851264 8.111577293364226e6 1.299121039828734e8 2.0791012568062494e9 3.3267650605152744e10 5.322905317913885e11 8.516680997221928e12 1.3626702590991364e14 2.1802729343761932e15 3.4884369029289516e16 5.5814991278571405e17; 0.0 0.0 -5825.422222222222 -497102.6962962963 -3.3935544067160495e7 -2.2063490237950697e9 -1.4176008786611145e11 -9.081514289729697e12 -5.813588478551652e14 -3.7209237334345224e16 -2.3814275270983977e18 -1.524119431397202e20 -9.754373663437728e21 -6.242800632999802e23 -3.995392643263833e25; 0.0 0.0 0.0 1.5149796458553793e6 5.1711305245196944e8 1.4120633752288446e11 3.6722613237697445e13 9.437855612336234e15 2.4184528070774876e18 6.19275108823952e20 1.585441046231299e23 4.058791010588005e25 1.0390544623887458e28 2.659981960471667e30 6.809555442332001e32; 0.0 0.0 0.0 0.0 -1.5574228403259315e9 -2.1264013179916716e12 -2.32259727959837e15 -2.416091036758076e18 -2.4837794852545453e21 -2.545876400322845e24 -2.607614058456583e27 -2.6703597816860604e30 -2.7344901414547876e33 -2.8001285864925646e36 -2.867334407071567e39; 0.0 0.0 0.0 0.0 0.0 6.385439734966193e12 3.4873014872562036e16 1.523625268458817e20 6.339828926585209e23 2.606977433930593e27 1.0688617672575583e31 4.379126921470521e34 1.7937998718897874e38 7.34751638946329e41 3.009554193662317e45; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0464459261392974e17 -2.2859960071821667e21 -3.995067582045079e25 -6.649415849064287e29 -1.093712607584068e34 -1.7936903870343255e38 -2.9394999814797044e42 -4.816370737596875e46 -7.891262227584487e50; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.858406625466511e21 5.99296715475431e26 4.1893876848424375e31 2.789137384775456e36 1.8350572689432526e41 1.2037987164586916e46 7.89114181647872e50 5.171854365786813e55; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.7979175804714053e27 -6.284177416201281e32 -1.7571833648988464e38 -4.679467438811868e43 -1.2315048782104076e49 -3.231471882195848e54 -8.47317829790887e59; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.8852604165581403e33 2.6357851020704912e39 2.948075732467912e45 3.140349418918881e51 3.305808346144411e57 3.4697797491530044e63; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -7.907362847260331e39 -4.4221178159620535e46 -1.978422020689163e53 -8.429819321970427e59 -3.549588068534498e66; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.3266356610832054e47 2.967633738572764e54 5.310787439031764e61 9.051451732797232e68; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -8.902901746372587e54 -7.966181633369073e62 -5.702414931552672e70; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.3898545256223295e63 8.553622524787915e71; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.5660867669957927e72]
    extrapolation_scalars = T[-1.0, 0.25, -0.015625, 0.000244140625, -9.5367431640625e-7, 9.313225746154785e-10, -2.2737367544323206e-13, 1.3877787807814457e-17, -2.117582368135751e-22, 8.077935669463161e-28, -7.703719777548943e-34, 1.8367099231598242e-40, -1.0947644252537633e-47, 1.6313261169996311e-55, -6.077163357286271e-64, 5.659799424266695e-73]
    extrapolation_scalars_2 = T[-0.25, 0.015625, -0.000244140625, 9.5367431640625e-7, -9.313225746154785e-10, 2.2737367544323206e-13, -1.3877787807814457e-17, 2.117582368135751e-22, -8.077935669463161e-28, 7.703719777548943e-34, -1.8367099231598242e-40, 1.0947644252537633e-47, -1.6313261169996311e-55, 6.077163357286271e-64, -5.659799424266695e-73]
  else # sequence == :bulirsch
    subdividing_sequence = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    extrapolation_weights = T[-1.0 -1.3333333333333333 -1.5 -1.6 -1.6457142857142857 -1.6718367346938776 -1.6835279006707577 -1.6901299708694666 -1.693069327340544 -1.6947243315705933 -1.6954602083971546 -1.695874240194077 -1.6960582742950205 -1.6961617997954963 -1.6962078123772122 -1.6962336948493626; 0.0 5.333333333333333 38.4 204.8 921.6 3932.16 16178.029714285714 65739.29534693877 264796.0427960611 1.063337834600653e6 4.260748471165052e6 1.7059653702729277e7 6.82682451256418e7 2.7313966499109036e8 1.0926772230310967e9 4.3709756753076935e9; 0.0 0.0 -72.9 -1499.6571428571428 -17995.885714285716 -188466.0031168831 -1.809273629922078e6 -1.687678722000189e7 -1.5430205458287445e8 -1.4010322512667694e9 -1.2658738458504457e10 -1.1417952888042778e11 -1.0286202719081353e12 -9.262670584090748e12 -8.338439277458397e13 -7.505626090600242e14; 0.0 0.0 0.0 1872.4571428571428 53926.76571428571 1.1504376685714286e6 2.0707878034285713e7 3.5341445178514284e8 5.816192120806923e9 9.4536202090576e10 1.5231567106062036e12 2.4466077986835336e13 3.9213804300291206e14 6.280341834369219e15 1.0052910177255184e17 1.608858416060063e18; 0.0 0.0 0.0 0.0 -57586.834285714285 -4.738573792653061e6 -2.2745154204734695e8 -9.528151870492496e9 -3.6588103182691187e11 -1.36516582563434e13 -4.992606448034158e14 -1.8132753113333124e16 -6.553390301665806e17 -2.3644157580681015e19 -8.52021725370699e20 -3.0689640436338757e22; 0.0 0.0 0.0 0.0 0.0 5.09977563903525e6 5.874941536168609e8 5.013283444197213e10 3.609564079821993e12 2.464129078491814e14 1.6221009705271826e16 1.0546231071871968e18 6.7967878012847596e19 4.36700279749998e21 2.7997424543832918e23 1.7935867203368857e25; 0.0 0.0 0.0 0.0 0.0 0.0 -5.70060368320064e8 -1.8763129837277533e11 -3.602520928757287e13 -6.036515068986755e15 -9.272087145963656e17 -1.3838308524243085e20 -2.0243468469749884e22 -2.9409072775127474e24 -4.251513956009017e26 -6.1356617753586065e28; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.9561237425512534e11 9.013818205676177e13 3.0767166142041348e16 8.860943848907908e18 2.419628400341786e21 6.371227239299971e23 1.6569236045823926e26 4.271386836316456e28 1.0977631674699422e31; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -8.554159840110055e13 -1.1262162440922038e17 -8.649340754628125e19 -5.797259955974749e22 -3.561836516950886e25 -2.126373139447408e28 -1.2442320541680835e31 -7.230324405099858e33; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1651370210969312e17 2.1475805572858636e20 2.9321633208809658e23 3.377852145654872e26 3.689515303627295e29 3.8860083472261895e32 4.0424356038700874e35; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0269737017532947e20 -1.0674622648776207e24 -3.279244077704051e27 -8.791712994944155e30 -2.1606513856374756e34 -5.159530537984771e37; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1022516982215851e24 8.126681320648103e27 4.438251558583284e31 2.0451463181951773e35 8.935380607282608e38; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -7.659880204931148e27 -1.6135647078547534e32 -1.9827483130119208e36 -2.126313710861715e40; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.665360734159382e32 4.911348648324117e36 1.0729004833885644e41; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.62766803500927e36 -3.8992995301161535e41; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.023990816545532e41]
    extrapolation_weights_2 = T[-4.0 -28.8 -153.6 -691.2 -2949.12 -12133.522285714285 -49304.47151020408 -198597.0320970458 -797503.3759504899 -3.1955613533737888e6 -1.279474027704696e7 -5.120118384423134e7 -2.0485474874331778e8 -8.195079172733225e8 -3.27823175648077e9; 0.0 64.8 1333.0285714285715 15996.342857142858 167525.33610389612 1.6082432265974027e6 1.5001588640001683e7 1.3715738185144395e8 1.2453620011260173e9 1.1252211963115074e10 1.0149291456038025e11 9.143291305850092e11 8.233484963636221e12 7.411946024407464e13 6.671667636089105e14; 0.0 0.0 -1755.4285714285713 -50556.34285714286 -1.0785353142857142e6 -1.941363565714286e7 -3.313260485485714e8 -5.45268011325649e9 -8.862768945991501e10 -1.427959416193316e12 -2.2936948112658125e13 -3.6762941531523006e14 -5.887820469721143e15 -9.424603291176736e16 -1.508304765056309e18; 0.0 0.0 0.0 55987.2 4.606946742857143e6 2.2113344365714285e8 9.263480985201038e9 3.557176698317199e11 1.327244552700053e13 4.853922935588765e14 1.7629065526851648e16 6.37135168217509e17 2.29873754256621e19 8.283544552215128e20 2.9837150424218233e22; 0.0 0.0 0.0 0.0 -5.020091644675325e6 -5.783145574665974e8 -4.9349508903816315e10 -3.5531646410747744e12 -2.4256270616403794e14 -1.5967556428626954e16 -1.0381446211373969e18 -6.690587991889685e19 -4.2987683787890434e21 -2.755996478533553e23 -1.765561927831622e25; 0.0 0.0 0.0 0.0 0.0 5.661016157622857e8 1.863283032451866e11 3.577503422307583e13 5.994594825452124e15 9.207697651894463e17 1.3742209159491396e20 2.0102888827598844e22 2.920484310307798e24 4.221989553536732e26 6.093053013029727e28; 0.0 0.0 0.0 0.0 0.0 0.0 -1.9484826341819125e11 -8.978607978310253e13 -3.0646981899298996e16 -8.826330786998111e18 -2.410176726902951e21 -6.346339632896457e23 -1.6504512467519928e26 -4.254701731487095e28 -1.0934750300970127e31; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 8.539308868165419e13 1.1242610075573214e17 8.63432453804023e19 5.787195268551182e22 3.555652772997846e25 2.1226815194136454e28 1.2420719290740417e31 7.217771758563227e33; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.1639991919747662e17 -2.145483310647889e20 -2.929299880137918e23 -3.3745534619188814e26 -3.685912261338597e29 -3.882213417199601e32 -4.0384879128506835e35; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0260939388619086e20 1.0669989566029343e24 3.277820794684214e27 8.787897147290099e30 2.1597136029180146e34 5.157291158410993e37; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.1019825938030739e24 -8.124697267591304e27 -4.437168001073864e31 -2.0446470148948364e35 -8.933199117876534e38; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 7.6590490547353e27 1.6133896248786405e32 1.9825331710508737e36 2.126082991058019e40; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.66525908860676e32 -4.911048883391968e36 -1.07283499873992e41; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.627542501479674e36 3.8991937548467816e41; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.023929415318473e41]
    extrapolation_scalars = T[-1.0, 0.25, -0.027777777777777776, 0.001736111111111111, -4.8225308641975306e-5, 7.535204475308642e-7, -5.232780885631001e-9, 2.0440550334496098e-11, -3.548706655294462e-14, 3.465533843060998e-17, -1.5041379527174468e-20, 3.672211798626579e-24, -3.9846048162180767e-28, 2.4320097755237285e-32, -6.597248740027475e-37, 1.0066602691692314e-41]
    extrapolation_scalars_2 = T[-0.25, 0.027777777777777776, -0.001736111111111111, 4.8225308641975306e-5, -7.535204475308642e-7, 5.232780885631001e-9, -2.0440550334496098e-11, 3.548706655294462e-14, -3.465533843060998e-17, 1.5041379527174468e-20, -3.672211798626579e-24, 3.9846048162180767e-28, -2.4320097755237285e-32, 6.597248740027475e-37, -1.0066602691692314e-41]
  end
  extrapolation_coefficients(subdividing_sequence,
      extrapolation_weights, extrapolation_scalars,
      extrapolation_weights_2, extrapolation_scalars_2)
end

@cache mutable struct ExtrapolationMidpointDeuflhardConstantCache{QType, extrapolation_coefficients} <: OrdinaryDiffEqConstantCache
  # Values that are mutated
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.n_min - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

  # Constant values
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n + alg.n_min - 1)
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType),alg.n_max - alg.n_min + 1)
    n_curr = alg.n_init
    n_old = alg.n_init

    coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits),alg)
    stage_number = Vector{Int}(undef, alg.n_max - alg.n_min + 1)
    for n in 1:length(stage_number)
      s = zero(eltype(coefficients.subdividing_sequence))
      for i in alg.n_min:(alg.n_min + n)
        s += coefficients.subdividing_sequence[i]
      end
      stage_number[n] = 2 * Int(s) - alg.n_min - n + 1
    end

    # Initialize cache
    ExtrapolationMidpointDeuflhardConstantCache(Q, n_curr, n_old, coefficients, stage_number)
end



@cache mutable struct ExtrapolationMidpointDeuflhardCache{uType,uNoUnitsType,rateType,QType,extrapolation_coefficients} <: OrdinaryDiffEqMutableCache
  # Values that are mutated
  utilde::uType
  u_temp1::uType
  u_temp2::uType
  u_temp3::Array{uType,1}
  u_temp4::Array{uType,1}
  tmp::uType # for get_tmp_cache()
  T::Array{uType,1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
  res::uNoUnitsType # Storage for the scaled residual of u and utilde

  fsalfirst::rateType
  k::rateType
  k_tmps::Array{rateType,1}

  # Constant values
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.n_min - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # Stage_number[n] contains information for extrapolation order (n + alg.n_min - 1)
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  # Initialize cache's members
  utilde = zero(u)
  u_temp1 = zero(u)
  u_temp2 = zero(u)
  u_temp3 = Array{typeof(u),1}(undef, Threads.nthreads())
  u_temp4 = Array{typeof(u),1}(undef, Threads.nthreads())

  for i=1:Threads.nthreads()
      u_temp3[i] = zero(u)
      u_temp4[i] = zero(u)
  end

  tmp = zero(u)
  T = Vector{typeof(u)}(undef,alg.n_max + 1)
  for i in 1:alg.n_max+1
    T[i] = zero(u)
  end
  res = uEltypeNoUnits.(zero(u))

  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  for i=1:Threads.nthreads()
      k_tmps[i] = zero(rate_prototype)
  end

  cc =  alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val(false))
  # Initialize cache
  ExtrapolationMidpointDeuflhardCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp, T, res, fsalfirst, k, k_tmps, cc.Q, cc.n_curr, cc.n_old, cc.coefficients,cc.stage_number)
end

@cache mutable struct ImplicitDeuflhardExtrapolationConstantCache{QType,extrapolation_coefficients,TF,UF} <: OrdinaryDiffEqConstantCache
  # Values that are mutated
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.n_min - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

  # Constant values
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n + alg.n_min - 1)

  tf::TF
  uf::UF
end

@cache mutable struct ImplicitDeuflhardExtrapolationCache{uType,QType,extrapolation_coefficients,rateType,JType,WType,F,JCType,GCType,uNoUnitsType,TFType,UFType} <: OrdinaryDiffEqMutableCache
  # Values that are mutated
  utilde::uType
  u_temp1::uType
  u_temp2::uType
  u_temp3::Array{uType,1}
  u_temp4::Array{uType,1}
  tmp::uType # for get_tmp_cache()
  T::Array{uType,1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
  res::uNoUnitsType # Storage for the scaled residual of u and utilde

  fsalfirst::rateType
  k::rateType
  k_tmps::Array{rateType,1}

  # Constant values
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.n_min - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # Stage_number[n] contains information for extrapolation order (n + alg.n_min - 1)


  du1::rateType
  du2::rateType
  J::JType
  W::WType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

function alg_cache(alg::ImplicitDeuflhardExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  # Initialize cache's members
  QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

  Q = fill(zero(QType),alg.n_max - alg.n_min + 1)
  n_curr = alg.n_init
  n_old = alg.n_init

  coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits),alg)
  stage_number = Vector{Int}(undef, alg.n_max - alg.n_min + 1)
  for n in 1:length(stage_number)
    s = zero(eltype(coefficients.subdividing_sequence))
    for i in alg.n_min:(alg.n_min + n)
      s += coefficients.subdividing_sequence[i]
    end
    stage_number[n] = 2 * Int(s) - alg.n_min - n + 1
  end

  tf = TimeDerivativeWrapper(f,u,p)
  uf = UDerivativeWrapper(f,t,p)
  ImplicitDeuflhardExtrapolationConstantCache(Q,n_curr,n_old,coefficients,stage_number,tf,uf)
end

function alg_cache(alg::ImplicitDeuflhardExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  utilde = zero(u)
  u_temp1 = zero(u)
  u_temp2 = zero(u)
  u_temp3 = Array{typeof(u),1}(undef, Threads.nthreads())
  u_temp4 = Array{typeof(u),1}(undef, Threads.nthreads())

  for i=1:Threads.nthreads()
      u_temp3[i] = zero(u)
      u_temp4[i] = zero(u)
  end

  tmp = zero(u)
  T = Vector{typeof(u)}(undef,alg.n_max + 1)
  for i in 1:alg.n_max+1
    T[i] = zero(u)
  end
  res = uEltypeNoUnits.(zero(u))

  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  for i=1:Threads.nthreads()
      k_tmps[i] = zero(rate_prototype)
  end

  cc =  alg_cache(alg::ImplicitDeuflhardExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val(false))

  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)

  if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
    W = WOperator(f, dt, true)
    J = nothing # is J = W.J better?
  else
    J = false .* vec(rate_prototype) .* vec(rate_prototype)' # uEltype?
    W = similar(J)
  end
  tf = TimeGradientWrapper(f,uprev,p)
  uf = UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,du1,du2)


  ImplicitDeuflhardExtrapolationCache(utilde,u_temp1,u_temp2,u_temp3,u_temp4,tmp,T,res,fsalfirst,k,k_tmps,cc.Q,cc.n_curr,cc.n_old,cc.coefficients,cc.stage_number,
    du1,du2,J,W,tf,uf,linsolve_tmp,linsolve,jac_config,grad_config)
end

@cache mutable struct ExtrapolationMidpointHairerWannerConstantCache{QType,extrapolation_coefficients} <: OrdinaryDiffEqConstantCache
  # Values that are mutated
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

  # Constant values
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
  sigma::Rational{Int} # Parameter for order selection
end

function alg_cache(alg::ExtrapolationMidpointHairerWanner,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  # Initialize cache's members
  QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

  Q = fill(zero(QType),alg.n_max + 1)
  n_curr = alg.n_init
  n_old = alg.n_init

  coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits),alg)
  stage_number = Vector{Int}(undef, alg.n_max + 1)
  for n in 1:length(stage_number)
    s = zero(eltype(coefficients.subdividing_sequence))
    for i in 1:n
      s += coefficients.subdividing_sequence[i]
    end
    stage_number[n] = 2 * Int(s) - n + 1
  end
  sigma = 9//10

  # Initialize the constant cache
  ExtrapolationMidpointHairerWannerConstantCache(Q, n_curr, n_old, coefficients, stage_number, sigma)
end

@cache mutable struct ExtrapolationMidpointHairerWannerCache{uType,uNoUnitsType,rateType,QType,extrapolation_coefficients} <: OrdinaryDiffEqMutableCache
  # Values that are mutated
  utilde::uType
  u_temp1::uType
  u_temp2::uType
  u_temp3::Array{uType,1}
  u_temp4::Array{uType,1}
  tmp::uType # for get_tmp_cache()
  T::Array{uType,1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
  res::uNoUnitsType # Storage for the scaled residual of u and utilde

  fsalfirst::rateType
  k::rateType
  k_tmps::Array{rateType,1}

  # Constant values
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
  sigma::Rational{Int} # Parameter for order selection
end


function alg_cache(alg::ExtrapolationMidpointHairerWanner,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  # Initialize cache's members
  utilde = zero(u)
  u_temp1 = zero(u)
  u_temp2 = zero(u)
  u_temp3 = Array{typeof(u),1}(undef, Threads.nthreads())
  u_temp4 = Array{typeof(u),1}(undef, Threads.nthreads())

  for i=1:Threads.nthreads()
    u_temp3[i] = zero(u)
    u_temp4[i] = zero(u)
  end
  tmp = zero(u)
  T = Vector{typeof(u)}(undef,alg.n_max + 1)
  for i in 1:alg.n_max+1
    T[i] = zero(u)
  end
  res = uEltypeNoUnits.(zero(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  for i=1:Threads.nthreads()
      k_tmps[i] = zero(rate_prototype)
  end

  cc = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val(false))

  # Initialize the cache
  ExtrapolationMidpointHairerWannerCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp, T, res, fsalfirst, k, k_tmps,
      cc.Q, cc.n_curr, cc.n_old, cc.coefficients, cc.stage_number, cc.sigma)
end

@cache mutable struct ImplicitHairerWannerExtrapolationConstantCache{QType,extrapolation_coefficients,TF,UF} <: OrdinaryDiffEqConstantCache
  # Values that are mutated
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

  # Constant values
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
  sigma::Rational{Int} # Parameter for order selection

  tf::TF
  uf::UF
end

function alg_cache(alg::ImplicitHairerWannerExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  # Initialize cache's members
  QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

  Q = fill(zero(QType),alg.n_max + 1)
  n_curr = alg.n_init
  n_old = alg.n_init

  coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits),alg)
  stage_number = Vector{Int}(undef, alg.n_max + 1)
  for n in 1:length(stage_number)
    s = zero(eltype(coefficients.subdividing_sequence))
    for i in 1:n
      s += coefficients.subdividing_sequence[i]
    end
    stage_number[n] = 2 * Int(s) - n + 1
  end
  sigma = 9//10

  # Initialize the constant cache
  tf = TimeDerivativeWrapper(f,u,p)
  uf = UDerivativeWrapper(f,t,p)
  ImplicitHairerWannerExtrapolationConstantCache(Q, n_curr, n_old, coefficients, stage_number, sigma, tf, uf)
end

@cache mutable struct ImplicitHairerWannerExtrapolationCache{uType,uNoUnitsType,rateType,QType,extrapolation_coefficients,JType,WType,F,JCType,GCType,TFType,UFType} <: OrdinaryDiffEqMutableCache
  # Values that are mutated
  utilde::uType
  u_temp1::uType
  u_temp2::uType
  u_temp3::Array{uType,1}
  u_temp4::Array{uType,1}
  tmp::uType # for get_tmp_cache()
  T::Array{uType,1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
  res::uNoUnitsType # Storage for the scaled residual of u and utilde

  fsalfirst::rateType
  k::rateType
  k_tmps::Array{rateType,1}

  # Constant values
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
  n_curr::Int # Storage for the current extrapolation order
  n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
  sigma::Rational{Int} # Parameter for order selection

  du1::rateType
  du2::rateType
  J::JType
  W::WType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end


function alg_cache(alg::ImplicitHairerWannerExtrapolation,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  # Initialize cache's members
  utilde = zero(u)
  u_temp1 = zero(u)
  u_temp2 = zero(u)
  u_temp3 = Array{typeof(u),1}(undef, Threads.nthreads())
  u_temp4 = Array{typeof(u),1}(undef, Threads.nthreads())

  for i=1:Threads.nthreads()
    u_temp3[i] = zero(u)
    u_temp4[i] = zero(u)
  end
  tmp = zero(u)
  T = Vector{typeof(u)}(undef,alg.n_max + 1)
  for i in 1:alg.n_max+1
    T[i] = zero(u)
  end
  res = uEltypeNoUnits.(zero(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  for i=1:Threads.nthreads()
      k_tmps[i] = zero(rate_prototype)
  end

  cc = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val(false))

  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)

  if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
    W = WOperator(f, dt, true)
    J = nothing # is J = W.J better?
  else
    J = false .* vec(rate_prototype) .* vec(rate_prototype)' # uEltype?
    W = similar(J)
  end
  tf = TimeGradientWrapper(f,uprev,p)
  uf = UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,du1,du2)

  # Initialize the cache
  ImplicitHairerWannerExtrapolationCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp, T, res, fsalfirst, k, k_tmps,
      cc.Q, cc.n_curr, cc.n_old, cc.coefficients, cc.stage_number, cc.sigma, du1, du2, J, W, tf, uf, linsolve_tmp,
      linsolve, jac_config, grad_config)
end
