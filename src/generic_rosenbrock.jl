abstract type RosenbrockTableau end
struct RosenbrockFixedTableau{T,T2}<:RosenbrockTableau
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

struct RosenbrockAdaptiveTableau{T,T2}<:RosenbrockTableau
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    btilde::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

function _common_nonzero_vals(tab::RosenbrockTableau)
    nzvals=[]
    push!(nzvals,[Symbol(:a,ind[1],ind[2]) for ind in findall(!iszero,tab.a)])
    push!(nzvals,[Symbol(:C,ind[1],ind[2]) for ind in findall(!iszero,tab.C)])
    push!(nzvals,[Symbol(:b,ind) for ind in findall(!iszero,tab.b)])
    push!(nzvals,:gamma)
    push!(nzvals,[Symbol(:d,ind) for ind in findall(!iszero,tab.d)])
    push!(nzvals,[Symbol(:c,ind) for ind in findall(!iszero,tab.c)])
    nzvals
end

function _nonzero_vals(tab::RosenbrockFixedTableau)
    nzvals=_common_nonzero_vals(tab)
    vcat(nzvals...)
end

function _nonzero_vals(tab::RosenbrockAdaptiveTableau)
    nzvals=_common_nonzero_vals(tab)
    push!(nzvals,[Symbol(:btilde,ind) for ind in findall(!iszero,tab.btilde)])
    vcat(nzvals...)
end

function _push_assigns!(assignexprs,inds,name,T,arr)
    for ind in inds
        if length(size(arr))==2
            indstr="$(ind[1])$(ind[2])"
        elseif length(size(arr))==1
            indstr="$(ind[1])"
        end
        push!(assignexprs,:($(Symbol(name,indstr))=convert($T,$(arr[ind]))))
    end
end

function gen_tableau(tab::RosenbrockTableau,tabname::Symbol,tabstructname::Symbol,skipstructdef::Bool=false)
    ainds=findall(!iszero,tab.a)
    Cinds=findall(!iszero,tab.C)
    assignexprs=Array{Expr,1}()
    _push_assigns!(assignexprs,ainds,:a,:T,tab.a)
    _push_assigns!(assignexprs,Cinds,:C,:T,tab.C)
    _push_assigns!(assignexprs,eachindex(tab.b),:b,:T,tab.b)
    push!(assignexprs,:(gamma=convert(T2,$(tab.gamma))))
    _push_assigns!(assignexprs,eachindex(tab.d),:d,:T,tab.d)
    _push_assigns!(assignexprs,findall(!iszero,tab.c),:c,:T2,tab.c)
    structexprs=[:($(assignexpr.args[1])::$(assignexpr.args[2].args[2])) for assignexpr in assignexprs]
    valsymbols=[assignexpr.args[1] for assignexpr in assignexprs]
    tabinitexpr=:($tabstructname($(valsymbols...)))
    tabstructexpr=quote
        struct $tabstructname{T,T2} <: OrdinaryDiffEqConstantCache
            $(structexprs...)
        end
    end
    if skipstructdef
        tabstructexpr=:nothing
    end
    tabexpr=quote
        function $tabname(T::Type,T2::Type)
            $(assignexprs...)
            $tabinitexpr
        end
    end
    quote
        $tabstructexpr
        $tabexpr
    end
end

function gen_cache(tab::RosenbrockTableau,algname::Symbol,cachename::Symbol,constcachename::Symbol,tabname::Symbol)
    kstype=[:($(Symbol(:k,i))::rateType) for i in 1:length(tab.b)]
    ksinit=[:($(Symbol(:k,i))=zero(rate_prototype)) for i in 1:length(tab.b)]
    ks=[Symbol(:k,i) for i in 1:length(tab.b)]
    quote    
        struct $constcachename{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
            tf::TF
            uf::UF
            tab::Tab
        end  
        function alg_cache(alg::$algname,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
            tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
            uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
            $constcachename(tf,uf,$tabname(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)))
        end
        @cache mutable struct $cachename{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
            u::uType
            uprev::uType
            du::rateType
            du1::rateType
            du2::rateType
            $(kstype...)
            fsalfirst::rateType
            fsallast::rateType
            dT::rateType
            J::JType
            W::WType
            tmp::rateType
            tab::TabType
            tf::TFType
            uf::UFType
            linsolve_tmp::rateType
            linsolve::F
            jac_config::JCType
            grad_config::GCType
        end
        function alg_cache(alg::$algname,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
            du = zero(rate_prototype)
            du1 = zero(rate_prototype)
            du2 = zero(rate_prototype)
            $(ksinit...)
            fsalfirst = zero(rate_prototype)
            fsallast = zero(rate_prototype)
            dT = zero(rate_prototype)
            if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype !== nothing
              W = WOperator(f, dt, true)
              J = nothing # is J = W.J better?
            else
              J = false .* vec(rate_prototype) .* vec(rate_prototype)' # uEltype?
              W = similar(J)
            end
            tmp = zero(rate_prototype)
            tab = $tabname(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
          
            tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
            uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
            linsolve_tmp = zero(rate_prototype)
            linsolve = alg.linsolve(Val{:init},uf,u)
            grad_config = build_grad_config(alg,f,tf,du1,t)
            jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
            $cachename(u,uprev,du,du1,du2,$(ks...),
                              fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                              linsolve,jac_config,grad_config)
        end
    end
end

function gen_initialize(cachename::Symbol,constcachename::Symbol)
    quote
        function initialize!(integrator, cache::$constcachename)
            integrator.kshortsize = 2
            integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
            integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
            integrator.destats.nf += 1
          
            # Avoid undefined entries if k is an array of arrays
            integrator.fsallast = zero(integrator.fsalfirst)
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
          end
          
          function initialize!(integrator, cache::$cachename)
            integrator.kshortsize = 2
            @unpack fsalfirst,fsallast = cache
            integrator.fsalfirst = fsalfirst
            integrator.fsallast = fsallast
            resize!(integrator.k, integrator.kshortsize)
            integrator.k .= [fsalfirst,fsallast]
            integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
            integrator.destats.nf += 1
          end
    end
end

function gen_constant_perform_step(tab::RosenbrockTableau,cachename::Symbol)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=_nonzero_vals(tab)
    dtCijexpr=quote end
    for Cind in findall(!iszero,tab.C)
        push!(dtCijexpr.args,:($(Symbol(:dtC,"$(Cind[1])$(Cind[2])"))=$(Symbol(:C,"$(Cind[1])$(Cind[2])"))/dt))
    end
    dtdiexpr=quote end
    for dind in eachindex(tab.d)
        push!(dtdiexpr.args,:($(Symbol(:dtd,dind[1]))=dt*$(Symbol(:d,dind[1]))))
    end
    iterexpr=quote linsolve_tmp = integrator.fsalfirst + dtd1*dT end
    n=length(tab.d)
    for i in 2:n
        uexpr=:(u=(+ uprev))
        for aind in findall(!iszero,tab.a[i,:])
            push!(uexpr.args[2].args,:($(Symbol(:a,i,aind))*$(Symbol(:k,aind))))
        end
        linsolveexpr=:(linsolve_tmp=du+$(Symbol(:dtd,i))*dT)
        for Cind in findall(!iszero,tab.C[i,:])
            push!(linsolveexpr.args[2].args,:($(Symbol(:dtC,i,Cind))*$(Symbol(:k,Cind))))
        end
        push!(iterexpr.args,
        quote
            $(Symbol(:k,i-1)) = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
            integrator.destats.nsolve += 1
            $uexpr
            du = f(u, p, t+$(Symbol(:c,i))*dt)
            integrator.destats.nf += 1
            $linsolveexpr
        end)
    end
    uexpr=:(u=(+ uprev))
    for i in 1:n
        push!(uexpr.args[2].args,:($(Symbol(:b,i))*$(Symbol(:k,i))))
    end
    push!(iterexpr.args,
    quote
        $(Symbol(:k,n))=_reshape(W\-_vec(linsolve_tmp), axes(uprev))
        integrator.destats.nsolve += 1
        $uexpr
    end)
    quote
        @muladd function perform_step!(integrator, cache::$cachename, repeat_step=false)
            @unpack t,dt,uprev,u,f,p = integrator
            @unpack tf,uf = cache
            $unpacktabexpr

            $dtCijexpr
            $dtdiexpr
            dtgamma = dt*gamma

            # Time derivative
            tf.u = uprev
            dT = ForwardDiff.derivative(tf, t)

            W = calc_W!(integrator, cache, dtgamma, repeat_step, true)

            $iterexpr

            integrator.fsallast = f(u, p, t + dt)
            integrator.destats.nf += 1
          
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
            integrator.u = u
        end
    end

end

function gen_perform_step(tab::RosenbrockTableau,cachename::Symbol)
    n=length(tab.b)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=_nonzero_vals(tab)
    ks=[Symbol(:k,i) for i in 1:n]
    dtCij=[:($(Symbol(:dtC,"$(Cind[1])$(Cind[2])"))=$(Symbol(:C,"$(Cind[1])$(Cind[2])"))/dt) for Cind in findall(!iszero,tab.C)]
    dtdi=[:($(Symbol(:dtd,dind[1]))=dt*$(Symbol(:d,dind[1]))) for dind in eachindex(tab.d)]
    iterexprs=[]
    for i in 1:(n-1)
        ki=Symbol(:k,i)
        dtdj=Symbol(:dtd,i+1)
        aijkj=[:($(Symbol(:a,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tab.a[i+1,:])]
        dtCijkj=[:($(Symbol(:dtC,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tab.C[i+1,:])]
        repeatstepexpr=[]
        if i==1
            repeatstepexpr=[:(!repeat_step)]
        end
        push!(iterexprs,quote
            if DiffEqBase.has_invW(f)
                mul!(vec($ki), W, vec(linsolve_tmp))
            else
                cache.linsolve(vec($ki), W, vec(linsolve_tmp), $(repeatstepexpr...))
                @.. $ki = -$ki
            end
            integrator.destats.nsolve += 1
            @.. u = +(uprev,$(aijkj...))
            f( du,  u, p, t+$(Symbol(:c,i+1))*dt)
            integrator.destats.nf += 1
            if mass_matrix == I
                @.. linsolve_tmp = +(du,$dtdj*dT,$(dtCijkj...))
            else
                @.. du1 = +($(dtCijkj...))
                mul!(du2,mass_matrix,du1)
                @.. linsolve_tmp = du + $dtdj*dT + du2
            end
        end)
    end
    klast=Symbol(:k,n)
    biki=[:($(Symbol(:b,i))*$(Symbol(:k,i))) for i in 1:n]
    push!(iterexprs,quote
        if DiffEqBase.has_invW(f)
            mul!(vec($klast), W, vec(linsolve_tmp))
        else
            cache.linsolve(vec($klast), W, vec(linsolve_tmp))
            @.. $klast = -$klast
        end
        integrator.destats.nsolve += 1
        @.. u = +(uprev,$(biki...))
        f( fsallast,  u, p, t + dt)
        integrator.destats.nf += 1
    end)
    quote
        @muladd function perform_step!(integrator, cache::$cachename, repeat_step=false)
            @unpack t,dt,uprev,u,f,p = integrator
            @unpack du,du1,du2,fsallast,dT,J,W,uf,tf,$(ks...),linsolve_tmp,jac_config = cache
            $unpacktabexpr

            # Assignments
            sizeu  = size(u)
            uidx = eachindex(integrator.uprev)
            mass_matrix = integrator.f.mass_matrix
            atmp = du # does not work with units - additional unitless array required!
          
            # Precalculations
            $(dtCij...)
            $(dtdi...)
            dtgamma = dt*gamma
            calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

            $(iterexprs...)
        end
    end
end

function RosenbrockW6S4OSTableau()
    a=[0                  0                  0                  0                  0;
       0.5812383407115008 0                  0                  0                  0;
       0.9039624413714670 1.8615191555345010 0                  0                  0;
       2.0765797196750000 0.1884255381414796 1.8701589674910320 0                  0;
       4.4355506384843120 5.4571817986101890 4.6163507880689300 3.1181119524023610 0;
       10.791701698483260 -10.05691522584131 14.995644854284190 5.2743399543909430 1.4297308712611900]
    C=[0                  0                  0                  0                  0;
       -2.661294105131369 0                  0                  0                  0;
       -3.128450202373838 0.0000000000000000 0                  0                  0;
       -6.920335474535658 -1.202675288266817 -9.733561811413620 0                  0;
       -28.09530629102695 20.371262954793770 -41.04375275302869 -19.66373175620895 0;
       9.7998186780974000 11.935792886603180 3.6738749290132010 14.807828541095500 0.8318583998690680]
    b=[6.4562170746532350,-4.853141317768053,9.7653183340692600,2.0810841772787230,0.6603936866352417,0.6000000000000000]
    gamma=0.2500000000000000
    d=[0.2500000000000000,0.0836691184292894,0.0544718623516351,-0.3402289722355864,0.0337651588339529,-0.0903074267618540]
    c=[0                 ,0.1453095851778752,0.3817422770256738,0.6367813704374599,0.7560744496323561,0.9271047239875670]
    RosenbrockFixedTableau(a,C,b,gamma,d,c)
end

macro RosenbrockW6S4OS(part)
    tab=RosenbrockW6S4OSTableau()
    algname=:RosenbrockW6S4OS
    tabname=:RosenbrockW6S4OSConstantCache
    tabstructname=:RosenbrockW6S4OSConstantCache
    cachename=:RosenbrockW6SCache
    constcachename=:RosenbrockW6SConstantCache
    if part.value==:constcache
        println("Generating const cache")
        return esc(gen_tableau(tab,tabname,tabstructname))
    elseif part.value==:cache
        println("Generating cache")
        return esc(gen_cache(tab,algname,cachename,constcachename,tabname))
    elseif part.value==:init
        println("Generating initialize")
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        println("Generating perform_step")
        constperformstepexpr=gen_constant_perform_step(tab,constcachename)
        performstepexpr=gen_perform_step(tab,cachename)
        return esc(quote $([constperformstepexpr,performstepexpr]...) end)
    else
        println("Unknown parameter!")
    end
end

