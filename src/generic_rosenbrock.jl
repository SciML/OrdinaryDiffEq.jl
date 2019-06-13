abstract type RosenbrockTableau{T,T2} end
struct RosenbrockFixedTableau{T,T2}<:RosenbrockTableau{T,T2}
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

struct RosenbrockAdaptiveTableau{T,T2}<:RosenbrockTableau{T,T2}
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    btilde::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

macro _bitarray2boolarray(expr)
    args=[:([i for i in $arg]) for arg in expr.args[2:end]]
    args[end-2]=:(tab.gamma!=0)
    esc(:($(expr.args[1])($(args...))))
end
_masktab(tab::RosenbrockFixedTableau)=@_bitarray2boolarray RosenbrockFixedTableau(tab.a.!=0,tab.C.!=0,tab.b.!=0,tab.gamma!=0,tab.d.!=0,tab.c.!=0)
_masktab(tab::RosenbrockAdaptiveTableau)=@_bitarray2boolarray RosenbrockAdaptiveTableau(tab.a.!=0,tab.C.!=0,tab.b.!=0,tab.btilde.!=0,tab.gamma!=0,tab.d.!=0,tab.c.!=0)

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

function _push_assigns!(assignexprs,inds,name,T::Symbol)
    for ind in inds
        push!(assignexprs,:($(Symbol(name,"$(Tuple(ind)...)"))::$T))
    end
end

function gen_tableau_struct(tab::RosenbrockTableau,tabstructname::Symbol)
    valexprs=Array{Expr,1}()
    _push_assigns!(valexprs,findall(!iszero,tab.a),:a,:T)
    _push_assigns!(valexprs,findall(!iszero,tab.C),:C,:T)
    _push_assigns!(valexprs,eachindex(tab.b),:b,:T)
    if typeof(tab)<:RosenbrockAdaptiveTableau
        _push_assigns!(valexprs,eachindex(tab.btilde),:btilde,:T)
    end
    push!(valexprs,:(gamma::T2))
    _push_assigns!(valexprs,eachindex(tab.d),:d,:T)
    _push_assigns!(valexprs,findall(!iszero,tab.c),:c,:T2)
    quote struct $tabstructname{T,T2} <: OrdinaryDiffEqConstantCache
        $(valexprs...)
        end
    end
end

function gen_tableau(tab::RosenbrockTableau,tabstructexpr::Expr,tabname::Symbol)
    tabstructname=tabstructexpr.args[2].args[2].args[1].args[1]
    valexprs=[valexpr for valexpr in tabstructexpr.args[2].args[3].args if typeof(valexpr)!=LineNumberNode]
    valsym2tabdict=Dict("a"=>tab.a,"C"=>tab.C,"gamma"=>tab.gamma,"c"=>tab.c,"d"=>tab.d,"b"=>tab.b)
    if typeof(tab)<:RosenbrockAdaptiveTableau
        valsym2tabdict["btilde"]=tab.btilde
    end
    pattern=r":([a-zA-Z]+)([1-9]?[1-9]?)"
    assignexprs=Array{Expr,1}()
    for valexpr in valexprs
        valsym=valexpr.args[1]
        valtype=valexpr.args[2]
        m=match(pattern,repr(valsym))
        val=valsym2tabdict[m[1]][CartesianIndex([parse(Int,i) for i in m[2]]...)]
        push!(assignexprs,:($valsym=convert($valtype,$val)))
    end
    quote function $tabname(T::Type,T2::Type)
            $(assignexprs...)
            $tabstructname($([valexpr.args[1] for valexpr in valexprs]...))
        end
    end
end

function gen_cache_struct(tab::RosenbrockTableau,cachename::Symbol,constcachename::Symbol)
    kstype=[:($(Symbol(:k,i))::rateType) for i in 1:length(tab.b)]
    constcacheexpr=quote struct $constcachename{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
            tf::TF
            uf::UF
            tab::Tab
        end
    end
    cacheexpr=quote
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
    end
    constcacheexpr,cacheexpr
end

function gen_algcache(cacheexpr::Expr,cachename::Symbol,constcachename::Symbol,algname::Symbol,tabname::Symbol)
    #cachename=cacheexpr.args[2].args[3].args[2].args[1].args[1]
    valsyms=[valexpr.args[1] for valexpr in cacheexpr.args[2].args[3].args[3].args if typeof(valexpr)!=LineNumberNode]
    ksinit=[:($valsym=zero(rate_prototype)) for valsym in valsyms if match(r":k[1-9]",repr(valsym))!=nothing]
    quote
        function alg_cache(alg::$algname,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
            tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
            uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
            $constcachename(tf,uf,$tabname(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)))
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
            $cachename($(valsyms...))
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

function gen_constant_perform_step(tabmask::RosenbrockTableau{Bool,Bool},cachename::Symbol,n_normalstep::Int,specialstepexpr=:nothing)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=_nonzero_vals(tabmask)
    dtCijexprs=[:($(Symbol(:dtC,Cind[1],Cind[2]))=$(Symbol(:C,Cind[1],Cind[2]))/dt) for Cind in findall(!iszero,tabmask.C)]
    dtdiexprs=[:($(Symbol(:dtd,dind))=dt*$(Symbol(:d,dind))) for dind in eachindex(tabmask.d)]
    iterexprs=[]
    for i in 1:n_normalstep
        aijkj=[:($(Symbol(:a,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.a[i+1,:])]
        Cijkj=[:($(Symbol(:dtC,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.C[i+1,:])]
        push!(iterexprs,
        quote
            $(Symbol(:k,i)) = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
            integrator.destats.nsolve += 1
            u=+(uprev,$(aijkj...))
            du = f(u, p, t+$(Symbol(:c,i+1))*dt)
            integrator.destats.nf += 1
            linsolve_tmp=+(du,$(Symbol(:dtd,i+1))*dT,$(Cijkj...))
        end)
    end
    push!(iterexprs,specialstepexpr)
    n=length(tabmask.b)
    biki=[:($(Symbol(:b,i))*$(Symbol(:k,i))) for i in 1:n]
    push!(iterexprs,
    quote
        $(Symbol(:k,n))=_reshape(W\-_vec(linsolve_tmp), axes(uprev))
        integrator.destats.nsolve += 1
        u=+(uprev,$(biki...))
    end)

    adaptiveexpr=[]
    if typeof(tabmask)<:RosenbrockAdaptiveTableau
        btildeiki=[:($(Symbol(:btilde,i))*$(Symbol(:k,i))) for i in findall(!iszero,tabmask.btilde)]
        push!(adaptiveexpr,quote
            if integrator.opts.adaptive
                utilde =  +($(btildeiki...))
                atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                                        integrator.opts.reltol,integrator.opts.internalnorm,t)
                integrator.EEst = integrator.opts.internalnorm(atmp,t)
            end
        end)
    end
    quote
        @muladd function perform_step!(integrator, cache::$cachename, repeat_step=false)
            @unpack t,dt,uprev,u,f,p = integrator
            @unpack tf,uf = cache
            $unpacktabexpr

            $(dtCijexprs...)
            $(dtdiexprs...)
            dtgamma = dt*gamma

            # Time derivative
            tf.u = uprev
            dT = ForwardDiff.derivative(tf, t)

            W = calc_W!(integrator, cache, dtgamma, repeat_step, true)
            linsolve_tmp = integrator.fsalfirst + dtd1*dT #calc_rosenbrock_differentiation!

            $(iterexprs...)

            integrator.fsallast = f(u, p, t + dt)
            integrator.destats.nf += 1
          
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
            integrator.u = u

            $(adaptiveexpr...)
        end
    end

end

function gen_perform_step(tabmask::RosenbrockTableau{Bool,Bool},cachename::Symbol,n_normalstep::Int,specialstepexpr=:nothing)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=_nonzero_vals(tabmask)
    dtCij=[:($(Symbol(:dtC,"$(Cind[1])$(Cind[2])"))=$(Symbol(:C,"$(Cind[1])$(Cind[2])"))/dt) for Cind in findall(!iszero,tabmask.C)]
    dtdi=[:($(Symbol(:dtd,dind[1]))=dt*$(Symbol(:d,dind[1]))) for dind in eachindex(tabmask.d)]
    iterexprs=[]
    for i in 1:n_normalstep
        ki=Symbol(:k,i)
        dtdj=Symbol(:dtd,i+1)
        aijkj=[:($(Symbol(:a,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.a[i+1,:])]
        dtCijkj=[:($(Symbol(:dtC,i+1,j))*$(Symbol(:k,j))) for j in findall(!iszero,tabmask.C[i+1,:])]
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
    push!(iterexprs,specialstepexpr)
    n=length(tabmask.b)
    ks=[Symbol(:k,i) for i in 1:n]
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
    end)

    adaptiveexpr=[]
    if typeof(tabmask)<:RosenbrockAdaptiveTableau
        btildeiki=[:($(Symbol(:btilde,i))*$(Symbol(:k,i))) for i in findall(!iszero,tabmask.btilde)]
        push!(adaptiveexpr,quote
            utilde=du
            if integrator.opts.adaptive
                @.. utilde = +($(btildeiki...))
                calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                                    integrator.opts.reltol,integrator.opts.internalnorm,t)
                integrator.EEst = integrator.opts.internalnorm(atmp,t)
            end
        end)
    end
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

            f( fsallast,  u, p, t + dt)
            integrator.destats.nf += 1

            $(adaptiveexpr...)
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
    tabmask=_masktab(tab)
    algname=:RosenbrockW6S4OS
    tabname=:RosenbrockW6S4OSConstantCache
    tabstructname=:RosenbrockW6S4OSConstantCache
    cachename=:RosenbrockW6SCache
    constcachename=:RosenbrockW6SConstantCache
    n_normalstep=length(tab.b)-1
    if part.value==:tableau
        #println("Generating const cache")
        tabstructexpr=gen_tableau_struct(tabmask,tabstructname)
        tabexpr=gen_tableau(tab,tabstructexpr,tabname)
        return esc(quote $([tabstructexpr,tabexpr]...) end)
    elseif part.value==:cache
        #println("Generating cache")
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        algcacheexpr=gen_algcache(cacheexpr,cachename,constcachename,algname,tabname)
        return esc(quote $([constcacheexpr,cacheexpr,algcacheexpr]...) end)
    elseif part.value==:init
        #println("Generating initialize")
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        #println("Generating perform_step")
        constperformstepexpr=gen_constant_perform_step(tabmask,constcachename,n_normalstep)
        performstepexpr=gen_perform_step(tabmask,cachename,n_normalstep)
        return esc(quote $([constperformstepexpr,performstepexpr]...) end)
    else
        println("Unknown parameter!")
    end
end

function Ros4dummyTableau()#create a tabmask for all ROS4 methods where false->0,true->non-0
    a=[false false false;
       true  false false;
       true  true  false]
    C=[false false false false;
       true  false false false;
       true  true  false false;
       true  true  true  false]
    b=[true,true,true,true]
    btilde=[true,true,true,true]
    gamma=true
    c=[false,true,true]
    d=[true,true,true,true]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function RosShamp4Tableau()
    a=[0      0     0;
       2      0     0;
       48//25 6//25 0]
    C=[ 0         0        0    0;
       -8         0        0    0;
        372//25   12//5    0    0;
       -112//125 -54//125 -2//5 0]
    b=[19//9,1//2,25//108,125//108]
    btilde=[17//54,7//36,0,125//108]
    gamma=1//2
    c=[0,1,3//5]
    d=[1//2,-3//2,242//100,116//1000]#2.42,0.116
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function Veldd4Tableau()
    a=[0                 0                 0;
       2                 0                 0;
       4.812234362695436 4.578146956747842 0]
    C=[ 0                  0                  0                 0;
       -5.333333333333331  0                  0                 0;
        6.100529678848254  1.804736797378427  0                 0;
       -2.540515456634749 -9.443746328915205 -1.988471753215993 0]
    b=[4.289339254654537,5.036098482851414,0.6085736420673917,1.355958941201148]
    btilde=[2.175672787531755,2.950911222575741,-0.7859744544887430,-1.355958941201148]
    gamma=0.2257081148225682
    c=[0,0.4514162296451364,0.8755928946018455]
    d=[0.2257081148225682,-0.04599403502680582,0.5177590504944076,-0.03805623938054428]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function Velds4Tableau()
    a=[0    0    0;
       2    0    0;
       7//4 1//4 0]
    C=[ 0     0    0 0;
       -8     0    0 0;
       -8    -1    0 0;
        1//2 -1//2 2 0]
    b=[4//3,2//3,-4//3,4//3]
    btilde=[-1//3,-1//3,0,-4//3]
    gamma=1//2
    c=[0,1,1//2]
    d=[1//2,-3//2,-3//4,1//4]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function GRK4TTableau()
    a=[0                 0                 0;
       2                 0                 0;
       4.524708207373116 4.163528788597648 0]
    C=[ 0                  0                   0                 0;
       -5.071675338776316  0                   0                 0;
        6.020152728650786  0.1597506846727117  0                 0;
       -1.856343618686113 -8.505380858179826  -2.084075136023187 0]
    b=[3.957503746640777,4.624892388363313,0.6174772638750108,1.282612945269037]
    btilde=[2.302155402932996,3.073634485392623,-0.8732808018045032,-1.282612945269037]
    gamma=0.231
    c=[0,0.462,0.8802083333333334]
    d=[0.2310000000000000,-0.03962966775244303,0.5507789395789127,-0.05535098457052764]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function GRK4ATableau()
    a=[0                 0                  0;
       1.108860759493671 0                  0;
       2.377085261983360 0.1850114988899692 0]
    C=[ 0                 0                  0                 0;
       -4.920188402397641 0                  0                 0;
        1.055588686048583 3.351817267668938  0                 0;
        3.846869007049313 3.427109241268180 -2.162408848753263 0]
    b=[1.845683240405840,0.1369796894360503,0.7129097783291559,0.6329113924050632]
    btilde=[0.04831870177201765,-0.6471108651049505,0.2186876660500240,-0.6329113924050632]
    gamma=0.395
    c=[0,0.438,0.87]
    d=[0.395,-0.3726723954840920,0.06629196544571492,0.4340946962568634]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function Ros4LSTableau()
    a=[0                 0                  0;
       2                 0                  0;
       1.867943637803922 0.2344449711399156 0]
    C=[ 0                  0                   0                  0;
       -7.137615036412310  0                   0                  0;
        2.580708087951457  0.6515950076447975  0                  0;
       -2.137148994382534 -0.3214669691237626 -0.6949742501781779 0]
    b=[2.255570073418735,0.2870493262186792,0.4353179431840180,1.093502252409163]
    btilde=[-0.2815431932141155,-0.07276199124938920,-0.1082196201495311,-1.093502252409163]
    gamma=0.5728200000000000
    c=[0,1.145640000000000,0.6552168638155900]
    d=[0.5728200000000000,-1.769193891319233,0.7592633437920482,-0.1049021087100450]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

macro Rosenbrock4(part)
    tabmask=Ros4dummyTableau()#_masktab(tab)
    cachename=:Rosenbrock4Cache
    constcachename=:Rosenbrock4ConstantCache
    RosShamp4tabname=:RosShamp4ConstantCache
    Veldd4tabname=:Veldd4ConstantCache
    Velds4tabname=:Velds4ConstantCache
    GRK4Ttabname=:GRK4TConstantCache
    GRK4Atabname=:GRK4AConstantCache
    Ros4LStabname=:Ros4LStabConstantCache
    n_normalstep=2 #for the third step a4j=a3j which reduced one function call
    if part.value==:tableau
        #println("Generating tableau for Rosenbrock4")
        tabstructexpr=gen_tableau_struct(tabmask,:Ros4ConstantCache)
        tabexprs=Array{Expr,1}()
        push!(tabexprs,tabstructexpr)
        push!(tabexprs,gen_tableau(RosShamp4Tableau(),tabstructexpr,RosShamp4tabname))
        push!(tabexprs,gen_tableau(Veldd4Tableau(),tabstructexpr,Veldd4tabname))
        push!(tabexprs,gen_tableau(Velds4Tableau(),tabstructexpr,Velds4tabname))
        push!(tabexprs,gen_tableau(GRK4TTableau(),tabstructexpr,GRK4Ttabname))
        push!(tabexprs,gen_tableau(GRK4ATableau(),tabstructexpr,GRK4Atabname))
        push!(tabexprs,gen_tableau(Ros4LSTableau(),tabstructexpr,Ros4LStabname))
        return esc(quote $(tabexprs...) end)
    elseif part.value==:cache
        #println("Generating cache for Rosenbrock4")
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        cacheexprs=Array{Expr,1}([constcacheexpr,cacheexpr])
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:RosShamp4,RosShamp4tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:Veldd4,Veldd4tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:Velds4,Velds4tabname))
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:GRK4T,GRK4Ttabname))
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:GRK4A,GRK4Atabname))
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:Ros4LStab,Ros4LStabname))
        return esc(quote $(cacheexprs...) end)
    elseif part.value==:performstep
        #println("Generating perform_step for Rosenbrock4")
        specialstepconst=quote
            k3 = _reshape(W\-_vec(linsolve_tmp), axes(uprev))
            integrator.destats.nsolve += 1
            #u = uprev  + a31*k1 + a32*k2 #a4j=a3j
            #du = f(u, p, t+c3*dt) #reduced function call
            linsolve_tmp =  du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3          
        end
        specialstep=quote
            if DiffEqBase.has_invW(f)
                mul!(vec(k3), W, vec(linsolve_tmp))
            else
                cache.linsolve(vec(k3), W, vec(linsolve_tmp))
                @.. k3 = -k3
            end
            integrator.destats.nsolve += 1
            #@.. u = uprev + a31*k1 + a32*k2 #a4j=a3j
            #f( du,  u, p, t+c3*dt) #reduced function call
            if mass_matrix == I
                @.. linsolve_tmp = du + dtd4*dT + dtC41*k1 + dtC42*k2 + dtC43*k3
            else
                @.. du1 = dtC41*k1 + dtC42*k2 + dtC43*k3
                mul!(du2,mass_matrix,du1)
                @.. linsolve_tmp = du + dtd4*dT + du2
            end
        end
        constperformstepexpr=gen_constant_perform_step(tabmask,constcachename,n_normalstep,specialstepconst)
        performstepexpr=gen_perform_step(tabmask,cachename,n_normalstep,specialstep)
        return esc(quote $([constperformstepexpr,performstepexpr]...) end)
    end
end

#ROS34PW methods (Rang and Angermann, 2005)
function Ros34dummyTableau()
    a=[false false false false;
       true  false false false;
       true  true  false false;
       true  true  true  false]
    C=[false false false false;
       true  false false false;
       true  true  false false;
       true  true  true  false]
    b=[true,true,true,true]
    btilde=[true,true,true,true]
    gamma=true
    c=[false,true,true,true]
    d=[true,true,true,true]
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

function _transformtab(Alpha,Gamma,B,Bhat)
    invGamma=inv(Gamma)
    a=Alpha*invGamma
    C=diagm(0=>diag(invGamma))-invGamma
    b=[(transpose(B)*invGamma)...]# [2Darray...]=>1Darray
    btilde=[(transpose(Bhat)*invGamma)...]
    gamma=Gamma[1,1]#Gamma11==Gamma22==...==Gammass
    d=[sum(Gamma,dims=2)...]#di=sum_j Gamma_ij
    c=[sum(Alpha,dims=2)...]#ci=sum_j Alpha_ij
    (a,C,b,btilde,d,c)
end

function ROS34PW1aTableau()
    gamma=4.358665215084590e-1
    Alpha=[0                 0                    0   0;
           2.218787467653286 0                    0   0;
           0                 0                    0   0; # can reduce one function call with specialized perform_step
           1.208587690772214 7.511610241919324e-2 0.5 0]
    Gamma=[ gamma                 0                    0     0;
           -2.218787467653286     gamma                0     0;
           -9.461966143940745e-2 -7.913526735718213e-3 gamma 0;
           -1.870323744195384    -9.624340112825115e-2 2.726301276675511e-1 gamma]
    B=[3.285609536316354e-1,-5.785609536316354e-1,0.25,1]
    Bhat=[-0.25,0,0.25,1]
    a,C,b,btilde,d,c=_transformtab(Alpha,Gamma,B,Bhat)
    RosenbrockAdaptiveTableau(a,C,b,btilde,gamma,d,c)
end

macro ROS34PW(part)
    tabmask=Ros34dummyTableau()
    cachename=:ROS34PWCache
    constcachename=:ROS34PWConstantCache
    ROS34PW1atabname=:ROS34PW1aConstantCache
    n_normalstep=length(tabmask.b)-1
    if part.value==:tableau
        tabstructexpr=gen_tableau_struct(tabmask,:Ros34ConstantCache)
        tabexprs=Array{Expr,1}()
        push!(tabexprs,tabstructexpr)
        push!(tabexprs,gen_tableau(ROS34PW1aTableau(),tabstructexpr,ROS34PW1atabname))
        return esc(quote $(tabexprs...) end)
    elseif part.value==:cache
        constcacheexpr,cacheexpr=gen_cache_struct(tabmask,cachename,constcachename)
        cacheexprs=Array{Expr,1}([constcacheexpr,cacheexpr])
        push!(cacheexprs,gen_algcache(cacheexpr,cachename,constcachename,:ROS34PW1a,ROS34PW1atabname))
        return esc(quote $(cacheexprs...) end)
    elseif part.value==:init
        return esc(gen_initialize(cachename,constcachename))
    elseif part.value==:performstep
        constperformstepexpr=gen_constant_perform_step(tabmask,constcachename,n_normalstep)
        performstepexpr=gen_perform_step(tabmask,cachename,n_normalstep)
        return esc(quote $([constperformstepexpr,performstepexpr]...) end)
    end
end
