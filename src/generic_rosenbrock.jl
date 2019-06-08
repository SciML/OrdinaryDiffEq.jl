struct RosenbrockTableau{T,T2}
    a::Array{T,2}
    C::Array{T,2}
    b::Array{T,1}
    gamma::T2
    d::Array{T,1}
    c::Array{T2,1}
end

function _gen_constant_cache!(cachevalexpr,inds,name,T,arr)
    for ind in inds
        if length(size(arr))==2
            indstr="$(ind[1])$(ind[2])"
        elseif length(size(arr))==1
            indstr="$(ind[1])"
        end
        push!(cachevalexpr.args,:($(Symbol(name,indstr))=convert($T,$(arr[ind]))))
    end
end

function gen_constant_cache(tab::RosenbrockTableau,cachename::Symbol)
    ainds=findall(!iszero,tab.a)
    Cinds=findall(!iszero,tab.C)
    cachevalexpr=quote end
    _gen_constant_cache!(cachevalexpr,ainds,:a,:T,tab.a)
    _gen_constant_cache!(cachevalexpr,Cinds,:C,:T,tab.C)
    _gen_constant_cache!(cachevalexpr,eachindex(tab.b),:b,:T,tab.b)
    push!(cachevalexpr.args,:(gamma=convert(T2,$(tab.gamma))))
    _gen_constant_cache!(cachevalexpr,eachindex(tab.d),:d,:T,tab.d)
    _gen_constant_cache!(cachevalexpr,findall(!iszero,tab.c),:c,:T2,tab.c)
    
    cachestructexpr=quote end
    cacheinitexpr=:($cachename())
    for lineexpr in cachevalexpr.args
        if typeof(lineexpr)==LineNumberNode
            continue
        end
        #e.g :(c6=convert(T2,0.92))
        varsym=lineexpr.args[1]#head(=)->args1
        typesym=lineexpr.args[2].args[2]#head(=)->args2 head(call)->args2
        push!(cachestructexpr.args,:($varsym::$typesym))
        push!(cacheinitexpr.args,varsym)
    end
    push!(cachevalexpr.args,cacheinitexpr)
    cacheexpr=quote
        struct $cachename{T,T2} <: OrdinaryDiffEqConstantCache
            $cachestructexpr
        end
        function $cachename(T::Type,T2::Type)
            $cachevalexpr
        end
    end
    cacheinitexpr,cacheexpr
end

function gen_constant_perform_step(tab,cachename,cacheinitexpr)
    unpacktabexpr=:(@unpack ()=cache.tab)
    unpacktabexpr.args[3].args[1].args=copy(cacheinitexpr.args[2:end])
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
#=
    for aind in ainds
        push!(cachevalexpr.args,:($(Symbol("a$(aind[1])$(aind[2])"))=convert(T,$(tab.a[aind]))))
    end
    for Cind in cinds
        push!(cachevalexpr.args,:($(Symbol("C$(cind[1])$(cind[2])"))=convert(T,$(tab.C[Cind]))))
    end
=#
#---------------Constant Cache--------------

macro RosenbrockW6S4OSConstant(part)
    tab=RosenbrockW6S4OSTableau()
    cacheinitexpr,cacheexpr=gen_constant_cache(tab,:RosenbrockW6S4OSConstantCache)
    if part.value==:cache
        println("Generating cache")
        return esc(cacheexpr)
    elseif part.value==:performstep
        println("Generating perform_step")
        performstepexpr=gen_constant_perform_step(tab,:RosenbrockWConstantCache,cacheinitexpr)
        return esc(performstepexpr)
    else
        println("Unknown parameter!")
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
    RosenbrockTableau(a,C,b,gamma,d,c)
end