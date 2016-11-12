"""
constructDormandPrince()

Constructs the tableau object for the Dormand-Prince Order 4/5 method.
"""
function constructDormandPrince(T::Type = Float64)
  A = [0 0 0 0 0 0 0
      1//5 0 0 0 0 0 0
      3//40 9//40 0 0 0 0 0
      44//45 -56//15 32//9 0 0 0 0
      19372//6561 -25360//2187 64448//6561 -212//729 0 0 0
      9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
      35//384 0 500//1113 125//192 -2187//6784 11//84 0]
  c = [0;1//5;3//10;4//5;8//9;1;1]
  α = [35//384;0;500//1113;125//192;-2187//6784;11//84;0]
  αEEst = [5179//57600;0;7571//16695;393//640;-92097//339200;187//2100;1//40]
  A = map(T,A)
  α = map(T,α)
  αEEst = map(T,αEEst)
  c = map(T,c)
  return(ExplicitRKTableau(A,c,α,5,αEEst=αEEst,adaptiveorder=4,fsal=true))
end

"""
ODE_DEFAULT_TABLEAU

Sets the default tableau for the ODE solver. Currently Dormand-Prince 4/5.
"""
const ODE_DEFAULT_TABLEAU = constructDormandPrince()

"""
ODE_DEFAULT_CALLBACK

All it does is call the saving functionality.
"""
@inline function ODE_DEFAULT_CALLBACK(alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,dtprev,dt,
  saveat,cursaveat,saveiter,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,
  dense,kshortsize,issimple_dense,fsal,fsalfirst,cache,calck,T,Ts)
  @ode_savevalues
  reeval_fsal = false
  cursaveat,saveiter,dt,t,T,reeval_fsal
end

@inline ODE_DEFAULT_NORM(u) = sqrt(sumabs2(u) / length(u))
@inline ODE_DEFAULT_ISOUTOFDOMAIN(t,u) = false
