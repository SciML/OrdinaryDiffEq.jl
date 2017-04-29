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

@inline UNITLESS_ABS2(x) = abs2(x)/(typeof(x)(one(x))*typeof(x)(one(x)))
@inline ODE_DEFAULT_NORM(u) = sqrt(sum(UNITLESS_ABS2,u) / length(u))
@inline ODE_DEFAULT_NORM(u::Number) = abs(u)
@inline ODE_DEFAULT_ISOUTOFDOMAIN(t,u) = false
@inline ODE_DEFAULT_PROG_MESSAGE(dt,t,u) = "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,t,u) = false
(p::typeof(ODE_DEFAULT_UNSTABLE_CHECK)){T<:AbstractFloat}(dt,t,u::AbstractArray{T}) = any(isnan,u)
(p::typeof(ODE_DEFAULT_UNSTABLE_CHECK))(dt,t,u::ArrayPartition) = any(isnan,chain(u.x...))
