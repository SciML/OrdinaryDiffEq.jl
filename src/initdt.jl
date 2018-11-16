@muladd function ode_determine_initdt(u0,t,tdir,dtmax,abstol,reltol,internalnorm,prob::DiffEqBase.AbstractODEProblem{uType,tType,true},integrator) where {tType,uType}
  nothing
end

@muladd function ode_determine_initdt(u0,t,tdir,dtmax,abstol,reltol,internalnorm,prob::DiffEqBase.AbstractODEProblem{uType,tType,false},integrator) where {uType,tType}
  nothing
end
