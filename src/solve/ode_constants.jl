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

@inline function ODE_DEFAULT_NORM(u)
  for i in eachindex(u)
    u[i] = u[i]*u[i]
  end
  sqrt( sum(u) / length(u))
end

const ODEJL_OPTION_LIST = Set([:tout,:tstop,:reltol,:abstol,:minstep,:maxstep,:initstep,:norm,:maxiters,:isoutofdomain])
const ODEJL_ALIASES = Dict{Symbol,Symbol}(:minstep=>:dtmin,:maxstep=>:dtmax,:initstep=>:dt,:tstop=>:T,:maxiters=>:maxiters)
const ODEJL_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEJL_ALIASES])

const ODEINTERFACE_OPTION_LIST = Set([:RTOL,:ATOL,:OUTPUTFCN,:OUTPUTMODE,
                                :MAXSTEPS,:STEST,:EPS,:RHO,:SSMINSEL,
                                :SSMAXSEL,:SSBETA,:MAXSS,:INITIALSS,
                                :MAXEXCOLUMN,:STEPSIZESEQUENCE,:MAXSTABCHECKS,
                                :MAXSTABCHECKLINE,:DENSEOUTPUTWOEE,:INTERPOLDEGRE,
                                :SSREDUCTION,:SSSELECTPAR1,:SSSELECTPAR2,
                                :ORDERDECFRAC,:ORDERINCFRAC,:OPT_RHO,:OPT_RHO2,
                                :RHSAUTONOMOUS,:M1,:M2,:LAMBDADENSE,:TRANSJTOH,
                                :STEPSIZESEQUENCE,:JACRECOMPFACTOR,:MASSMATRIX,
                                :JACOBIMATRIX,:JACOBIBANDSSTRUCT,:WORKFORRHS,
                                :WORKFORJAC,:WORKFORDEC,:WORKFORSOL,:MAXNEWTONITER,
                                :NEWTONSTARTZERO,:NEWTONSTOPCRIT,:DIMFIND1VAR,
                                :MAXSTAGES,:MINSTAGES,:INITSTAGES,:STEPSIZESTRATEGY,
                                :FREEZESSLEFT,:FREEZESSRIGHT,:ORDERDECFACTOR,
                                :ORDERINCFACTOR,:ORDERDECCSTEPFAC1,:ORDERDECSTEPFAC2,
                                :RHS_CALLMODE
                                ])

const ODEINTERFACE_ALIASES = Dict{Symbol,Symbol}(:RTOL=>:reltol,
                                                 :ATOL=>:abstol,
                                                 :MAXSTEPS=> :maxiters,
                                                 :MAXSS=>:dtmax,
                                                 :INITIALSS=>:dt,
                                                 #:SSMINSEL=>:qmin,
                                                 :SSBETA=>:beta2,
                                                 :SSMAXSEL=>:qmax)

const ODEINTERFACE_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEINTERFACE_ALIASES])

const ODEINTERFACE_STRINGS = Dict{Symbol,String}(
  :LOGIO            => "logio",
  :LOGLEVEL         => "loglevel",
  :RHS_CALLMODE     => "RightHandSideCallMode",

  :RTOL             => "RelTol",
  :ATOL             => "AbsTol",
  :MAXSTEPS         => "MaxNumberOfSteps",
  :EPS              => "eps",

  :OUTPUTFCN        => "OutputFcn",
  :OUTPUTMODE       => "OutputFcnMode",

  :STEST            => "StiffTestAfterStep",
  :RHO              => "rho",
  :SSMINSEL         => "StepSizeMinSelection",
  :SSMAXSEL         => "StepSizeMaxSelection",
  :SSBETA           => "StepSizeBeta",
  :MAXSS            => "MaxStep",
  :INITIALSS        => "InitialStep",


  :MAXEXCOLUMN      => "MaxExtrapolationColumn",
  :MAXSTABCHECKS    => "MaxNumberOfStabilityChecks",
  :MAXSTABCHECKLINE => "MaxLineForStabilityCheck",
  :INTERPOLDEGREE   => "DegreeOfInterpolation",
  :ORDERDECFRAC     => "OrderDecreaseFraction",
  :ORDERINCFRAC     => "OrderIncreaseFraction",
  :STEPSIZESEQUENCE => "StepSizeSequence",
  :SSREDUCTION      => "StepSizeReduction",
  :SSSELECTPAR1     => "StepSizeSelectionParam1",
  :SSSELECTPAR2     => "StepSizeSelectionParam2",
  :RHO2             => "rho2",
  :DENSEOUTPUTWOEE  => "DeactivateErrorEstInDenseOutput",

  :TRANSJTOH        => "TransfromJACtoHess",
  :MAXNEWTONITER    => "MaxNewtonIterations",
  :NEWTONSTARTZERO  => "StartNewtonWithZeros",
  :DIMOFIND1VAR     => "DimensionOfIndex1Vars",
  :DIMOFIND2VAR     => "DimensionOfIndex2Vars",
  :DIMOFIND3VAR     => "DimensionOfIndex3Vars",
  :STEPSIZESTRATEGY => "StepSizeStrategy",
  :M1               => "M1",
  :M2               => "M2",
  :JACRECOMPFACTOR  => "RecomputeJACFactor",
  :NEWTONSTOPCRIT   => "NewtonStopCriterion",
  :FREEZESSLEFT     => "FreezeStepSizeLeftBound",
  :FREEZESSRIGHT    => "FreezeStepSizeRightBound",
  :MASSMATRIX       => "MassMatrix",
  :JACOBIMATRIX     => "JacobiMatrix",
  :JACOBIBANDSTRUCT => "JacobiBandStructure",

  :MAXSTAGES        => "MaximalNumberOfStages",
  :MINSTAGES        => "MinimalNumberOfStages",
  :INITSTAGES       => "InitialNumberOfStages",
  :ORDERINCFACTOR   => "OrderIncreaseFactor",
  :ORDERDECFACTOR   => "OrderDecreaseFactor",
  :ORDERDECSTEPFAC1 => "OrderDecreaseStepFactor1",
  :ORDERDECSTEPFAC2 => "OrderDecreaseStepFactor2",

  :RHSAUTONOMOUS    => "AutonomousRHS",
  :LAMBDADENSE      => "LambdaForDenseOutput",
  :WORKFORRHS       => "WorkForRightHandSide",
  :WORKFORJAC       => "WorkForJacobimatrix",
  :WORKFORDEC       => "WorkForLuDecomposition",
  :WORKFORSOL       => "WorkForSubstitution",

  :BVPCLASS         => "BoundaryValueProblemClass",
  :SOLMETHOD        => "SolutionMethod",
  :IVPOPT           => "OptionsForIVPsolver")
