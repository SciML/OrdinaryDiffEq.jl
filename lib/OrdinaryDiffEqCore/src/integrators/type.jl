mutable struct DEOptions{
        absType, relType, QT, tType, Controller, F1, F2, F3, F4, F5, F6,
        F7, tstopsType, discType, ECType, SType, MI, tcache, savecache,
        disccache, verbType,
    }
    maxiters::MI
    save_everystep::Bool
    adaptive::Bool
    abstol::absType
    reltol::relType
    # TODO vvv remove this block as these are controller and not integrator parameters vvv
    gamma::QT
    qmax::QT
    qmin::QT
    qsteady_max::QT
    qsteady_min::QT
    qoldinit::QT
    # TODO ^^^ remove this block as these are controller and not integrator parameters ^^^
    failfactor::QT
    dtmax::tType
    dtmin::tType
    controller::Controller
    internalnorm::F1
    internalopnorm::F2
    save_idxs::SType
    tstops::tstopsType
    saveat::tstopsType
    d_discontinuities::discType
    tstops_cache::tcache
    saveat_cache::savecache
    d_discontinuities_cache::disccache
    userdata::ECType
    progress::Bool
    progress_steps::Int
    progress_name::String
    progress_message::F6
    progress_id::Symbol
    timeseries_errors::Bool
    dense_errors::Bool
    dense::Bool
    save_on::Bool
    save_start::Bool
    save_end::Bool
    save_discretes::Bool
    save_end_user::F3
    callback::F4
    isoutofdomain::F5
    unstable_check::F7
    verbose::verbType
    calck::Bool
    force_dtmin::Bool
    advance_to_tstop::Bool
    stop_at_next_tstop::Bool
end

"""
    ODEIntegrator

Fundamental `struct` allowing interactively stepping through the numerical solving of a differential equation.
The full documentation is hosted here:
[https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/).
This docstring describes basic functionality only!

Initialize using `integrator = init(prob::ODEProblem, alg; kwargs...)`. The keyword args which are accepted are the same
[common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
used by `solve`.

For reference, relevant fields of the `ODEIntegrator` are:

  - `t` - time of the proposed step
  - `u` - value at the proposed step
  - `opts` - common solver options
  - `alg` - the algorithm associated with the solution
  - `f` - the function being solved
  - `sol` - the current state of the solution
  - `tprev` - the last timepoint
  - `uprev` - the value at the last timepoint

`opts` holds all of the common solver options, and can be mutated to change the solver characteristics.
For example, to modify the absolute tolerance for the future timesteps, one can do:

```julia
integrator.opts.abstol = 1e-9
```

For more info see the linked documentation page.
"""
mutable struct ODEIntegrator{
        algType <: Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, IIP,
        uType, duType, tType, pType, eigenType, EEstT, QT, tdirType,
        ksEltype, SolType, F, CacheType, O, FSALType, EventErrorType,
        CallbackCacheType, IA, DV, CC, RNGType,
    } <:
    SciMLBase.AbstractODEIntegrator{algType, IIP, uType, tType}
    sol::SolType
    u::uType
    du::duType
    k::ksEltype
    t::tType
    dt::tType
    f::F
    p::pType
    uprev::uType
    uprev2::uType
    duprev::duType
    tprev::tType
    alg::algType
    dtcache::tType
    dtchangeable::Bool
    dtpropose::tType
    tdir::tdirType
    eigen_est::eigenType
    EEst::EEstT
    # TODO vvv remove these
    qold::QT
    q11::QT
    erracc::QT
    dtacc::tType
    # TODO ^^^ remove these
    controller_cache::CC
    success_iter::Int
    iter::Int
    saveiter::Int
    saveiter_dense::Int
    cache::CacheType
    callback_cache::CallbackCacheType
    kshortsize::Int
    force_stepfail::Bool
    last_stepfail::Bool
    just_hit_tstop::Bool
    next_step_tstop::Bool
    tstop_target::tType
    do_error_check::Bool
    event_last_time::Int
    vector_event_last_time::Int
    last_event_error::EventErrorType
    accept_step::Bool
    isout::Bool
    reeval_fsal::Bool
    u_modified::Bool
    reinitialize::Bool
    isdae::Bool
    opts::O
    stats::SciMLBase.DEStats
    initializealg::IA
    differential_vars::DV
    fsalfirst::FSALType
    fsallast::FSALType
    rng::RNGType
end
