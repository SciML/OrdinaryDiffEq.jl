mutable struct SDEIntegrator{
        algType, IIP, uType, uEltype, tType, tdirType, P2, eigenType, tTypeNoUnits,
        uEltypeNoUnits, randType, randType2, rateType, solType, cacheType, F4,
        F6, OType, noiseType, EventErrorType, CallbackCacheType, RCs, IA, RNGType,
        ksEltype,
    } <:
    AbstractSDEIntegrator{algType, IIP, uType, tType}
    f::F4
    c::F6
    noise::noiseType
    uprev::uType
    tprev::tType
    t::tType
    u::uType
    p::P2
    dt::tType
    dtpropose::tType
    dtcache::tType
    T::tType
    tdir::tdirType
    just_hit_tstop::Bool
    do_error_check::Bool
    isout::Bool
    event_last_time::Int
    vector_event_last_time::Int
    last_event_error::EventErrorType
    accept_step::Bool
    last_stepfail::Bool
    force_stepfail::Bool
    dtchangeable::Bool
    u_modified::Bool
    reeval_fsal::Bool
    saveiter::Int
    saveiter_dense::Int
    kshortsize::Int
    k::ksEltype
    alg::algType
    sol::solType
    cache::cacheType
    callback_cache::CallbackCacheType
    sqdt::tType
    W::randType
    P::randType2
    rate_constants::RCs
    opts::OType
    iter::Int
    success_iter::Int
    eigen_est::eigenType
    EEst::tTypeNoUnits
    q::tTypeNoUnits
    qold::tTypeNoUnits
    q11::tTypeNoUnits
    stats::DiffEqBase.Stats
    initializealg::IA
    rng::RNGType
    isdae::Bool
    user_provided_noise::Bool
end
