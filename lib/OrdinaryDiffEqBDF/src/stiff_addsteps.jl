function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::Union{
            QNDFConstantCache, QNDFCache,
            FBDFConstantCache, FBDFCache,
            DFBDFConstantCache, DFBDFCache,
        },
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    # Backward differences are precomputed in perform_step! and saved via copyat_or_push!,
    # so no additional computation is needed here.
    return nothing
end
