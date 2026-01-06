function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::FunctionMapCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    return nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::FunctionMapConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    return nothing
end
