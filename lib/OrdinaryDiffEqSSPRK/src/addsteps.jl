function _ode_addsteps!(k, t, uprev, u, dt, f, p,
    cache::Union{SSPRK22ConstantCache, SSPRK33ConstantCache,
        SSPRK43ConstantCache, SSPRK432ConstantCache},
    always_calc_begin = false, allow_calc_end = true,
    force_calc_end = false)
    if length(k) < 1 || always_calc_begin
        copyat_or_push!(k, 1, f(uprev, p, t))
    end
    nothing
end

function _ode_addsteps!(k, t, uprev, u, dt, f, p,
        cache::Union{SSPRK22Cache, SSPRK33Cache, SSPRK43Cache,
            SSPRK432Cache}, always_calc_begin = false,
        allow_calc_end = true, force_calc_end = false)
    if length(k) < 1 || always_calc_begin
        f(cache.k, uprev, p, t)
        copyat_or_push!(k, 1, cache.k)
    end
    nothing
end