function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(increment_nf!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(fixed_t_for_tstop_error!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(increment_accept!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(increment_reject!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(check_error!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(log_step!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(final_progress), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(ode_determine_initdt), args...
    )
    return true
end

# Discrete time-grid index helpers have zero derivative; block inlining so
# Enzyme treats the index as constant.
function EnzymeCore.EnzymeRules.inactive(
        ::typeof(_searchsortedfirst), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive(
        ::typeof(_searchsortedlast), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive(
        ::typeof(_ts_grid_kind), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive(
        ::typeof(ts_hint_start), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive(
        ::typeof(reprobe_ts_hint!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive(
        ::typeof(maybe_reprobe_ts_hint!), args...
    )
    return true
end
