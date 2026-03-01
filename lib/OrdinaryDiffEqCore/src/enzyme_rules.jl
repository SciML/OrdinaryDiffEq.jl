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
