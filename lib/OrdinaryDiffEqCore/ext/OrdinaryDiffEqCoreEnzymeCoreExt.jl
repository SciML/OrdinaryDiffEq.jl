module OrdinaryDiffEqCoreEnzymeCoreExt
import OrdinaryDiffEqCore, EnzymeCore

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_nf!), args...
    )
    return true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!), args...
    )
    return true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_accept!), args...
    )
    return true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_reject!), args...
    )
    return true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.check_error!), args...
    )
    return true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.log_step!), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.final_progress), args...
    )
    return true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.ode_determine_initdt), args...
    )
    return true
end

end
