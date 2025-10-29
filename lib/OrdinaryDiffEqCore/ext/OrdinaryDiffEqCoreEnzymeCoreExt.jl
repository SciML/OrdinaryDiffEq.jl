module OrdinaryDiffEqCoreEnzymeCoreExt
import OrdinaryDiffEqCore, EnzymeCore

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_nf!), args...)
    true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.fixed_t_for_tstop_error!), args...)
    true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_accept!), args...)
    true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_reject!), args...)
    true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.check_error!), args...)
    true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.log_step!), args...)
    true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.final_progress), args...)
    true
end

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.ode_determine_initdt), args...)
    true
end

end
