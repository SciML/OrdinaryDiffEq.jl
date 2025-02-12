module OrdinaryDiffEqCoreEnzymeCoreExt
import OrdinaryDiffEqCore, EnzymeCore

function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.increment_nf!), args...)
    true
end
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!), args...)
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

end
