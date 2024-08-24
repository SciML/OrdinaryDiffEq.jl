module OrdinaryDiffEqCoreEnzymeCoreExt
import OrdinaryDiffEqCore, EnzymeCore

EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_nf!), args...) = true
EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!), args...) = true
EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_accept!), args...) = true
EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_reject!), args...) = true
EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_nf_perform_step!), args...) = true
EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.check_error!), args...) = true
EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.log_step!), args...) = true

end
