module OrdinaryDiffEqCoreEnzymeCoreExt
import OrdinaryDiffEqCore, EnzymeCore

Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_nf!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_nf_from_initdt!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_accept!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_reject!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_nf_perform_step!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.increment_nf_fsal!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.check_error!), args...) = true
Enzyme.EnzymeCore.EnzymeRules.inactive(::typeof(OrdinaryDiffEqCore.log_step!), args...) = true

end