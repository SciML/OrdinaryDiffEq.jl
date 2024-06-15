abstract type AbstractPerformStepCallback end

struct NoPerformStepCallback <:AbstractPerformStepCallback end

has_poststep_callback(::AbstractPerformStepCallback) = false
has_postfsal_callback(::AbstractPerformStepCallback) = false
has_postEEst_callback(::AbstractPerformStepCallback) = false

struct PerformStepCallback <: AbstractPerformStepCallback
    poststep_cb
    has_poststep_cb::Bool
    postfsal_cb
    has_postfsal_cb::Bool
    postEEst_cb
    has_postEEst_cb::Bool

    function PerformStepCallback(;poststep = nothing, postfsal = nothing, postEEst = nothing)
        has_poststep_cb = poststep isa Nothing ? false : true
        has_postfsal_cb = postfsal isa Nothing ? false : true
        has_postEEst_cb = postEEst isa Nothing ? false : true
        new(poststep, has_poststep_cb, postfsal, has_postfsal_cb, postEEst, has_postEEst_cb)
    end
end

has_poststep_callback(pscb::PerformStepCallback) = pscb.has_poststep_cb
has_postfsal_callback(pscb::PerformStepCallback) = pscb.has_postfsal_cb
has_postEEst_callback(pscb::PerformStepCallback) = pscb.has_postEEst_cb

