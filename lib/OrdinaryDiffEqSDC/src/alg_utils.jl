isfsal(::SDC) = false

"""
    sdc_iteration_order(alg)

Order reached by the SDC iteration itself, before the collocation cap.

One order is gained per sweep, two on the first sweep for a second-order
sweeper, and one more if the step update is the quadrature rule rather than the
last node value. This is `qmat.utils.sdc.getOrderSDC` without its table of
empirically found bonus cases, so for some node/sweeper combinations (notably
`SDCSweeper.LU` and `SDCSweeper.MIN_SR_NS`) the observed order is higher than
this predicts.
"""
function sdc_iteration_order(alg::SDC)
    K = alg.num_sweeps
    order = 0
    if K > 0
        order += alg.sweeper in SDC_SECOND_ORDER_SWEEPERS ? 2 : 1
        order += K - 1
    end
    alg.step_update === SDCStepUpdate.Quadrature && (order += 1)
    return order
end

function alg_order(alg::SDC)
    coll = sdc_collocation_order(alg.num_nodes, alg.node_type, alg.quad_type)
    return max(1, min(sdc_iteration_order(alg), coll))
end

"""
    sdc_validate(num_nodes, quad_type, num_sweeps, step_update)

Reject unusable parameter combinations at algorithm construction time.

This deliberately does not go through `prepare_alg`: the generic
`prepare_alg(::OrdinaryDiffEqImplicitAlgorithm, ::AbstractArray, …)` in
`OrdinaryDiffEqDifferentiation` does the autodiff preparation every implicit
solver needs, and adding an `SDC` method here would make the two ambiguous.
"""
function sdc_validate(
        num_nodes::Int, quad_type::SDCQuadrature.T,
        num_sweeps::Int, step_update::SDCStepUpdate.T
    )
    num_sweeps >= 0 ||
        throw(ArgumentError("SDC: `num_sweeps` must be ≥ 0, got $(num_sweeps)"))
    endpoint_on_left = quad_type in (SDCQuadrature.Lobatto, SDCQuadrature.RadauLeft)
    minimum_nodes = endpoint_on_left ? 2 : 1
    num_nodes >= minimum_nodes || throw(
        ArgumentError(
            "SDC: `num_nodes` must be ≥ $(minimum_nodes) for `quad_type = $(quad_type)`, " *
                "got $(num_nodes)"
        )
    )
    endpoint_on_right = quad_type in (SDCQuadrature.RadauRight, SDCQuadrature.Lobatto)
    if step_update === SDCStepUpdate.LastNode
        endpoint_on_right || throw(
            ArgumentError(
                "SDC: `step_update = SDCStepUpdate.LastNode` needs the last node to be " *
                    "the right endpoint, so `quad_type` must be `SDCQuadrature.RadauRight` " *
                    "or `SDCQuadrature.Lobatto`, got $(quad_type)"
            )
        )
        # With no sweep the node values are still the copy initialisation, so the
        # last of them is uₙ and the step would return the initial value forever.
        num_sweeps >= 1 || throw(
            ArgumentError(
                "SDC: `step_update = SDCStepUpdate.LastNode` needs `num_sweeps ≥ 1`, " *
                    "got $(num_sweeps)"
            )
        )
    end
    return nothing
end
