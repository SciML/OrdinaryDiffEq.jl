struct DummyController <: AbstractController
end

function default_controller(alg::Union{QNDF, FBDF}, args...)
    DummyController()
end