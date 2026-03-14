# unwrap_cache(::SDEIntegrator, is_stiff) is now handled by ODE's generic
# SciMLBase.unwrap_cache(::ODEIntegrator, is_stiff) which uses the
# is_composite_algorithm trait instead of checking `alg isa CompositeAlgorithm`.
