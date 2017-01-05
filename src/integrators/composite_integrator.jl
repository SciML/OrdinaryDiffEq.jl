@inline function initialize!(integrator,cache::CompositeCache)
  cache.current = cache.choice_function(integrator)
  initialize!(integrator,cache.caches[cache.current])
  resize!(integrator.k,integrator.kshortsize)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::CompositeCache)
  perform_step!(integrator,cache.caches[cache.current])
end
