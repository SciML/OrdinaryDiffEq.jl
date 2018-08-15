is_constant_cache(cache::OrdinaryDiffEqConstantCache) = true
is_constant_cache(cache::OrdinaryDiffEqCache) = false
is_constant_cache(cache::CompositeCache) = is_constant_cache(cache.caches[1])
