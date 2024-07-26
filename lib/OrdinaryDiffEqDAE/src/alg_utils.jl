isadaptive(alg::DImplicitEuler) = true
isadaptive(alg::DABDF2) = true
isadaptive(alg::DFBDF) = true

alg_extrapolates(alg::DImplicitEuler) = true
alg_extrapolates(alg::DABDF2) = true

alg_order(alg::DImplicitEuler) = 1
alg_order(alg::DABDF2) = 2
alg_order(alg::DFBDF) = 1#dummy value

isfsal(alg::DImplicitEuler) = false