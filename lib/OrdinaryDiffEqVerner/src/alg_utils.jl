isfsal(alg::Vern7) = false
isfsal(alg::Vern8) = false
isfsal(alg::Vern9) = false

alg_order(alg::Vern6) = 6
alg_order(alg::Vern7) = 7
alg_order(alg::Vern8) = 8
alg_order(alg::Vern9) = 9

alg_stability_size(alg::Vern6) = 4.8553
alg_stability_size(alg::Vern7) = 4.6400
alg_stability_size(alg::Vern8) = 5.8641
alg_stability_size(alg::Vern9) = 4.4762

has_lazy_interpolation(alg::Union{Vern6,Vern7,Vern8,Vern9}) = true
