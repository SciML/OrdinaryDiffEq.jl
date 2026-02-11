isfsal(alg::Vern7) = false
isfsal(alg::Vern8) = false
isfsal(alg::Vern9) = false
isfsal(alg::RKV76IIa) = false

alg_order(alg::Vern6) = 6
alg_order(alg::Vern7) = 7
alg_order(alg::Vern8) = 8
alg_order(alg::Vern9) = 9
alg_order(alg::RKV76IIa) = 7

alg_stability_size(alg::Vern6) = 4.8553
alg_stability_size(alg::Vern7) = 4.64
alg_stability_size(alg::Vern8) = 5.8641
alg_stability_size(alg::Vern9) = 4.4762
alg_stability_size(alg::RKV76IIa) = 4.910807773  # From the file: Real Stability Interval is nearly [ -4.910807773, 0]

SciMLBase.has_lazy_interpolation(alg::Union{Vern6, Vern7, Vern8, Vern9, RKV76IIa}) = _unwrap_val(alg.lazy)
