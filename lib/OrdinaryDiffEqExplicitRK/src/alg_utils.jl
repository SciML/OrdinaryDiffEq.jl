isfsal(tab::DiffEqBase.ExplicitRKTableau) = tab.fsal

alg_order(alg::ExplicitRK) = alg.tableau.order

alg_adaptive_order(alg::ExplicitRK) = alg.tableau.adaptiveorder

alg_stability_size(alg::ExplicitRK) = alg.tableau.stability_size