alg_order(alg::Euler) = 1
alg_order(alg::SplitEuler) = 1
alg_order(alg::Heun) = 2
alg_order(alg::Ralston) = 2
alg_order(alg::Midpoint) = 2
alg_order(alg::BS3) = 3
alg_order(alg::OwrenZen3) = 3
alg_order(alg::BS5) = 5
alg_order(alg::OwrenZen4) = 4
alg_order(alg::OwrenZen5) = 5
alg_order(alg::DP5) = 5
alg_order(alg::Anas5) = 5
alg_order(alg::RKO65) = 5
alg_order(alg::FRK65) = 6
alg_order(alg::RK4) = 4
alg_order(alg::RKM) = 4
alg_order(alg::MSRK5) = 5
alg_order(alg::MSRK6) = 6
alg_order(alg::PSRK4p7q6) = 4
alg_order(alg::PSRK3p6q5) = 3
alg_order(alg::PSRK3p5q4) = 3
alg_order(alg::Stepanov5) = 5
alg_order(alg::SIR54) = 5
alg_order(alg::Alshina2) = 2
alg_order(alg::Alshina3) = 3
alg_order(alg::Alshina6) = 6

isfsal(alg::FRK65) = true
isfsal(alg::RKO65) = false
isfsal(alg::PSRK3p5q4) = false
isfsal(alg::PSRK3p6q5) = false
isfsal(alg::PSRK4p7q6) = false

beta2_default(alg::DP5) = 4 // 100

beta1_default(alg::DP5, beta2) = typeof(beta2)(1 // alg_order(alg)) - 3beta2 / 4

alg_stability_size(alg::DP5) = 3.3066

ssp_coefficient(alg::Euler) = 1

function DiffEqBase.prepare_alg(
        alg::SplitEuler,
        u0::AbstractArray,
        p, prob)
    alg
end