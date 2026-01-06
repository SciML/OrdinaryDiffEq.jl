@doc explicit_rk_docstring(
    "Tanaka-Yamashita 7 Runge-Kutta method.",
    "TanYam7",
    references = "Tanaka M., Muramatsu S., Yamashita S., (1992), On the Optimization of Some Nine-Stage
    Seventh-order Runge-Kutta Method, Information Processing Society of Japan,
    33 (12), pp. 1512-1526."
)
Base.@kwdef struct TanYam7{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function TanYam7(stage_limiter!, step_limiter! = trivial_limiter!)
    return TanYam7(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Tsitouras-Papakostas 8/7 Runge-Kutta method.", "TsitPap8",
    references = """@article{tsitouras1999cheap,
    title={Cheap error estimation for Runge--Kutta methods},
    author={Tsitouras, Ch and Papakostas, SN},
    journal={SIAM Journal on Scientific Computing},
    volume={20},
    number={6},
    pages={2067--2088},
    year={1999},
    publisher={SIAM}}"""
)
Base.@kwdef struct TsitPap8{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function TsitPap8(stage_limiter!, step_limiter! = trivial_limiter!)
    return TsitPap8(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Hairer's 8/5/3 adaption of the Dormand-Prince Runge-Kutta method.",
    "DP8",
    references = "E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
    Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
    Springer-Verlag."
)
Base.@kwdef struct DP8{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function DP8(stage_limiter!, step_limiter! = trivial_limiter!)
    return DP8(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Phase-fitted Runge-Kutta of 8th order.", "PFRK87",
    references = """@article{tsitouras2017phase,
    title={Phase-fitted Runge--Kutta pairs of orders 8 (7)},
    author={Tsitouras, Ch and Famelis, I Th and Simos, TE},
    journal={Journal of Computational and Applied Mathematics},
    volume={321},
    pages={226--231},
    year={2017},
    publisher={Elsevier}}""",
    extra_keyword_description = """- `omega`: a periodicity phase estimate,
                   when accurate this method results in zero numerical dissipation.
    """,
    extra_keyword_default = "omega = 0.0"
)
Base.@kwdef struct PFRK87{StageLimiter, StepLimiter, Thread, T} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    omega::T = 0.0
end
# for backwards compatibility
function PFRK87(stage_limiter!, step_limiter! = trivial_limiter!; omega = 0.0)
    return PFRK87(stage_limiter!, step_limiter!, False(), omega)
end
