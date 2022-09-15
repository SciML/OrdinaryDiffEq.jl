using Documenter, OrdinaryDiffEq

makedocs(sitename = "OrdinaryDiffEq.jl",
         authors = "Chris Rackauckas et al.",
         clean = true,
         doctest = false,
         modules = [OrdinaryDiffEq],
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://ordinarydiffeq.sciml.ai/stable/"),
         pages = [
             "OrdinaryDiffEq.jl: ODE solvers and utilities" => "index.md",
             "Usage" => "usage.md",
             "Standard Non-Stiff ODEProblem Solvers" => [
                "explicitrk.md",
                "lowstorage_ssprk.md"
                "explicit_extrapolation.md",
                "nonstiff_multistep.md",
             ],
             "Standard Stiff ODEProblem Solvers" => [
                "firk.md",
                "rosenbrock.md",
                "stabilized_rk.md",
                "sdirk.md",
                "stiff_multistep.md",
                "implicit_extrapolation.md"
             ],
             "Second Order and Dynamical ODE Solvers" => [
                "nystrom.md",
                "symplectic.md"
             ],
             "IMEX Solvers" => [
                "imex_multistep.md",
                "imex_sdirk.md",
             ],
             "Semilinear ODE Solvers" => [
                "exponential_rk.md",
                "magnus.md"
             ],
             "DAEProblem Solvers" => [
                "fully_implicit.md"
             ],
             "Misc Solvers" => [
                "misc.md"
             ]
         ])

deploydocs(repo = "github.com/SciML/OrdinaryDiffEq.jl";
           push_preview = true)
