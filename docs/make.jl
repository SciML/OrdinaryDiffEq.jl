using Documenter, OrdinaryDiffEq

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

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
                 "nonstiff/explicitrk.md",
                 "nonstiff/lowstorage_ssprk.md",
                 "nonstiff/explicit_extrapolation.md",
                 "nonstiff/nonstiff_multistep.md",
             ],
             "Standard Stiff ODEProblem Solvers" => [
                 "stiff/firk.md",
                 "stiff/rosenbrock.md",
                 "stiff/stabilized_rk.md",
                 "stiff/sdirk.md",
                 "stiff/stiff_multistep.md",
                 "stiff/implicit_extrapolation.md",
             ],
             "Second Order and Dynamical ODE Solvers" => [
                 "dynamical/nystrom.md",
                 "dynamical/symplectic.md",
             ],
             "IMEX Solvers" => [
                 "imex/imex_multistep.md",
                 "imex/imex_sdirk.md",
             ],
             "Semilinear ODE Solvers" => [
                 "semilinear/exponential_rk.md",
                 "semilinear/magnus.md",
             ],
             "DAEProblem Solvers" => [
                 "dae/fully_implicit.md",
             ],
             "Misc Solvers" => [
                 "misc.md",
             ],
         ])

deploydocs(repo = "github.com/SciML/OrdinaryDiffEq.jl";
           push_preview = true)
