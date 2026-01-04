using Markdown
function first_steps(name, solver)
    return Markdown.parse(
        """## Installation
        To be able to access the solvers in `$name`, you must first install them use the Julia package manager:

        ```julia
        using Pkg
        Pkg.add("$name")
        ```
        This will only install the solvers listed at the bottom of this page.
        If you want to explore other solvers for your problem,
        you will need to install some of the other libraries listed in the navigation bar on the left.

        ## Example usage

        ```julia
        using $name

        function lorenz!(du, u, p, t)
            du[1] = 10.0 * (u[2] - u[1])
            du[2] = u[1] * (28.0 - u[3]) - u[2]
            du[3] = u[1] * u[2] - (8 / 3) * u[3]
        end
        u0 = [1.0; 0.0; 0.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(lorenz!, u0, tspan)
        sol = solve(prob, $solver())
        ```"""
    )
end
