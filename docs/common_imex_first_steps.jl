using Markdown
function imex_first_steps(name, solver)
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

        IMEX methods require a `SplitODEProblem` that separates the differential equation into two parts:
        one for stiff (implicit) terms and one for non-stiff (explicit) terms.

        ```julia
        using $name

        # Example: A 1D reaction-diffusion system
        # du/dt = f₁(u,p,t) + f₂(u,p,t)
        # where f₁ contains stiff linear diffusion (implicit)
        # and f₂ contains non-stiff reaction terms (explicit)

        function f1!(du, u, p, t)
            # Stiff diffusion term: D * ∂²u/∂x²
            # Using simple finite differences
            D = 0.01
            N = length(u)
            du[1] = D * (u[2] - 2u[1] + u[end])
            for i in 2:N-1
                du[i] = D * (u[i+1] - 2u[i] + u[i-1])
            end
            du[N] = D * (u[1] - 2u[N] + u[N-1])
        end

        function f2!(du, u, p, t)
            # Non-stiff reaction term: u * (1 - u)
            @. du = u * (1 - u)
        end

        u0 = 0.5 .+ 0.1 * sin.(2π * (0:0.1:0.9))  # Initial condition with perturbation
        tspan = (0.0, 10.0)
        prob = SplitODEProblem(f1!, f2!, u0, tspan)
        sol = solve(prob, $solver(), dt = 0.01)
        ```

        ```"""
    )
end
