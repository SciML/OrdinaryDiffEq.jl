using Symbolics
using SymbolicUtils
using Hungarian
using JuMP
using GLPK

export compute_signature_matrix, compute_hvt, compute_offsets, compute_jacobian
"""
    signature_matrix(eqs, vars, t_ind) -> Matrix{Float64}

Construct the signature matrix Σ for the system of equations `eqs` with respect to
the variables `vars`. Each entry Σ[i,j] is the (highest) derivative order of `vars[j]`
appearing in eqs[i], or -Inf if it does not appear.
"""
function signature_matrix(eqs::Vector, vars::Vector, t_ind)
    # eqs[i] is something like f_i(x₁(t), x₂'(t), x₃''(t), ...)
    n_eqs   = length(eqs)
    n_vars  = length(vars)
    Σ = fill(-Inf, n_eqs, n_vars) # Initialize with -Inf

    for i in 1:n_eqs
        # For each equation
        for j in 1:n_vars
            # Check for each variable in the equation what the highest derivative order is
            order = max_derivative_order(eqs[i], vars[j], t_ind)
            Σ[i,j] = order
        end
    end
    return Σ
end
"""
    max_derivative_order(ex, var, t_ind)

Returns the highest derivative order of `var` that appears in the symbolic expression `ex`,
using `t_ind` as the independent variable (e.g., time). If `var` does not appear, returns -Inf.
"""
function max_derivative_order(ex, var, t_ind)
    # Base case: if ex is exactly var(t_ind), order is 0.
    if isequal(ex, var(t_ind))
        return 0
    end

    # If it's a number or an unrelated symbol we ignore
    if ex isa Number || (ex isa Symbol && ex != var)
        return -Inf
    end

    # Check if ex is a derivative expression.
    if iscall(ex) && operation(ex) isa Symbolics.Differential
        inner = arguments(ex)[1]  # The function being differentiated.
        # Recursively check the order of the inner expression.
        sub_order = max_derivative_order(inner, var, t_ind)
        # If the inner expression is not related to var, return -Inf otherwise add 1 to derivative order.
        return sub_order == -Inf ? -Inf : 1 + sub_order
    end

    # For composite expressions (e.g., sums, products), traverse their arguments.
    if iscall(ex)
        best = -Inf
        for arg in arguments(ex)
            # Recursively check the order of each component
            best = max(best, max_derivative_order(arg, var, t_ind))
        end
        return best
    end
    return -Inf
end

"""
    highest_value_transversal(Σ) -> (Vector{Tuple{Int, Int}}, Float64)

Finds the highest value transversal (HVT) of the signature matrix `Σ` using the Hungarian algorithm.
Returns the transversal as a vector of (row, column) indices and its value.
"""
function highest_value_transversal(Σ::Matrix{Float64})
    n = size(Σ, 1)
    
    # The Hungarian algorithm minimizes so multiply by -1 to max
    cost_matrix = -Σ
    assignment = hungarian(cost_matrix)[1]
    # Extract transversal and its value.
    transversal = [(i, assignment[i]) for i in 1:n]
    value = sum(Σ[i, assignment[i]] for i in 1:n)
    
    return transversal, value
end


"""
    find_offsets(Σ::Matrix{Float64}, T::Vector{Tuple{Int, Int}}) -> (Vector{Int}, Vector{Int})

Finds the canonical offsets `c` and `d` for the signature matrix `Σ` and the highest value transversal `T`.
Returns the vectors `c` and `d`.
"""
function find_offsets(Σ::Matrix{Float64}, T::Vector{Tuple{Int, Int}})
    n = size(Σ, 1)
    
    # Create JuMP model with GLPK solver
    model = Model(GLPK.Optimizer)
    
    # Define variables for c and d
    @variable(model, c[1:n] >= 0, Int)
    @variable(model, d[1:n] >= 0, Int)
    
    # Add constraints for all i, j: d[j] - c[i] >= Σ[i, j]
    for i in 1:n
        for j in 1:n
            if Σ[i, j] != -Inf
                @constraint(model, d[j] - c[i] >= Σ[i, j])
            end
        end
    end
    
    # Add constraints for equality over transversal
    for (i, j) in T
        @constraint(model, d[j] - c[i] == Σ[i, j])
    end
    
    # min sum c and d to find canonical offsets
    @objective(model, Min, sum(c) + sum(d))
    optimize!(model)
    c_values = value.(c)
    d_values = value.(d)
    return c_values, d_values
end


"""
    system_jacobian(eqs::Vector{SymbolicUtils.BasicSymbolic{Number}}, vars::Vector{SymbolicUtils.BasicSymbolic{SymbolicUtils.FnType{Tuple{Number}, Number, Nothing}}}, t_ind::SymbolicUtils.BasicSymbolic{Number}, c::Vector{Int}, d::Vector{Int}, Σ::Matrix{Float64}) -> Matrix{SymbolicUtils.BasicSymbolic{Number}}

Constructs the System Jacobian matrix J for the system of equations `eqs` with respect to the variables `vars`.
The offsets `c` and `d` and the signature matrix `Σ` are used to determine the structure of the Jacobian.
"""
function system_jacobian(
    eqs::Vector{SymbolicUtils.BasicSymbolic{Number}},
    vars::Vector{SymbolicUtils.BasicSymbolic{SymbolicUtils.FnType{Tuple{Number}, Number, Nothing}}},
    t_ind::SymbolicUtils.BasicSymbolic{Number},
    c::Vector{Int},
    d::Vector{Int},
    Σ::Matrix{Float64}
)
    n = length(eqs)
    J = zeros(Symbolics.Num, n, n)

    for i in 1:n
        for j in 1:n
            if d[j] - c[i] == Σ[i, j]
                f_i = eqs[i]
                x_j = vars[j]
                σ_ij = Int(Σ[i, j])
                # Compute the σ_ij-th derivative of x_j
                x_j_deriv = x_j(t_ind)
                for _ in 1:σ_ij
                    x_j_deriv = Differential(t_ind)(x_j_deriv)
                end
                # Compute the partial derivative ∂f_i / ∂x_j^(σ_ij)
                J[i, j] = expand_derivatives(Differential(x_j_deriv)(f_i))
            else
                # Set J[i, j] = 0 if d[j] - c[i] != Σ[i, j]
                J[i, j] = 0
            end
        end
    end

    return J
end