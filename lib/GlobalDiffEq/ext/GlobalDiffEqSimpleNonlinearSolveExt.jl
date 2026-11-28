module GlobalDiffEqSimpleNonlinearSolveExt

using SimpleNonlinearSolve: SimpleNewtonRaphson

# The solver of the algebraic projection in the DAE saltation products
# (`GlobalDiffEq._algebraic_solver`).
algebraic_solver() = SimpleNewtonRaphson()

end
