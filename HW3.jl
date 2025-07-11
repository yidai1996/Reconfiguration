using JuMP
using HiGHS

# Define the primal problem's parameters
# c = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # Objective coefficients of the primal
W = [1 -1 -1 -1 0 0;
     0 1 0 0 1 0;
     0 0 1 0 0 1;
     -1 0 0 0 0 0;
     0 -1 0 0 0 0;
     0 0 -1 0 0 0;
     0 0 0 -1 0 0;
     0 0 0 0 -1 0;
     0 0 0 0 0 -1]     # Coefficient matrix of the primal
h = [1, 2, 7, 0, 0, 0, 0, 0, 0]       # Right-hand side of the primal
T = [-1 0 0;
    0 -1 0;
    0 0 -1;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0]
S = [0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     -1 0 0 0 0 0;
     0 -1 0 0 0 0;
     0 0 -1 0 0 0;
     0 0 0 -1 0 0;
     0 0 0 0 -1 0;
     0 0 0 0 0 -1]


# Define primal variables count
n_x = 6  # Length of x
n_t = 3  # Length of t
n_s = 6  # Length of s
n_constraints = 9  # Number of constraints

# Create the dual model
dual_model = Model(HiGHS.Optimizer)

# Define dual variables for constraints
@variable(dual_model, y[1:n_constraints] >= 0)

# Dual objective function
@objective(dual_model, Max, h' * y)

# Dual constraint related to primal t-variables (sum(t_i))
for i in 1:n_t
    @constraint(dual_model, T[:, i]' * y <= 1)
end

# Dual constraint related to primal s-variables (sum(s_j))
for j in 1:n_s
    @constraint(dual_model, S[:, j]' * y <= 1)
end

# Dual constraint related to primal x-variables
for k in 1:n_x
    @constraint(dual_model, W[:, k]' * y <= 0)
end

# Solve the dual problem
optimize!(dual_model)

# Retrieve dual optimal value and solution
dual_optimal_value = objective_value(dual_model)
dual_optimal_solution = value.(y)

println("Dual optimal value: ", dual_optimal_value)
println("Dual optimal solution: ", dual_optimal_solution)