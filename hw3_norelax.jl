using JuMP
using CPLEX

x_hat = 0
# Example data for W, T, S, and h
W = [1 -1 -1 -1 0 0;
     0 1 0 0 1 0;
     0 0 1 0 0 1]
# q = [1, 0, 0, 0, 0, 0]
q = [1.5, 0, 2/7, 1, 0, 0]   
# h = [-1, 2, 7] # kesi 1
h = [0, 2, 7] # kesi 2

# Define primal model
primal_model = Model(CPLEX.Optimizer)

# Define variables x, t, s
@variable(primal_model, x[1:6] >= 0)

@objective(primal_model, Min, q' * x)


JuMP.@constraints primal_model begin
   c1, x[1] - x[2] - x[3] -x[4] == h[1] - x_hat
   c2, x[2] + x[5] == 2
   c3, x[3] + x[6] == 7
end

# Solve the primal problem
optimize!(primal_model)

dual_vars1 = dual(c1)
dual_vars2 = dual(c2)
dual_vars3 = dual(c3)

# Retrieve primal optimal solution and value
primal_optimal_value = objective_value(primal_model)
primal_x = value.(x)

println("Primal optimal value: ", primal_optimal_value)
println("Primal optimal solution x: ", primal_x)

# Define the dual model

dual_model = Model(HiGHS.Optimizer)

# Define dual variables for constraints
@variable(dual_model, y[1:3] >= 0)

# Dual objective function: max h^T y
@objective(dual_model, Max, h' * y)

# Dual constraints: W^T y <= q
for j in 1:6
    @constraint(dual_model, W[:, j]' * y <= q[j])
end

# Solve the dual problem
optimize!(dual_model)

# Retrieve dual optimal value and dual solution
dual_optimal_value = objective_value(dual_model)
dual_y = value.(y)

println("Dual optimal value: ", dual_optimal_value)
println("Dual optimal solution y: ", dual_y)


# new problem
p = Model(CPLEX.Optimizer)

# Define variables x, t, s
@variable(p, x, lower_bound=-20, upper_bound=20)
@variable(p, t )
JuMP.@constraints p begin
    a1, t >= -5/4*x-1/2
    a2, t >=1/2*(x-7)
    a3, t >= 0
 end

@objective(p, Min, t)




# Solve the primal problem
optimize!(p)


# Retrieve primal optimal solution and value
p_optimal_value = objective_value(p)
p = value.(x)