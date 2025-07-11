# Case study from: https://ieeexplore.ieee.org/document/7180413
# Case study from: https://ieeexplore.ieee.org/abstract/document/25627

using JuMP, CPLEX

# ----------------------------
# 1. Define Sets and Data
# ----------------------------

buses = 1:33

# Branch connectivity (from, to)
branches = [
    (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9), (9,10),
    (10,11), (11,12), (12,13), (13,14), (14,15), (15,16), (16,17),
    (17,18), (2,19), (19,20), (20,21), (21,22), (3,23), (23,24),
    (24,25), (6,26), (26,27), (27,28), (28,29), (29,30), (30,31),
    (31,32), (32,33)
]

# Define tie switches (operable branches)
tie_switches = [(9,15), (12,22), (18,33), (25,29), (32,20)]

# Define R (resistance) values for branches (pu) -- here just dummy values for demo
# Resistance (R) in p.u.
R = Dict(
    (1,2)=>0.0922, (2,3)=>0.4930, (3,4)=>0.3660, (4,5)=>0.3811, (5,6)=>0.8190,
    (6,7)=>0.1872, (7,8)=>0.7114, (8,9)=>1.0300, (9,10)=>1.0440, (10,11)=>0.1966,
    (11,12)=>0.3744, (12,13)=>1.4680, (13,14)=>0.5416, (14,15)=>0.5910, (15,16)=>0.7463,
    (16,17)=>1.2890, (17,18)=>0.7320, (2,19)=>0.1640, (19,20)=>1.5040, (20,21)=>0.4095,
    (21,22)=>0.7089, (3,23)=>0.4512, (23,24)=>0.8980, (24,25)=>0.8960, (6,26)=>0.2030,
    (26,27)=>0.2842, (27,28)=>1.0590, (28,29)=>0.8042, (29,30)=>0.5075, (30,31)=>0.9744,
    (31,32)=>0.3105, (32,33)=>0.3410
)

Pmax = Dict((i,j) => 1000.0 for (i,j) in branches)  # Simplified: assume 1000kW per branch

# 5. Define Load Data (example placeholders)
PL = Dict(i => 100.0 for i in buses)  # Each bus demands 100kW
QL = Dict(i => 50.0 for i in buses)   # Each bus demands 50kVAR

# Define maximum DG installation at each bus -- here assume only at selected buses
DG_max = Dict(i => (i in [6, 18, 33] ? 500.0 : 0.0) for i in buses)

# ----------------------------
# 2. Build JuMP Model
# ----------------------------

model=Model(CPLEX.Optimizer)
# MOI.set(model, MOI.RawOptimizerAttribute("print_level"), 1)

# DG power generation decision variable
@variable(model, 0 <= P_DG[i in buses] <= DG_max[i])

# Power flow on each branch
@variable(model, P[(i,j) in branches])

# Line switch decision variables only for tie switches
@variable(model, x[(i,j) in tie_switches], Bin)

# Fix all non-tie-switch branches as closed (x=1)
for (i,j) in branches
    if (i,j) ∉ tie_switches
        println("i,j= ",i, " ", j)
        println(Pmax[(i,j)])
        @constraint(model, P[(i,j)] <= Pmax[(i,j)])
        @constraint(model, P[(i,j)] >= -Pmax[(i,j)])
    else
        # For tie switches, tie flow to switch status
        @constraint(model, -Pmax[(i,j)] * x[(i,j)] <= P[i,j])
        @constraint(model, P[i,j] <= Pmax[(i,j)] * x[(i,j)])
    end
end
# Power balance at each bus
for i in buses
    incoming = [(k,i) for (k,j) in branches if j == i]
    outgoing = [(i,k) for (i,k) in branches if i == i]
    @constraint(model,
        sum(P[(k,i)] for (k,i) in incoming) - sum(P[(i,k)] for (i,k) in outgoing)
        == PL[i] - P_DG[i]
    )
end

# # Line switching logic
# for (i,j) in branches
#     if (i,j) ∉ tie_switches
#         @constraint(model, x[(i,j)] == 1)  # Non-tie lines must stay closed
#     end
# end

# Radiality constraint: total number of closed branches = n - 1
@constraint(model, sum(x[(i,j)] for (i,j) in tie_switches) + (length(branches) - length(tie_switches)) == length(buses) - 1)

# ----------------------------
# 4. Objective
# ----------------------------

# Minimize total line losses
@objective(model, Min, sum(R[(i,j)] * P[(i,j)]^2 for (i,j) in branches))

# ----------------------------
# 5. Solve
# ----------------------------

optimize!(model)

# ----------------------------
# 6. Output results
# ----------------------------
println("Objective value (Total Loss): ", objective_value(model))
for (i,j) in tie_switches
    println("Branch ($i, $j): Switch status = ", value(x[(i,j)]))
end


