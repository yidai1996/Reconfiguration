# ========================================================================
# Reconfiguration with renewable generation
# ========================================================================

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
tie_switches = [(9,15), (12,22), (18,33), (25,29), (8,21)]

# Define R (resistance) values for branches (pu) -- here just dummy values for demo
# Resistance (R) in p.u.
# r (ohm)
r = Dict(
    (1,2)=>0.0922, (2,3)=>0.4930, (3,4)=>0.3660, (4,5)=>0.3811, (5,6)=>0.8190,
    (6,7)=>0.1872, (7,8)=>0.7114, (8,9)=>1.0300, (9,10)=>1.0440, (10,11)=>0.1966,
    (11,12)=>0.3744, (12,13)=>1.4680, (13,14)=>0.5416, (14,15)=>0.5910, (15,16)=>0.7463,
    (16,17)=>1.2890, (17,18)=>0.7320, (2,19)=>0.1640, (19,20)=>1.5040, (20,21)=>0.4095,
    (21,22)=>0.7089, (3,23)=>0.4512, (23,24)=>0.8980, (24,25)=>0.8960, (6,26)=>0.2030,
    (26,27)=>0.2842, (27,28)=>1.0590, (28,29)=>0.8042, (29,30)=>0.5075, (30,31)=>0.9744,
    (31,32)=>0.3105, (32,33)=>0.3410,
    # tie lines
    (8,21)=>2.0000, (9,15)=>2.0000, (12,22)=>2.0000, (18,33)=>0.50000, (25,29)=>0.5000
)

# x(ohm) 
x = Dict(
    (1,2)=>0.0470, (2,3)=>0.2511, (3,4)=>0.1864, (4,5)=>0.1941, (5,6)=>0.7070,
    (6,7)=>0.6188, (7,8)=>0.2351, (8,9)=>0.7400, (9,10)=>0.7400, (10,11)=>0.0650,
    (11,12)=>0.1238, (12,13)=>1.1550, (13,14)=>0.7129, (14,15)=>0.5260, (15,16)=>0.5450,
    (16,17)=>1.7210, (17,18)=>0.5740, (2,19)=>0.1565, (19,20)=>1.3554, (20,21)=>0.4784,
    (21,22)=>0.9373, (3,23)=>0.3083, (23,24)=>0.7091, (24,25)=>0.7011, (6,26)=>0.1034,
    (26,27)=>0.1447, (27,28)=>0.9337, (28,29)=>0.7006, (29,30)=>0.2585, (30,31)=>0.9630,
    (31,32)=>0.3619, (32,33)=>0.5302,
    # tie lines
    (8,21)=>2.0000, (9,15)=>2.0000, (12,22)=>2.0000, (18,33)=>0.50000, (25,29)=>0.5000
)


# Define Load Data (example placeholders)
# PL(kW)
PL = Dict(
    2=>100, 3=>90, 4=>120, 5=>60, 6=>60,
    7=>200, 8=>200, 9=>60, 10=>60, 11=>45,
    12=>60, 13=>60, 14=>120, 15=>60, 16=>60,
    17=>60, 18=>90, 19=>90, 20=>90, 21=>90,
    22=>90, 23=>90, 24=>420, 25=>420, 26=>60,
    27=>60, 28=>60, 29=>120, 30=>200, 31=>150,
    32=>210, 33=>60
)
# QL(kvar)
QL = Dict(
    2=>60, 3=>40, 4=>80, 5=>30, 6=>20,
    7=>100, 8=>100, 9=>20, 10=>20, 11=>30,
    12=>35, 13=>35, 14=>80, 15=>10, 16=>20,
    17=>20, 18=>40, 19=>40, 20=>40, 21=>40,
    22=>40, 23=>50, 24=>200, 25=>200, 26=>25,
    27=>25, 28=>20, 29=>70, 30=>600, 31=>70,
    32=>100, 33=>40
)

Pmax = Dict((i,j) => 1000.0 for (i,j) in branches) 

# Define maximum DG installation at each bus -- here assume only at selected buses
DG_max = Dict(i => (i in [6, 18, 33] ? 500.0 : 0.0) for i in buses)

# Squared voltage magnitude at bus i 
V_sq = Dict(
    2 => 0.9927, 3=>0.9574, 4=>0.9374, 5=>0.9176, 6=>0.8707, 7=>0.8641, 
    8=>0.8550, 9=>0.8432, 10=>0.8324, 11=>0.8308, 12=>0.8280, 13=>0.8167, 
    14=>0.8125, 15=>0.8099, 16=>0.8074, 17=>0.8037, 18=>0.8026, 
    19=>0.9916, 20=>0.9845, 21=>0.9831, 22=>0.9818,
    23=>0.9504, 24=>0.9373, 25=>0.9309,
    26=>0.8643, 27=>0.8557, 28=>0.8201, 29=>0.7945, 30=>0.7816, 31=>0.7739, 32=>0.7723, 33=>0.7717
)

# ----------------------------
# 2. Build JuMP Model
# ----------------------------

model=Model(CPLEX.Optimizer)
# MOI.set(model, MOI.RawOptimizerAttribute("print_level"), 1)

# DG power generation decision variable
@variable(model, 0 <= PG[i in buses] <= DG_max[i])

# Power flow on each branch
@variable(model, P[(i,j) in branches])

@variable(model, alpha[(i,j) in branches], Bin)

@constraints MPC begin
    # Active and reactive power balance at all buses 
    pb1[i], PG[i] - PD[i] == sum(P[i,j]) - sum(P[k,i]-r[k,i]*l[k,i] for (k,i) in E) + g[i]*v[i]
    pb2[i], QG[i] - QD[i] == sum(Q[i,j]) - sum(Q[k,i]-x[k,i]*l[k,i] for (k,i) in E) + b[i]*v[i]

    # voltage magnitude 
    vm1[i,j], v[j] <= M(1-a[i,j]) + v[i] - 2*(r[i,j]*P[i,j] + x[i,j]*Q[i,j]) + (r[i,j]^2 + x[i,j]^2)*l[i,j]
    vm2[i,j], v[j] >= - M(1-a[i,j]) + v[i] - 2*(r[i,j]*P[i,j] + x[i,j]*Q[i,j]) + (r[i,j]^2 + x[i,j]^2)*l[i,j]

    # radiality and connectivity of DNs 
    c14, sum(a[i,j]) == n - 1
    c15[i,j], b[i,j] + b[j,i] == a[i,j]
    c16, sum(b[i,j] for i in N for i not in mg for j in N_i) == 1

    # Power flow equations, cone constraints
    pf[i,j], l[i,j] >= (P[i,j]^2 + Q[i,j]^2)/v[i] # forall (i,j) \in E
    
    # limits of l[i,j]
    limit_of_l[i,j], l[i,j] <= a[i,j]*I_bar[i,j]^2
    # limits of v[i]
    limit_of_v1[i in N and not in mg], v[i]<= V_upper_bar[i]^2
    limit_of_v2[i in N and not in mg], v[i]>= V_lower_bar[i]^2

    # voltage equal to 1 pu 
    v_equal, v__mg == 1

    # Include PV mode in the reconfiguration problem
    c22[i in N(PV)], gama_lo[i] + gama_up[i] <= 1
    c23[i in N(PV)], v_i <= gama_up[i]*V[0,i]^2 + (1-gama_up[i])*V_upper_bar[i]^2
    c24[i in N(PV)], v_i >= gama_lo[i]*V[0,i]^2 + (1-gama_lo[i])*V_lower_bar[i]^2
    c25[i in N(PV)], v_i <= (1-gama_up[i]-gama_lo[i])*V[0,i]^2 + (gama_up[i]+gama_lo[i])*V_upper_bar[i]^2
    c26[i in N(PV)], v_i >= (1-gama_up[i]-gama_lo[i])*V[0,i]^2 + (gama_up[i]+gama_lo[i])*V_lower_bar[i]^2
    c27[i in N(PV)], QG[i] <= QG_upper[i] + gama_lo[i]*(QG_lower[i]-QG_upper[i])
    c28[i in N(PV)], QG[i] >= QG_lower[i] + gama_up[i]*(QG_upper[i]-QG_lower[i])
end


@objective(MPC,Min,sum(r[i,j]*l[i,j] for (i,j) in E))

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

# Adjust radiality constraint: still n-1 branches active
@constraint(model, sum(x[i,j] for (i,j) in branches) == length(buses) - 1)



# Reactance (X) in p.u.
X = Dict(
    (1,2)=>0.0470, (2,3)=>0.2511, (3,4)=>0.1864, (4,5)=>0.1941, (5,6)=>0.7070,
    (6,7)=>0.6188, (7,8)=>0.2351, (8,9)=>0.7400, (9,10)=>0.7400, (10,11)=>0.0650,
    (11,12)=>0.1238, (12,13)=>1.1550, (13,14)=>0.7129, (14,15)=>0.5260, (15,16)=>0.5450,
    (16,17)=>1.7210, (17,18)=>0.5740, (2,19)=>0.1565, (19,20)=>1.3550, (20,21)=>0.4784,
    (21,22)=>0.9373, (3,23)=>0.3083, (23,24)=>0.7091, (24,25)=>0.7011, (6,26)=>0.1034,
    (26,27)=>0.1447, (27,28)=>0.9337, (28,29)=>0.7006, (29,30)=>0.2585, (30,31)=>0.9630,
    (31,32)=>0.3619, (32,33)=>0.5302
)

@objective(model, Min, sum(R[(i,j)] * P[i,j]^2 for (i,j) in branches))


# below is the formulation from 10.1109/TPWRS.2015.2457954

function loadProcessData()
    # b = # susceptance of bus i 
    # E = # set of lines
end

function Problem_Formulation()

    MPC=Model(CPLEX.Optimizer)
    MOI.set(MPC, MOI.RawOptimizerAttribute("print_level"), 1)

    JuMP.@variables MPC begin
       
    end

    @constraints MPC begin
        # Active and reactive power balance at all buses 
        pb1[i], PG[i] - PD[i] == sum(P[i,j]) - sum(P[k,i]-r[k,i]*l[k,i] for (k,i) in E) + g[i]*v[i]
        pb2[i], QG[i] - QD[i] == sum(Q[i,j]) - sum(Q[k,i]-x[k,i]*l[k,i] for (k,i) in E) + b[i]*v[i]

        # voltage magnitude 
        vm1[i,j], v[j] <= M(1-a[i,j]) + v[i] - 2*(r[i,j]*P[i,j] + x[i,j]*Q[i,j]) + (r[i,j]^2 + x[i,j]^2)*l[i,j]
        vm2[i,j], v[j] >= - M(1-a[i,j]) + v[i] - 2*(r[i,j]*P[i,j] + x[i,j]*Q[i,j]) + (r[i,j]^2 + x[i,j]^2)*l[i,j]

        # radiality and connectivity of DNs 
        c14, sum(a[i,j]) == n - 1
        c15[i,j], b[i,j] + b[j,i] == a[i,j]
        c16, sum(b[i,j] for i in N for i not in mg for j in N_i) == 1

        # Power flow equations, cone constraints
        pf[i,j], l[i,j] >= (P[i,j]^2 + Q[i,j]^2)/v[i] # forall (i,j) \in E
        
        # limits of l[i,j]
        limit_of_l[i,j], l[i,j] <= a[i,j]*I_bar[i,j]^2
        # limits of v[i]
        limit_of_v1[i in N and not in mg], v[i]<= V_upper_bar[i]^2
        limit_of_v2[i in N and not in mg], v[i]>= V_lower_bar[i]^2

        # voltage equal to 1 pu 
        v_equal, v__mg == 1

        # Include PV mode in the reconfiguration problem
        c22[i in N(PV)], gama_lo[i] + gama_up[i] <= 1
        c23[i in N(PV)], v_i <= gama_up[i]*V[0,i]^2 + (1-gama_up[i])*V_upper_bar[i]^2
        c24[i in N(PV)], v_i >= gama_lo[i]*V[0,i]^2 + (1-gama_lo[i])*V_lower_bar[i]^2
        c25[i in N(PV)], v_i <= (1-gama_up[i]-gama_lo[i])*V[0,i]^2 + (gama_up[i]+gama_lo[i])*V_upper_bar[i]^2
        c26[i in N(PV)], v_i >= (1-gama_up[i]-gama_lo[i])*V[0,i]^2 + (gama_up[i]+gama_lo[i])*V_lower_bar[i]^2
        c27[i in N(PV)], QG[i] <= QG_upper[i] + gama_lo[i]*(QG_lower[i]-QG_upper[i])
        c28[i in N(PV)], QG[i] >= QG_lower[i] + gama_up[i]*(QG_upper[i]-QG_lower[i])
    end


    @objective(MPC,Min,sum(r[i,j]*l[i,j] for (i,j) in E))

    JuMP.optimize!(MPC)

    
    if st=="Infeasible_Problem_Detected"
        println()
        error("Solver infeasible, problem stopping")
    end
    # The error function throws a generic ErrorException. This will interrupt execution of the function or block immediately.
    # receive a solver specific string explaning why the optimization stopped
    obj=objective_value(MPC)
    return a, obj

end

