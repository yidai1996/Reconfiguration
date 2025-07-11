using JuMP
using CPLEX

function bus_system_reconfiguration()
    println("new simulation starts")
# -----------------------------
# 1. Define sets and structure
# -----------------------------
buses = 1:33 # number of buses
T = 1:24  # 24-hour time periods

branches = [
    (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9), (9,10),
    (10,11), (11,12), (12,13), (13,14), (14,15), (15,16), (16,17),
    (17,18), (2,19), (19,20), (20,21), (21,22), (3,23), (23,24),
    (24,25), (6,26), (26,27), (27,28), (28,29), (29,30), (30,31),
    (31,32), (32,33),
    (8,21), (9,15), (12,22), (18,33), (25,29)
]

tie_switches = [(8,21), (9,15), (12,22), (18,33), (25,29)]
optional_main_lines = [(14,15), (21,22), (16,17)]
fixed_branches = setdiff(branches, tie_switches ∪ optional_main_lines)

# -----------------------------
# 2. Parameters and DG info
# -----------------------------
R = Dict(
    (1,2)=>0.0922, (2,3)=>0.4930, (3,4)=>0.3660, (4,5)=>0.3811, (5,6)=>0.8190,
    (6,7)=>0.1872, (7,8)=>0.7114, (8,9)=>1.0300, (9,10)=>1.0440, (10,11)=>0.1966,
    (11,12)=>0.3744, (12,13)=>1.4680, (13,14)=>0.5416, (14,15)=>0.5910, (15,16)=>0.7463,
    (16,17)=>1.2890, (17,18)=>0.7320, (2,19)=>0.1640, (19,20)=>1.5040, (20,21)=>0.4095,
    (21,22)=>0.7089, (3,23)=>0.4512, (23,24)=>0.8980, (24,25)=>0.8960, (6,26)=>0.2030,
    (26,27)=>0.2842, (27,28)=>1.0590, (28,29)=>0.8042, (29,30)=>0.5075, (30,31)=>0.9744,
    (31,32)=>0.3105, (32,33)=>0.3410,
    (8,21)=>2.0000, (9,15)=>2.0000, (12,22)=>2.0000, (18,33)=>0.50000, (25,29)=>0.5000
)
X = Dict(
    (1,2)=>0.0470, (2,3)=>0.2511, (3,4)=>0.1864, (4,5)=>0.1941, (5,6)=>0.7070,
    (6,7)=>0.6188, (7,8)=>0.2351, (8,9)=>0.7400, (9,10)=>0.7400, (10,11)=>0.0650,
    (11,12)=>0.1238, (12,13)=>1.1550, (13,14)=>0.7129, (14,15)=>0.5260, (15,16)=>0.5450,
    (16,17)=>1.7210, (17,18)=>0.5740, (2,19)=>0.1565, (19,20)=>1.3554, (20,21)=>0.4784,
    (21,22)=>0.9373, (3,23)=>0.3083, (23,24)=>0.7091, (24,25)=>0.7011, (6,26)=>0.1034,
    (26,27)=>0.1447, (27,28)=>0.9337, (28,29)=>0.7006, (29,30)=>0.2585, (30,31)=>0.9630,
    (31,32)=>0.3619, (32,33)=>0.5302,
    (8,21)=>2.0000, (9,15)=>2.0000, (12,22)=>2.0000, (18,33)=>0.50000, (25,29)=>0.5000
)

G = Dict(
    (1,2)=>8.6089, (2,3)=>1.6106, (3,4)=>2.1695, (4,5)=>2.0835, (5,6)=>0.6996,
    (6,7)=>0.4479, (7,8)=>1.2673, (8,9)=>0.6403, (9,10)=>0.6375, (10,11)=>4.5853,
    (11,12)=>2.4077, (12,13)=>0.4207, (13,14)=>0.6757, (14,15)=>0.9442,
    (15,16)=>0.8739, (16,17)=>0.2788, (17,18)=>0.8459, (2,19)=>3.1914,
    (19,20)=>0.3669, (20,21)=>1.0326, (21,22)=>0.5133, (3,23)=>1.5109,
    (23,24)=>0.6859, (24,25)=>0.6922, (6,26)=>3.9113, (26,27)=>2.7943,
    (27,28)=>0.5313, (28,29)=>0.7069, (29,30)=>1.5645, (30,31)=>0.5131,
    (31,32)=>1.3655, (32,33)=>0.8581, (8,21)=>0.25, (9,15)=>0.25,
    (12,22)=>0.25, (18,33)=>1.0, (25,29)=>1.0
)

B = Dict(
    (1,2)=>-4.3885, (2,3)=>-0.8203, (3,4)=>-1.1049, (4,5)=>-1.0612, (5,6)=>-0.604,
    (6,7)=>-1.4805, (7,8)=>-0.4188, (8,9)=>-0.4601, (9,10)=>-0.4519, (10,11)=>-1.516,
    (11,12)=>-0.7961, (12,13)=>-0.331, (13,14)=>-0.8894, (14,15)=>-0.8403,
    (15,16)=>-0.6382, (16,17)=>-0.3722, (17,18)=>-0.6634, (2,19)=>-3.0454,
    (19,20)=>-0.3307, (20,21)=>-1.2064, (21,22)=>-0.6787, (3,23)=>-1.0324,
    (23,24)=>-0.5416, (24,25)=>-0.5417, (6,26)=>-1.9923, (26,27)=>-1.4227,
    (27,28)=>-0.4684, (28,29)=>-0.6159, (29,30)=>-0.7969, (30,31)=>-0.5131,
    (31,32)=>-1.5916, (32,33)=>-1.3342, (8,21)=>-0.25, (9,15)=>-0.25,
    (12,22)=>-0.25, (18,33)=>-1.0, (25,29)=>-1.0
)


Pmax = Dict((i,j) => 1000.0 for (i,j) in branches)
PL_base = Dict(
    2=>100, 3=>90, 4=>120, 5=>60, 6=>60,
    7=>200, 8=>200, 9=>60, 10=>60, 11=>45,
    12=>60, 13=>60, 14=>120, 15=>60, 16=>60,
    17=>60, 18=>90, 19=>90, 20=>90, 21=>90,
    22=>90, 23=>90, 24=>420, 25=>420, 26=>60,
    27=>60, 28=>60, 29=>120, 30=>200, 31=>150,
    32=>210, 33=>60
)
PL = Dict((i, t) => get(PL_base, i, 0.0) for i in buses for t in T)
QL_base = Dict(
    2=>60, 3=>40, 4=>80, 5=>30, 6=>20,
    7=>100, 8=>100, 9=>20, 10=>20, 11=>30,
    12=>35, 13=>35, 14=>80, 15=>10, 16=>20,
    17=>20, 18=>40, 19=>40, 20=>40, 21=>40,
    22=>40, 23=>50, 24=>200, 25=>200, 26=>25,
    27=>25, 28=>20, 29=>70, 30=>600, 31=>70,
    32=>100, 33=>40
)
QL = Dict((i, t) => get(QL_base, i, 0.0) for i in buses for t in T)

DG_max = Dict(i => 0.0 for i in buses)
DG_max[10] = 400.0  # wind
DG_max[33] = 500.0  # wind
DG_max[7]  = 350.0  # solar
DG_max[14] = 450.0  # solar

forecast_DG = Dict{Tuple{Int, Int}, Float64}()
for t in T
    forecast_DG[(10, t)] = 400 # wind pattern
    forecast_DG[(33, t)] = 500 
    forecast_DG[(7, t)]  = 350  # solar pattern
    forecast_DG[(14, t)] = 450 
end

V_sq_data = Dict(
    2 => 0.9927, 3=>0.9574, 4=>0.9374, 5=>0.9176, 6=>0.8707, 7=>0.8641, 
    8=>0.8550, 9=>0.8432, 10=>0.8324, 11=>0.8308, 12=>0.8280, 13=>0.8167, 
    14=>0.8125, 15=>0.8099, 16=>0.8074, 17=>0.8037, 18=>0.8026, 
    19=>0.9916, 20=>0.9845, 21=>0.9831, 22=>0.9818,
    23=>0.9504, 24=>0.9373, 25=>0.9309,
    26=>0.8643, 27=>0.8557, 28=>0.8201, 29=>0.7945, 30=>0.7816, 31=>0.7739, 32=>0.7723, 33=>0.7717
)

# -----------------------------
# 3. Model
# -----------------------------
model = Model(CPLEX.Optimizer)

@variable(model, x_tie[i in tie_switches], Bin)
@variable(model, x_opt[i in optional_main_lines], Bin)

x = Dict{Tuple{Int,Int}, Any}()
for l in tie_switches x[l] = x_tie[l] end
for l in optional_main_lines x[l] = x_opt[l] end
for l in fixed_branches x[l] = 1 end

@variable(model, P[i in branches, t in T])
@variable(model, Q[i in branches, t in T])
@variable(model, 0 <= P_DG[i in buses, t in T] <= DG_max[i])
@variable(model, V_sq[i in buses, t in T], lower_bound = 0.7000, upper_bound = 1.1025)
@variable(model, scm[i in branches, t in T] >= 0)

# -----------------------------
# 4. Constraints
# -----------------------------
for t in T
    for i in buses
        inflow = [(k,i) for (k,j) in branches if j == i]
        outflow = [(i,k) for (j,k) in branches if j == i]

        @constraint(model,
            sum(P[(k,i),t]-R[(k,i)] * scm[(k,i),t] for (k,i) in inflow) - sum(P[(i,k),t] for (i,k) in outflow)
            == PL[(i, t)] - P_DG[i,t])

        @constraint(model,
            sum(Q[(k,i),t]-X[(k,i)]*scm[(k,i),t] for (k,i) in inflow) - sum(Q[(i,k),t] for (i,k) in outflow)
            == QL[(i, t)])
    end

    for (i,j) in tie_switches ∪ optional_main_lines
        @constraint(model, -Pmax[(i,j)] * x[(i,j)] <= P[(i,j),t] )
        @constraint(model, P[(i,j),t] <= Pmax[(i,j)] * x[(i,j)])
        @constraint(model, -Pmax[(i,j)] * x[(i,j)] <= Q[(i,j),t] )
        @constraint(model, Q[(i,j),t] <= Pmax[(i,j)] * x[(i,j)])
    end

    for (i,j) in branches
        @constraint(model, V_sq[j,t] <= V_sq[i,t] - 2 * (R[(i,j)] * P[(i,j),t] + X[(i,j)] * Q[(i,j),t]) + (R[(i,j)]^2+X[(i,j)]^2)*scm[(i,j),t])
        @constraint(model, [2*P[(i,j),t], 2*Q[(i,j),t], V_sq[i,t] - scm[(i,j),t]] in SecondOrderCone())
        # @constraint(model, [sqrt(V_sq[i,t]), P[(i,j),t], Q[(i,j),t], sqrt(scm[(i,j),t])] in SecondOrderCone())
    end

    for i in [10, 33, 7, 14]  # only DG buses
        @constraint(model, P_DG[i,t] <= forecast_DG[(i,t)])
    end
end

for i in keys(V_sq_data), t in T
    set_start_value(V_sq[i, t], V_sq_data[i])
end

@constraint(model,
    sum(x_tie[(i,j)] for (i,j) in tie_switches) +
    sum(x_opt[(i,j)] for (i,j) in optional_main_lines) +
    length(fixed_branches) == 32)

# -----------------------------
# 5. Objective function
# -----------------------------
@objective(model, Min, sum(R[(i,j)] * scm[(i,j),t] for (i,j) in branches, t in T))

if termination_status(model) == MOI.OPTIMAL
    println("Objective value: ", objective_value(model))
else
    println("No solution found. Termination: ", termination_status(model))
end

# -----------------------------
# 6. Solve and print
# -----------------------------
optimize!(model)

end

# println("Objective value (Total loss): ", objective_value(model))
# for (i,j) in tie_switches ∪ optional_main_lines
#     println("Switch ($i,$j): status = ", value(x[(i,j)]))
# end
# for t in T, i in buses
#     if get(DG_max, i, 0.0) > 0
#         println("Bus $i at hour $t: DG output = ", value(P_DG[i,t]))
#     end
# end

bus_system_reconfiguration()