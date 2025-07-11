# reconfigure_33bus.jl

using JuMP
# choose your solver here; Gurobi is recommended for speed
# using Gurobi
# model = Model(Gurobi.Optimizer)
using CPLEX
model = Model(CPLEX.Optimizer)

# -----------------------------
# 1. Define the 33‐bus system data
# -----------------------------
# Bus indices 1..33 (1 is the substation/root)
buses = 1:33

# (P_load, Q_load) at each bus in kW, kVAr
# data from Baran & Wu (1989) IEEE 33‐node test feeder
P_load = Dict{Int,Float64}(
  2=>100.0, 3=>90.0, 4=>120.0, 5=>60.0, 6=>60.0,
  7=>200.0, 8=>200.0, 9=>60.0, 10=>60.0, 11=>45.0,
  12=>60.0, 13=>60.0, 14=>120.0, 15=>200.0, 16=>150.0,
  17=>200.0, 18=>210.0, 19=>60.0, 20=>60.0, 21=>45.0,
  22=>120.0, 23=>90.0, 24=>60.0, 25=>60.0, 26=>120.0,
  27=>200.0, 28=>200.0, 29=>200.0, 30=>200.0, 31=>210.0,
  32=>60.0, 33=>60.0
)
Q_load = Dict(i => 0.38*P_load[i] for i in keys(P_load))  # typical PF=0.92

# Optionally: DG at certain buses (in kW)
P_gen_max = Dict{Int,Float64}(18=>250.0, 33=>250.0)  # example

# Branch data: (from, to, r [Ω], x [Ω], is_switchable)
raw_branches = [
  (1,2,0.0922,0.0470,false),  (2,3,0.4930,0.2511,false),
  (3,4,0.3660,0.1864,false),(4,5,0.3811,0.1941,false),
  (5,6,0.8190,0.7070,false),(6,7,0.1872,0.6188,false),
  (7,8,1.7114,1.2351,false),(8,9,1.0300,0.7400,false),
  (9,10,1.0410,0.7400,false),(10,11,0.1966,0.0650,false),
  (11,12,0.3744,0.1238,false),(12,13,1.4680,1.1550,false),
  (13,14,0.5416,0.7129,false),(14,15,0.5910,0.5260,false),
  (15,16,0.7463,0.5450,false),(16,17,1.2890,1.7210,false),
  (17,18,0.7320,0.5740,false),(2,19,0.1640,0.1565,false),
  (19,20,1.5042,1.3554,false),(20,21,0.4095,0.4784,false),
  (21,22,0.7089,0.9373,false),(3,23,0.4512,0.3083,false),
  (23,24,0.8980,0.7091,false),(24,25,0.8960,0.7011,false),
  (6,26,0.2030,0.1034,false),(26,27,0.2842,0.1447,false),
  (27,28,1.2043,0.6166,false),(28,29,0.5075,0.2585,false),
  (29,30,0.0976,0.0400,false),(30,31,0.7114,0.2351,false),
  (31,32,1.0300,0.7400,false),(32,33,0.1640,0.1565,false),
  # Example switchable tie‐lines (not in original): 
  (8, 14, 0.1, 0.1, true), (15,33, 0.1, 0.1, true)
]

# -----------------------------
# 2. Per‐unit conversion
# -----------------------------
S_base = 1000.0   # 1 MVA = 1000 kW
V_base = 12.66e3  # line‐to‐line in volts

branches = []
for (i,j,r,x,sw) in raw_branches
    r_pu = r * S_base / (V_base^2)
    x_pu = x * S_base / (V_base^2)
    push!(branches, (i, j, r_pu, x_pu, sw))
end

# convert loads and gens to PU
P_L = Dict(i => P_load[i]/S_base for i in keys(P_load))
Q_L = Dict(i => Q_load[i]/S_base for i in keys(Q_load))
P_Gmax = Dict(i => P_gen_max[i]/S_base for i in keys(P_gen_max))

# -----------------------------
# 3. Build the JuMP model
# -----------------------------
@variable(model, 0.9 <= V[buses] <= 1.1)    # bus voltage magnitudes (PU)
@constraint(model, V[1] == 1.0)             # slack/substation

nb = length(buses)
nl = length(branches)

# indexing branches 1..nl
@variable(model, Pij[1:nl])    # real power on branch
@variable(model, Qij[1:nl])    # reactive power on branch
@variable(model, Iij[1:nl] >= 0)  # squared current magnitude
@variable(model, z[1:nl], Bin)    # branch status

# big‐M for flow limits if open
M_P = 1.0
M_Q = 1.0
M_I = 1.1

# -----------------------------
# 4. Constraints
# -----------------------------
# 4.1 Radiality: exactly (nb-1) branches closed
@constraint(model, sum(z) == nb-1)

# Precompute adjacency for flow balance
children = Dict(i=>Int[] for i in buses)
for (k,(i,j,_,_,_)) in enumerate(branches)
    push!(children[i], k)
end

# 4.2 Power‐balance at each bus
for i in buses
    if i == 1
        # at substation, injection = sum flows out
        @constraint(model,
          sum(Pij[k] for k in children[i]) == 0 )
    else
        # flow in from parent(s) minus flow out = load − DG
        inflow = sum(Pij[k] for (k,(ii,jj,_,_,_)) in enumerate(branches) if jj==i)
        outflow = sum(Pij[k] for k in children[i])
        Pd = get(P_L,i,0.0); Pg = get(P_Gmax,i,0.0)
        @constraint(model, inflow - outflow == Pd - Pg)
    end
end

# 4.3 Reactive‐power balance (same structure)
for i in buses
    inflow = sum(Qij[k] for (k,(ii,jj,_,_,_)) in enumerate(branches) if jj==i)
    outflow = sum(Qij[k] for k in children[i])
    Qd = get(Q_L,i,0.0)
    @constraint(model, inflow - outflow == Qd)
end

# 4.4 Voltage drop along each branch
for (k,(i,j,r,x,sw)) in enumerate(branches)
    @constraint(model,
      V[j] == V[i]
             - 2*( r*Pij[k] + x*Qij[k] )
             + (r^2 + x^2)*Iij[k] )
end

# 4.5 Current‐flow relation: P^2 + Q^2 ≤ V_i * I
# This is a second‐order cone; approximate with piecewise or if supported use CP‐Solver:
for k in 1:nl
    @constraint(model, Pij[k]^2 + Qij[k]^2 <= V[branches[k][1]] * Iij[k])
end

# 4.6 On/off logic: flows zero if branch open
for k in 1:nl
    @constraint(model, Pij[k] <= M_P * z[k])
    @constraint(model, -Pij[k] <= M_P * z[k])
    @constraint(model, Qij[k] <= M_Q * z[k])
    @constraint(model, -Qij[k] <= M_Q * z[k])
    @constraint(model, Iij[k] <= M_I * z[k])
end

# -----------------------------
# 5. Objective: minimize losses
# -----------------------------
@objective(model, Min, sum(branches[k][3] * Iij[k] for k=1:nl))

# -----------------------------
# 6. Solve and extract results
# -----------------------------
optimize!(model)

println("\nOptimal radial configuration:")
for k in 1:nl
    if value(z[k]) > 0.5
        println("  Branch ", branches[k][1], "→", branches[k][2], " CLOSED")
    end
end

println("\nBus voltages (PU):")
for i in buses
    println("  V[$i] = ", round(value(V[i]), digits=4))
end

println("\nTotal losses (PU): ", objective_value(model))
