# https://arxiv.org/abs/2411.11791

# ieee37_reconfiguration.jl
# Julia script to model and solve a flexible network reconfiguration problem
# on a modified IEEE 37‐node feeder using PowerModelsONM.jl and CPLEX.

using JuMP
using CPLEX
using PowerModelsDistribution
using PowerModelsONM
using XLSX, DataFrames
# at the top of your script or REPL session:
empty!(Base.ARGS)
push!(Base.ARGS, "--network", "network.json", "--output", "result.json")
# 1. --- Read and assemble network data ---------------------------------------

# Base system values (placeholder guesses)
const BASE_MVA = 1.0          # guess: per-unit base power
const BASE_V_KV = 4.8         # guess: nominal feeder voltage in kV
const BASE_Z = BASE_V_KV^2 / BASE_MVA  # computed base impedance in ohm
const SUBSTATIONS = [799]

# 1.1 Load data: aggregate per-phase spot loads from Excel (data from Spot Loads.xls)
load_file = "Spot Loads.xlsx"
spot = XLSX.readtable(load_file, "Sheet1")
# convert to DataFrame
spot_df = DataFrame(spot)
# Create Dict of active/reactive loads in MW and MVAr (converted from kW/kVAr)
pd = Dict{Int,Float64}()
qd = Dict{Int,Float64}()
for row in eachrow(spot_df)
    bus = parse(Int, string(row[:Node]))
    # columns are from file: "Ph-1 kW", "Ph-1 kVAr", etc.
    p1 = row[:"Ph-1 kW"]
    q1 = row[:"Ph-1 kVAr"]
    p2 = row[:"Ph-2 kW"]
    q2 = row[:"Ph-2 kVAr"]
    p3 = row[:"Ph-3 kW"]
    q3 = row[:"Ph-3 kVAr"]
    # convert kW→MW, kVAr→MVAr
    pd[bus] = (p1 + p2 + p3) / 1000
    qd[bus] = (q1 + q2 + q3) / 1000
end

# 1.2 Line segment data: lengths and configuration codes (data from Line Data.xls)
line_df = DataFrame(XLSX.readtable("Line Data.xlsx", "Sheet1"))
# Self-phase impedance (R + jX) in ohm/mile for each config (from standard feeder data)
const CONFIG_Z = Dict(
    721 => (R=0.2926, X=0.1973),  # from manufacturer/IEEE doc
    722 => (R=0.4751, X=0.2973),
    723 => (R=1.2936, X=0.6713),
    724 => (R=2.0952, X=0.7758),
)

# 2. --- Build POWERMODELS data structure -------------------------------------
data = Dict{String,Any}()
data["bus"] = Dict{String, Dict{String,Any}}()
data["baseMVA"] = BASE_MVA
data["gen"] = Dict{String, Dict{String,Any}}()
data["branch"] = Dict{String,Any}()

# Assemble bus list from line data and spot loads
all_buses = unique(vcat(Int.(line_df[!,:NodeA]), Int.(line_df[!,:NodeB])))
for bus in keys(pd)
    push!(all_buses, bus)
end
# Substations (slack buses) from system document\const SUBSTATIONS = [798, 799]

# Populate bus entries
for bus in sort(unique(all_buses))
    pd_val = get(pd, bus, 0.0)  # from Spot Loads.xls or zero if none
    qd_val = get(qd, bus, 0.0)
    data["bus"][string(bus)] = Dict(
        "bus_i" => bus,
        "type"  => (bus in SUBSTATIONS ? 3 : 1),  # 3=slack (from doc), 1=PQ
        "pd"    => pd_val,
        "qd"    => qd_val,
    )
end

# 2.1 Generators: solar PV units (buses from IEEE doc, capacity guess)
const PV_BUSES = [741, 740, 735, 725, 724]  # from system description
const PV_CAP_MW = 0.1  # placeholder guess: 100 kW each
gen_id = 1
for bus in PV_BUSES
    data["gen"][string(gen_id)] = Dict(
        "gen_id" => gen_id,
        "bus"    => bus,
        "pg"     => 0.0,
        "qg"     => 0.0,
        "pmax"   => PV_CAP_MW,
        "pmin"   => 0.0,
        "qmax"   => 0.0,
        "qmin"   => 0.0,
        "status" => 1,
    )
    gen_id += 1
end

# 2.2 Branches: switchable vs non-switchable (from provided lists)
const SWITCHABLE = Set([
    (737,734), (734,710), (736,732), (708,733), (709,730),
    (731,718), (704,720), (703,702), (702,705), (702,701), (742,744)
])

branch_id = 1
for row in eachrow(line_df)
    f, t = Int(row[:NodeA]), Int(row[:NodeB])
    length_ft = Float64(row[:Length])  # from Line Data.xls
    cfg = Int(row[:Config])               # from Line Data.xls
    # Lookup impedance per-mile from CONFIG_Z (IEEE doc values)
    zcfg = CONFIG_Z[cfg]
    miles = length_ft / 5280.0
    # convert to per-unit using BASE_Z
    r_pu = zcfg.R * miles / BASE_Z
    x_pu = zcfg.X * miles / BASE_Z
    data["branch"][string(branch_id)] = Dict(
        "branch_i"  => branch_id,
        "f_bus"     => f,
        "t_bus"     => t,
        "r"         => r_pu,
        "x"         => x_pu,
        "b"         => 0.0,
        "rate_a"    => 100.0,            # guess: line MVA rating
        "status"    => 1,
        "switchable"=> ((f,t) in SWITCHABLE ? 1 : 0), # =1 means it's a switchable line
    )
    branch_id += 1
end


using JSON
open("network.json", "w") do io
  JSON.print(io, data)
end

open("onm_settings.json", "w") do io
    JSON.print(io, Dict(
      "switch_model"  => "LPUBFDiagPowerModel",
      "switch_only"   => true,
      "dispatch_only" => false
    ))
end
# 3. --- Configure and run PowerModelsONM --------------------------------------

# from https://github.com/lanl-ansi/PowerModelsONM.jl/blob/main/test/args.jl
    # @test args == Dict{String,Any}(
    #     "network" => "../test/data/ieee13_feeder.dss",
    #     "output" => "../test/data/test_output.json",
    #     "faults" => "../test/data/ieee13_faults.json",
    #     "inverters" => "../test/data/ieee13_inverters.json",
    #     "settings" => "../test/data/ieee13_settings.json",
    #     "events" => "../test/data/ieee13_events.json",
    #     "log-level" => "debug",
    #     "gurobi" => true,
    #     "opt-disp-formulation" => "acr",
    #     "opt-disp-solver" => "misocp_solver",
    #     "skip" => String["faults", "stability"],
    #     "pretty-print" => true,
    #     "disable-presolver" => true,
    #     "disable-isolation-constraint" => true,
    #     "disable-radial-constraint" => true,
    #     "disable-inverter-constraint" => true,
    #     "nprocs" => 1,
    # )


args = Dict{String,Any}()
# args["base_network"] = data
args["network"] = "network.json"
args["output"]       = "result.json"
args["solvers"] = Dict(
  "misocp_solver" => () -> CPLEX.Optimizer(),
  "nlp_solver"    => () -> CPLEX.Optimizer()
)
# here “settings” is *already* the Dict you want
# args["settings"] = "onm_settings.json"

args = Dict{String,Any}(
  "network"      => "network.json",
  "output"       => "result.json",
  # these three *are* runtime args, not network settings:
  "switch_model"  => "LPUBFDiagPowerModel",
  "switch_only"   => true,
  "dispatch_only" => false,
  "solvers"      => Dict(
    "misocp_solver" => () -> CPLEX.Optimizer(),
    "nlp_solver"    => () -> CPLEX.Optimizer()
  )
)

# Execute reconfiguration optimization
result = PowerModelsONM.entrypoint(args)

# # 
# args["network_data"] = JSON.parsefile(args["network"])
# args["solvers"] = Dict(
#     "misocp_solver" => () -> CPLEX.Optimizer(),
#     "nlp_solver"    => () -> CPLEX.Optimizer(),
# )
# args["settings"] = Dict(
#     "switch_model"  => "LPUBFDiagPowerModel",  # linearized unbalanced branch flow
#     "switch_only"   => true,
#     "dispatch_only" => false,
# )


# println("--- ARGS DUMP ---")
# for (k,v) in args
#     println(lpad(k, 12), " => ", typeof(v))
# end
# error("stopping here to inspect args")

# Execute reconfiguration optimization
result = PowerModelsONM.entrypoint(args)

# 4. --- Output optimal switch operations -------------------------------------
println("Optimal switch changes:")
display(result["output_data"]["Switch changes"])
