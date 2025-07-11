# Case study from: https://ieeexplore.ieee.org/document/7180413
# Case study from: https://ieeexplore.ieee.org/abstract/document/25627

using JuMP
using CPLEX
using MathOptInterface
# using BARON


function bus_system_reconfiguration()
    println("new simulation starts") # Print a message indicating the simulation has begun
    # -----------------------------
    # 1. Define sets and structure
    # -----------------------------
    buses = 1:33 # Bus indices from 1 to 33
    T = 1:24  # 24-hour time periods

    # List of all branches (edges) in the network, including tie switches
    branches = [
        (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9), (9,10),
        (10,11), (11,12), (12,13), (13,14), (14,15), (15,16), (16,17),
        (17,18), (2,19), (19,20), (20,21), (21,22), (3,23), (23,24),
        (24,25), (6,26), (26,27), (27,28), (28,29), (29,30), (30,31),
        (31,32), (32,33)
    ]

    tie_switches = [(8,21), (9,15), (12,22), (18,33), (25,29)]

    M = 1000000 # Big-M constant for linear relaxations
    # -----------------------------
    # 2. Parameters and DG info
    # -----------------------------

    # Series resistance R_ij for each branch (active loss coefficient) https://ieeexplore.ieee.org/abstract/document/25627
    R = Dict(
        (1,2)=>0.0922, (2,3)=>0.4930, (3,4)=>0.3660, (4,5)=>0.3811, (5,6)=>0.8190,
        (6,7)=>0.1872, (7,8)=>0.7114, (8,9)=>1.0300, (9,10)=>1.0440, (10,11)=>0.1966,
        (11,12)=>0.3744, (12,13)=>1.4680, (13,14)=>0.5416, (14,15)=>0.5910, (15,16)=>0.7463,
        (16,17)=>1.2890, (17,18)=>0.7320, (2,19)=>0.1640, (19,20)=>1.5040, (20,21)=>0.4095,
        (21,22)=>0.7089, (3,23)=>0.4512, (23,24)=>0.8980, (24,25)=>0.8960, (6,26)=>0.2030,
        (26,27)=>0.2842, (27,28)=>1.0590, (28,29)=>0.8042, (29,30)=>0.5075, (30,31)=>0.9744,
        (31,32)=>0.3105, (32,33)=>0.3410
    )

    # Series reactance X_ij for each branch (reactive loss coefficient) https://ieeexplore.ieee.org/abstract/document/25627. 
    # X and R are used to calculate voltage drop along the line (i,j) and losses (active loss and reactive loss) on that branch.
    X = Dict(
        (1,2)=>0.0470, (2,3)=>0.2511, (3,4)=>0.1864, (4,5)=>0.1941, (5,6)=>0.7070,
        (6,7)=>0.6188, (7,8)=>0.2351, (8,9)=>0.7400, (9,10)=>0.7400, (10,11)=>0.0650,
        (11,12)=>0.1238, (12,13)=>1.1550, (13,14)=>0.7129, (14,15)=>0.5260, (15,16)=>0.5450,
        (16,17)=>1.7210, (17,18)=>0.5740, (2,19)=>0.1565, (19,20)=>1.3554, (20,21)=>0.4784,
        (21,22)=>0.9373, (3,23)=>0.3083, (23,24)=>0.7091, (24,25)=>0.7011, (6,26)=>0.1034,
        (26,27)=>0.1447, (27,28)=>0.9337, (28,29)=>0.7006, (29,30)=>0.2585, (30,31)=>0.9630,
        (31,32)=>0.3619, (32,33)=>0.5302
    )

    # Maximum branch flow limits Pmax_ij (kW)
    Pmax = Dict((i,j) => 500 for (i,j) in branches)
    Pmax[(1,2)] = 5084 # greatly enlarge capacity of branch (1,2) since it's the substation

    # PL_base: active power demand at bus i (https://ieeexplore.ieee.org/abstract/document/25627)
    PL_base = Dict(
        2=>100, 3=>90, 4=>120, 5=>60, 6=>60,
        7=>200, 8=>200, 9=>60, 10=>60, 11=>45,
        12=>60, 13=>60, 14=>120, 15=>60, 16=>60,
        17=>60, 18=>90, 19=>90, 20=>90, 21=>90,
        22=>90, 23=>90, 24=>420, 25=>420, 26=>60,
        27=>60, 28=>60, 29=>120, 30=>200, 31=>150,
        32=>210, 33=>60
    ) 

    PL = Dict((i, t) => get(PL_base, i, 0.0) for i in buses for t in T) # adding time dimension 

    # Maximum branch flow limits Qmax_ij (kW)
    Qmax = Dict((i,j) => 200 for (i,j) in branches)
    Qmax[(1,2)] = 2547 # greatly enlarge capacity of branch (1,2) since it's the substation
    # QL_base: reactive power demand at bus i (https://ieeexplore.ieee.org/abstract/document/25627)
    QL_base = Dict(
        2=>60, 3=>40, 4=>80, 5=>30, 6=>20,
        7=>100, 8=>100, 9=>20, 10=>20, 11=>30,
        12=>35, 13=>35, 14=>80, 15=>10, 16=>20,
        17=>20, 18=>40, 19=>40, 20=>40, 21=>40,
        22=>40, 23=>50, 24=>200, 25=>200, 26=>25,
        27=>25, 28=>20, 29=>70, 30=>600, 31=>70,
        32=>100, 33=>40
    )
    QL = Dict((i, t) => get(QL_base, i, 0.0) for i in buses for t in T) # adding time dimension 

    # Maximum DG output at each bus (kW)
    P_DG = Dict((i,t) => 0.0 for i in buses for t in T)
    Q_DG = Dict((i,t) => 0.0 for i in buses for t in T)
    for t in T
        P_DG[(1,t)] = 5084
        Q_DG[(1,t)] = 2547
    end 
    
    DG_max = Dict(i => 0.0 for i in buses)
    # # Additional DG for renewable buses
    # DG_max[10] = 400.0  # wind
    # DG_max[33] = 500.0  # wind
    # DG_max[7]  = 350.0  # solar
    # DG_max[14] = 450.0  # solar

    # # Renewable generation load variation (constant now, should be replaced by load variation data when no other bugs)
    # forecast_DG = Dict{Tuple{Int, Int}, Float64}()
    # for t in T
    #     forecast_DG[(10, t)] = 400 # wind pattern
    #     forecast_DG[(33, t)] = 500 
    #     forecast_DG[(7, t)]  = 350  # solar pattern
    #     forecast_DG[(14, t)] = 450 
    # end

    # Initial voltage‐squared guesses for warm start (from paper https://ieeexplore.ieee.org/abstract/document/25627)
    V_sq_data = Dict(
        1=>1, 2 => 0.9927, 3=>0.9574, 4=>0.9374, 5=>0.9176, 6=>0.8707, 7=>0.8641, 
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
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 1)


    # --- Decision variables ---
    @variable(model, P[(i,j) in branches, t in T])
    @variable(model, Q[(i,j) in branches, t in T])
    @variable(model, V_sq[i in buses, t in T], lower_bound=0.81, upper_bound=1.21)
    @variable(model, scm[i in branches, t in T] >= 0) 

    # -----------------------------
    # 4. Constraints
    # -----------------------------
    @variable(model, P_slack[t in T] >= 0)
    @variable(model, Q_slack[t in T])
    for t in T
        @constraint(model, V_sq[1,t] == 1.0)

        
        for i in buses
            @constraint(model,
                sum(P[(i,j),t] for (i,j) in branches if i == 1)
                - sum(P[(i,j),t] for (i,j) in branches if j == 1)
                == P_slack[t])

            # if i == 1
            #     @constraint(model,
            #         sum(P[(i,j),t] for (i,j) in branches if i == 1)
            #     - sum(P[(j,i),t] for (j,i) in branches if i == 1)
            #     == P_slack[t]
            #     )
            #     @constraint(model,
            #         sum(Q[(i,j),t] for (i,j) in branches if i == 1)
            #     - sum(Q[(j,i),t] for (j,i) in branches if i == 1)
            #     == Q_slack[t]
            #     )
            
            # else
                @constraint(model,
                sum(Q[(k,i),t] for (k,j) in branches if j==i) -
                sum(Q[(i,k),t] for (j,k) in branches if j==i) ==
                Q_DG[i,t] - QL[(i,t)])
                @constraint(model,
                sum(P[(k,i),t] for (k,j) in branches if j==i) -
                sum(P[(i,k),t] for (j,k) in branches if j==i) ==
                P_DG[i,t] - PL[(i,t)])
            # end
        end

        for (i,j) in branches
            @constraint(model,
                V_sq[j,t] == V_sq[i,t]
                            - 2*(R[(i,j)]*P[(i,j),t] + X[(i,j)]*Q[(i,j),t])
                            + (R[(i,j)]^2 + X[(i,j)]^2)*scm[(i,j),t])
        end
    end

    
    # 5. Objective function: minimize total active loss
    # -----------------------------
    @objective(model, Min, sum(R[(i,j)]*scm[(i,j),t] for (i,j) in branches for t in T))


    # -----------------------------
    # 6. Solve and print/report
    # -----------------------------
    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        println("Objective value: ", objective_value(model))
    else
        println("No solution found. Termination: ", termination_status(model))
    end

    # println("Objective value (Total loss): ", objective_value(model))
    # for t in T
        # for (i,j) in branches
        #     println("Switch ($i,$j): status = ", value(alpha[(i,j),1]))
        #     println("Switch ($i,$j): P = ", value(P[(i,j),1]))
        #     println("Switch ($i,$j): Q = ", value(Q[(i,j),1]))
        #     println("Switch ($i,$j): scm = ", value(scm[(i,j),1]))
        #     println("Switch ($i,$j): V_sq = ", value(V_sq[(i),1]))
        # end
    # # end
    # for t in T, i in buses
    #     if get(DG_max, i, 0.0) > 0
    #         println("Bus $i at hour $t: DG output = ", value(P_DG[i,t]))
    #     end
    # end

end



bus_system_reconfiguration()