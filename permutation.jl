# file with the functions for permutating the

using Plots, JuMP, DifferentialEquations, NLsolve, BenchmarkTools, Ipopt
using MathOptInterface, Printf, ProgressBars, DelimitedFiles, Profile, XLSX
using DataFrames
include("reactor_reconfigure_simulator.jl")

function permutate_weights(out_dir, disturbances)
    original_weights = [1,1e7,1e7,1e-5,1e7]
    powers_each_side = 2
    permutation_weights = zeros(Float64, (2*powers_each_side + 1,length(original_weights)))
    for i in 1:(2*powers_each_side + 1)
        for j in 1:length(original_weights)
            permutation_weights[i,j] = original_weights[j]*exp10(i - (powers_each_side + 1))
        end
    end
    display(permutation_weights)
    top_ten = fill(typemax(Float64), (10,12))

    num_permutations = 5^5 - 1 # iterate in base 5 through all possible permutations

    # normalizing constants make the different fields factor equally into the sums
    n = 0
    avg_xB = 0
    xB_norm_const_test_1 = 11.085 / 1.971e-6
    xB_norm_const_test_2 = 11.085 / 1.971e-6
    avg_T = 0
    avg_flow = 0
    avg_heat = 0
    avg_max_xB = 0
    avg_max_xB_const_test_2 = 10.91 / 0.111299


    # for i in ProgressBar(1755:1765)
    for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        base_five = string(i, base=5, pad=5)
        # println(base_five)
        current_weights = original_weights
        for j in 1:length(base_five)
            char = base_five[j]
            current_weights[j] = permutation_weights[parse(Int64, char) + 1, j]
        end
        q_T = current_weights[1]
        q_xA = current_weights[2]
        q_xB = current_weights[3]
        r_heat = current_weights[4]
        r_flow = current_weights[5]
        # println(current_weights)
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        discrepancies = MPC_tracking([0 0 1 1;0 0 1 1], disturbances,q_T,q_xA,q_xB,r_heat,r_flow,90,1000,[8 15];tmax=5000, print=false)
        n += 1
        avg_xB += discrepancies[1]
        avg_T += discrepancies[2]
        avg_flow += discrepancies[3]
        avg_heat += discrepancies[4]
        avg_max_xB += discrepancies[5]

        # println("Discrepancies: $discrepancies")
        # sum_discrepancies = sum(discrepancies) # basic sum
        # sum_discrepancies = xB_norm_const_test_1 * discrepancies[1] + discrepancies[2] # test 1
        # sum_discrepancies = xB_norm_const_test_2 * discrepancies[1] + discrepancies[2] + avg_max_xB_const_test_2 * discrepancies[5] # test 2
        # sum_discrepancies = discrepancies[5] # test 3
        # sum_discrepancies = xB_norm_const_test_1 * discrepancies[1] + discrepancies[2] + discrepancies[6] # test 4
        sum_discrepancies = discrepancies[1] # like test 3 but instead with max temp
        m = maximum(top_ten[:,12])
        if sum_discrepancies < m
            insert_index = findall(x -> x==m, top_ten[:,12])[1]
            top_ten[insert_index,1:5] = current_weights
            top_ten[insert_index,6:11] = discrepancies
            top_ten[insert_index,12] = sum_discrepancies
        end
        # for j in 1:10
        #     # println("$(top_ten[j,10]) $sum_discrepancies")
        #     is_full = maximum(top_ten[:,10]) != Inf
        #     # println("$(top_ten[:,10]) $is_full")
        #     if top_ten[j,10] > sum_discrepancies
        #         !is_full && top_ten[j,10] != Inf && continue
        #         top_ten[insert_index,1:5] = current_weights
        #         top_ten[insert_index,6:9] = discrepancies
        #         top_ten[insert_index,10] = sum_discrepancies
        #         break
        #     end
        # end
    end
    println("Average discrepancies")
    print("xB: $(avg_xB/n)\n")
    print("T: $(avg_T/n)\n")
    print("Heat: $(avg_heat/n)\n")
    print("Flow: $(avg_flow/n)\n")
    print("max xB: $(avg_max_xB/n)\n")

    display(top_ten)
    println("writing top ten configurations to top_ten.txt")
    top_ten_file = out_dir * "\\top_ten.txt"
    touch(top_ten_file)
    file = open(top_ten_file, "w")
    writedlm(file, top_ten)
    close(file)

    originalDiscrepancy = MPC_tracking([0 0 1 1;0 0 1 1], disturbances,1,1e7,1e7,1e-3,1e9,90,1000,[8 15];tmax=5000,print=true) # no disturbance
    return top_ten
end

function save_profile_images_permutations(inputMatrix, disturbances, out_dir)
    count = 1
    for row in eachrow(inputMatrix)
        println(row)
        q_T = row[1]
        q_xA = row[2]
        q_xB = row[3]
        r_heat = row[4]
        r_flow = row[5]
        image_name = join(row[1:5], "_")
        image_name = (out_dir * "\\Perm" * string(count) * "_" * image_name * ".pdf")
        MPC_tracking([0 0 1 1;0 0 1 1], disturbances,q_T,q_xA,q_xB,r_heat,r_flow,90,1000,[8 15];tmax=5000, print=false, save_plots=true, plot_name=image_name)
        count += 1
    end
end

function permutate_initial_conditions(out_dir, n1, n2, Dist_T0, SetChange_xB, SetChange_T)
    num_permutations_T=80
    num_permutations_xB=20
    # Nxm matrix where N is number of reactors and m is number of initial conditions
    original_values = [300 388.7 0.11]
    # original_values = repeat([300 388.7 0.11],N)
    # original_values = [300 370 0.055;300 380 0.08; 300 388.7 0.11] # 3R series
    # original_values = [300 370 0.055;300 388.7 0.11; 300 388.7 0.11] # 3R 2and1 parallel
    # original_values = [300 370 0.055;300 370 0.055; 300 388.7 0.11] # 3R mixing

    step_size_T=[0 1 0]
    step_size_xB=[0 0 0.01]
    
    

    for i in ProgressBar(0:num_permutations_T)
        permutation_conditions1=original_values
        if i>=round(Int,0.5*num_permutations_T)
            permutation_conditions1 = original_values -step_size_T*(i-round(Int,0.5*num_permutations_T)+1)
        else
            permutation_conditions1 = original_values +step_size_T*i
        end
        for j in ProgressBar(0:num_permutations_xB)
            if j>=0.5*num_permutations_xB
                permutation_conditions2 = permutation_conditions1 .-step_size_xB*(j-round(Int,0.5*num_permutations_xB)+1)
            else
                permutation_conditions2 = permutation_conditions1 .+step_size_xB*j
            end
            println("permutation_conditions=",permutation_conditions2)
            if permutation_conditions2[3]<=0.008
                continue
            else
                image_name = (out_dir * "\\Perm_initialCondtion" * string(permutation_conditions2) * ".pdf")
                MPC_tracking(n1,n2,Dist_T0,SetChange_xB,SetChange_T,
                1,1e7,1e7,1e-5,1e7,90,1000,[8 15],15,permutation_conditions2;tmax=200, print=false,
                save_plots=true, plot_name=image_name)
            end
            
        end
    end
end

function save_profile_images_initial_conditions(inputMatrix, adjacencies, disturbances, out_dir)
    count = 1
    for row in eachrow(inputMatrix)
        image_name = join(row[1:3], "_")
        image_name = (out_dir * "\\Perm" * string(count) * "_" * image_name * ".pdf")
        N = size(adjacencies)[1] - 1
        original_values = repeat([300 388.7 0.11],N)
        initial_values = original_values .+ transpose(row[1:3])
        MPC_tracking(adjacencies, adjacencies,disturbances,[0 0;0 0;0 0],[0 0;0 0;0 0],1,1e7,1e7,1e-3,1e9,90,1000,[8 15],15,initial_values
        ;tmax=5000, print=false, save_plots=true, plot_name=image_name)
        count += 1
    end
end

# changes the configuration in the middle of running, permutates only xBs' from 0.26-0.46
# reactors to permutate is the list of reactor numbers to permutate xBs' on
function permutate_setpoint(out_dir, n1, n2, Dist_T0, SetChange_xB, SetChange_T)#=,
    reactors_to_permutate)=#
    # TODO rewrite the code below for permuting setpoints
    # TODO hardcode SetChange_xB for 0.16-0.36
    N = size(n1)[1] - 1
    unique_permutations = 0
    # Nxm matrix where N is number of reactors and m is number of initial conditions


    num_permutations = 9
    # num_permutations = 29

    println("SetChange_xB=",SetChange_xB)
    println("SetChange_T=",SetChange_T)

    step_size = [0; 0; -0.01]
    # step_size = [0; 0; 0.01]
    
    # normalizing constants make the different fields factor equally into the sums
    n = 0

    for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        # println(base_five)

        unique_permutations += 1
        SetChange_xB = SetChange_xB .+ step_size
        println("SetChange_xB=",SetChange_xB)
        image_name = (out_dir * "\\Perm_SetChange_xB" * string(SetChange_xB[end]) * ".pdf")
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        MPC_tracking(n1,n2,Dist_T0,SetChange_xB,SetChange_T,
        1,1e7,1e7,1e-5,1e7,90,1000,[8 15],0,[300 388.7 0.11];tmax=200, print=false,
        save_plots=true, plot_name=image_name)
        n += 1


    end

    # TODO have MPC_tracking always write full data to file, in here copy that file
    # and rename to current input values

end

# changes the configuration in the middle of running, permutates only Dist_T0
# reactors to permutate is the list of reactor numbers to permutate xBs' on
function permutate_temp_in(out_dir, n1, n2, Dist_T0, SetChange_xB, SetChange_T)#=,
    reactors_to_permutate)=#
    N = size(n1)[1] - 1
    unique_permutations = 0
    # Nxm matrix where N is number of reactors and m is number of initial conditions


    # num_permutations = 19
    num_permutations = 39

    # Disturbance happens on th 15th time point
    # step_size = repeat([0 -1], N)
    step_size = repeat([0 1], N)
 
    n = 0

    # for i in ProgressBar(10:20)
    for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        # println(base_five)

        unique_permutations += 1
        Dist_T0 = Dist_T0 .+ step_size
        println("Dist_T0=",Dist_T0)
        # print("SetChange_xB: " * string(SetChange_xB))
        image_name = (out_dir * "\\Perm_Dist_TO" * string(Dist_T0[1,2]) * ".pdf")
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        MPC_tracking(n1,n2,Dist_T0,SetChange_xB,SetChange_T,
        1,1e7,1e7,1e-5,1e7,90,1000,[8 15],15,[300 388.7 0.11];tmax=5000, print=false,
        save_plots=true, plot_name=image_name)
        n += 1

    end

    # TODO have MPC_tracking always write full data to file, in here copy that file
    # and rename to current input values

end

function save_profile_images_permutation_setpoints(inputMatrix, n1,n2, reactors_to_permutate, disturbances, out_dir)
    # TODO rewrite the code below for permuting setpoint
    count = 1
    for row in eachrow(inputMatrix)
        image_name = join(row[1:3], "_")
        image_name = (out_dir * "\\Perm" * string(count) * "_" * image_name * ".pdf")
        N = size(adjacencies)[1] - 1
        original_values = repeat([300 388.7 0.11],N)
        initial_values = original_values .+ transpose(row[1:3])
        MPC_tracking(n1,n2,disturbances,1,1e7,1e7,1e-3,1e9,90,1000,[8 15],initial_values;tmax=5000, print=false, save_plots=true, plot_name=image_name)
        count += 1
    end
end

# format: 1-3,2-3 (mixing configuration)
# a single dash (-) denotes parallel reactors
# 1-3 shows reactor connection from reactor 1 to reactor 3
# comma separates connections
function configuration_text_to_matrix(configuration_text, num_reactors)
    m = zeros(Int64, num_reactors + 1, num_reactors + 1)
    if configuration_text != "-"  # if not parallel
        comma_split = split(configuration_text, ",")
        for connection in comma_split
            dash_split = split(connection, "-")
            m[parse(Int64, dash_split[1]),parse(Int64, dash_split[2])] = 1
        end
    end

    # assuming all reactors have inputs and outputs
    for i in 1:num_reactors
        if m[i,:] == zeros(Int64, num_reactors + 1)
            m[i,num_reactors+1] = 1
        end
        m[num_reactors+1, i] = 1
    end
    return m
end

function configuration_matrix_to_text(configuration_matrix)
    num_reactors = size(configuration_matrix)[1] - 1
    has_connections = false
    configuration_text = ""
    for i in 1:num_reactors
        for j in 1:num_reactors
            if configuration_matrix[i,j] == 1
                if has_connections
                    configuration_text *= ","
                end
                has_connections = true
                configuration_text *= string(i) * "-" * string(j)
            end
        end
    end
    if length(configuration_text) == 0
        configuration_text *= "-"
    end
    return configuration_text
end
