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

function permutate_initial_conditions(out_dir, adjacencies, disturbances; num_final_permutations=10)
    N = size(adjacencies)[1] - 1
    unique_permutations = 0
    # Nxm matrix where N is number of reactors and m is number of initial conditions
    original_values = repeat([300 388.7 0.11],N)
    # original_values = [300 370 0.055;300 380 0.08; 300 388.7 0.11] # 3R series
    # original_values = [300 370 0.055;300 388.7 0.11; 300 388.7 0.11] # 3R 2and1 parallel
    # original_values = [300 370 0.055;300 370 0.055; 300 388.7 0.11] # 3R mixing
    steps_each_side = 2
    step_size = [0,10,0.05]
    permutation_weights = zeros(Float64, (2*steps_each_side + 1,size(original_values)[2],N))
    for n in 1:N
        for i in 1:(2*steps_each_side + 1)
            for j in 1:size(original_values)[2]
                permutation_weights[i,j,n] = original_values[n,j] + ((i - (steps_each_side + 1)) * step_size[j])
            end
        end
    end
    display(permutation_weights)
# add

    num_permutations = 5^(size(original_values)[2]) - 1 # iterate in base 5 through all possible permutations

    # lists all previous displacements so in the event of a 0 step we don't repeat unneeded simulations
    # fill with random small number so that the original displacements of 0 gets put there too
    completed_permutations = fill(0.000001010110000001, (num_permutations+1,size(original_values)[2]))
    top = fill(0.0, (num_final_permutations,10))

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

    # for i in ProgressBar(10:20)
    for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        base_five = string(i, base=5, pad=size(original_values)[2])
        # println(base_five)
        current_values = deepcopy(original_values)
        for n in 1:N
            for j in 1:size(original_values)[2]
                char = base_five[j]
                current_values[n,j] = permutation_weights[parse(Int64, char) + 1, j,n]
            end
        end
        println("current_values=",current_values)
        current_displacements = current_values[1,:] .- original_values[1,:]

        # check to see if permutation already in top permutations
        already_exists = false
        for j in 1:num_permutations
            if completed_permutations[j,:] == current_displacements
                already_exists = true
            end
        end
        completed_permutations[i+1,:] = current_displacements
        if already_exists
            continue
        end
        unique_permutations += 1

        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        discrepancies = MPC_tracking(adjacencies, adjacencies,disturbances,[0 0;0 0;0 0],[0 0;0 0;0 0],1,1e7,1e7,1e-3,1e9,90,1000,[8 15],15,current_values
        ;tmax=5000, print=false)
        n += 1
        avg_xB += discrepancies[1]
        avg_T += discrepancies[2]
        avg_flow += discrepancies[3]
        avg_heat += discrepancies[4]
        avg_max_xB += discrepancies[5]

        # println("Discrepancies: $discrepancies")
        sum_discrepancies = sum(discrepancies) # basic sum
        # sum_discrepancies = xB_norm_const_test_1 * discrepancies[1] + discrepancies[2] # test 1
        # sum_discrepancies = xB_norm_const_test_2 * discrepancies[1] + discrepancies[2] + avg_max_xB_const_test_2 * discrepancies[5] # test 2
        # sum_discrepancies = discrepancies[5] # test 3
        # sum_discrepancies = xB_norm_const_test_1 * discrepancies[1] + discrepancies[2] + discrepancies[6] # test 4
        # sum_discrepancies = discrepancies[1] # just xB


        m = minimum(top[:,10]) # find the worst performing permutations
        if sum_discrepancies > m
            insert_index = findall(x -> x==m, top[:,10])[1]
            top[insert_index,1:3] .= current_values[1,:] .- original_values[1,:]
            top[insert_index,4:9] = discrepancies
            top[insert_index,10] = sum_discrepancies
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
    println("$(unique_permutations) unique permutations found!")
    for i in 1:num_permutations
        println(completed_permutations[i,:])
    end
    println("Average discrepancies")
    print("xB: $(avg_xB/n)\n")
    print("T: $(avg_T/n)\n")
    print("Heat: $(avg_heat/n)\n")
    print("Flow: $(avg_flow/n)\n")
    print("max xB: $(avg_max_xB/n)\n")

    top = round.(top,digits=9)
    display(top)
    println("writing top configurations to file")
    top_file = out_dir * "\\top_initial_conditions.txt"
    top_excel_file = out_dir * "\\top_initial_conditions.xlsx"
    touch(top_file)
    file = open(top_file, "w")
    column_names = ["T0", "Ts", "xBs", "xBtvt", "Tvt", "flowvt", "heatvt", "max_Tvt", "tt_stable", "PI"]
    # write to text file
    write(file, join(column_names, "\t") * "\n")
    writedlm(file, top)
    # write to excel file
    XLSX.writetable(top_excel_file, [top[:, i] for i in 1:size(top,2)], column_names, overwrite=true)
    close(file)

    return top
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
function permutate_setpoint(out_dir, n1, n2, Dist_T0, initial_conditions, reconfiguration_conditions)#=,
    reactors_to_permutate)=#
    # TODO rewrite the code below for permuting setpoints
    # TODO hardcode SetChange_xB for 0.16-0.36
    N = size(n1)[1] - 1
    unique_permutations = 0
    # Nxm matrix where N is number of reactors and m is number of initial conditions


    num_permutations = 9
    # num_permutations = 29
    SetChange_xB=zeros(1,N)
    SetChange_T=zeros(1,N)
    

    for i=1:N
        SetChange_xB[i] = reconfiguration_conditions[i,1] - initial_conditions[i,1]
        SetChange_T[i] = reconfiguration_conditions[i,2] - initial_conditions[i,2]
    end
 
    step_size = [0 0 -0.01]
    # step_size = [0 0 0.01]
    
    # normalizing constants make the different fields factor equally into the sums
    n = 0

    # for i in ProgressBar(10:20)
    for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        # println(base_five)

        unique_permutations += 1
        SetChange_xB = SetChange_xB .+ step_size
        println("SetChange_xB=",SetChange_xB)
        # print("SetChange_xB: " * string(SetChange_xB))
        image_name = (out_dir * "\\Perm_SetChange_xB" * string(SetChange_xB[end]) * ".pdf")
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        discrepancies = MPC_tracking(n1,n2,Dist_T0,SetChange_xB,SetChange_T,
        1,1e7,1e7,1e-3,1e9,90,1000,[8 15],15,initial_conditions;tmax=5000, print=false,
        save_plots=true, plot_name=image_name)
        n += 1


    end

    # TODO have MPC_tracking always write full data to file, in here copy that file
    # and rename to current input values

end

# changes the configuration in the middle of running, permutates only Dist_T0
# reactors to permutate is the list of reactor numbers to permutate xBs' on
function permutate_temp_in(out_dir, n1, n2, Dist_T0, initial_conditions, reconfiguration_conditions,
    reactors_to_permutate)
    N = size(n1)[1] - 1
    unique_permutations = 0
    # Nxm matrix where N is number of reactors and m is number of initial conditions


    num_permutations = 19
    SetChange_xB = [0.15 0.15 0.15] .* reactors_to_permutate
    step_size = repeat([0.01 0.01], N)
    SetChange_T = [reconfiguration_conditions[i,2] - initial_conditions[i,2] for i in 1:N]

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

    # for i in ProgressBar(10:20)
    for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        # println(base_five)

        unique_permutations += 1
        Dist_T0 = Dist_T0 .+ step_size
        # print("SetChange_xB: " * string(SetChange_xB))
        image_name = (out_dir * "\\Perm_Dist_TO" * string(Dist_T0[1]) * ".pdf")
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        discrepancies = MPC_tracking(n1,n2,Dist_T0,SetChange_xB,SetChange_T,
        1,1e7,1e7,1e-3,1e9,90,1000,[8 15],15,initial_conditions;tmax=5000, print=false,
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
