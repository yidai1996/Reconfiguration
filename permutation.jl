# file with the functions for permutating the

using Plots, JuMP, DifferentialEquations, NLsolve, BenchmarkTools, Ipopt, MathOptInterface, Printf, ProgressBars, DelimitedFiles, Profile
# include("more_reactor_reconfigure_combination_xB_compacted_deleteuselesssyntax_v7.jl")

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
        image_name = (out_dir * "\\Perm" * string(count) * "_" * image_name * ".png")
        MPC_tracking([0 0 1 1;0 0 1 1], disturbances,q_T,q_xA,q_xB,r_heat,r_flow,90,1000,[8 15];tmax=5000, print=false, save_plots=true, plot_name=image_name)
        count += 1
    end
end

function permutate_initial_conditions(out_dir, adjacencies, disturbances)
    original_values = [300,388.7,0.11]
    steps_each_side = 2
    step_size = [0,25,0.2]
    permutation_weights = zeros(Float64, (2*steps_each_side + 1,length(original_values)))
    for i in 1:(2*steps_each_side + 1)
        for j in 1:length(original_values)
            permutation_weights[i,j] = original_values[j] + ((i - (steps_each_side + 1)) * step_size[j])
        end
    end
    display(permutation_weights)
    top_ten = fill(0.0, (10,10))

    num_permutations = 5^(length(original_values)) - 1 # iterate in base 5 through all possible permutations

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
        base_five = string(i, base=5, pad=length(original_values))
        # println(base_five)
        current_values = original_values
        for j in 1:length(base_five)
            char = base_five[j]
            current_values[j] = permutation_weights[parse(Int64, char) + 1, j]
        end
        T0 = current_values[1]
        Ts = current_values[2]
        xBs = current_values[3]

        # println(current_weights)
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        discrepancies = MPC_tracking(adjacencies, disturbances,1,1e7,1e7,1e-3,1e9,90,1000,[8 15]
        ;tmax=5000, print=false,initial_values=current_values)
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
        sum_discrepancies = discrepancies[1] # just xB
        m = minimum(top_ten[:,10]) # find the worst performing permutations
        if sum_discrepancies > m
            insert_index = findall(x -> x==m, top_ten[:,10])[1]
            top_ten[insert_index,1:3] = current_values
            top_ten[insert_index,4:9] = discrepancies
            top_ten[insert_index,10] = sum_discrepancies
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
    top_ten_file = out_dir * "\\top_ten_initial_conditions.txt"
    touch(top_ten_file)
    file = open(top_ten_file, "w")
    writedlm(file, top_ten)
    close(file)

    return top_ten
end

function save_profile_images_initial_conditions(inputMatrix, adjacencies, disturbances, out_dir)
    count = 1
    for row in eachrow(inputMatrix)
        println(row)
        T0 = row[1]
        Ts = row[2]
        xBs = row[3]
        image_name = join(row[1:3], "_")
        image_name = (out_dir * "\\Perm" * string(count) * "_" * image_name * ".png")
        MPC_tracking(adjacencies, disturbances,1,1e7,1e7,1e-3,1e9,90,1000,[8 15];tmax=5000, print=false, save_plots=true, plot_name=image_name,initial_values=row[1:3])
        count += 1
    end
end
