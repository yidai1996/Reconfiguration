import Graphs, GraphNeuralNetworks
import GraphNeuralNetworks.GNNGraph

# use Graphs.jl to construction the undirected graph of the three-reactor system
# assumption: every variable (x_{i,t}) is a node, so does everyone has a relationship with others? Or some nodes may not have direct relationship, e.g. mole fraction at time = 1, and the flow rate at time = 12


# Test of demo code
using LightGraphs
using Flux
using GraphNeuralNetworks
using SparseArrays
using GraphPlot # https://github.com/JuliaGraphs/GraphPlot.jl


g = SimpleGraph(19)

# Define the edges based on your description
# Nodes: 1: Tin1, 2:Tin2, 3:Tin3, 4:T1, 5:T2, 6:T3, 7:x1, 8:x2, 9:x3, 10:xtot, 11:F1, 12:F2, 13:F3, 14:Q1, 15:Q2, 16:Q3, 17:Tset1, 18:Tset2, 19:Tset3
self_defined_edges = [
    (7, 11),  # x1 - F1
    (7, 10),  # x1 - xtot
    (7, 14),  # x1 - Q1
    (7, 1),  # x1 - Tin1
    (7, 4),  # x1 - T1
    (14, 4),  # Q1 - T1
    (4, 1),   # T1 - Tin1
    (11, 10), # F1 - xtot
    (8, 11), # x2 - xtot
    (8, 12),  # x2 - F2
    (8, 2), # x2 - Tin2
    (8, 5), # x2 - T2 
    (8, 15),  # x2 - Q2
    (15, 5), # Q2 - T2 
    (5, 2), # T2 - Tin2 
    (12, 10), # F2 - xtot 
    (9, 10), # x3 - xtot
    (9, 13), # x3 - F3 
    (9, 3), # x3 - Tin3 
    (9, 6), # x3 - T3 
    (9, 16), # x3 - Q3 
    (16, 6), # Q3 - T3 
    (6, 3), # T3 - Tin3 
    (13, 10),  # F3 - xtot
    # Add Tset
    (4, 17), # T1 - T1set
    (5, 18), # T2 - T2set
    (6, 19), # T3 - T3set
    (17, 7), # T1set - x1
    (18, 8), # T2set - x2 
    (19, 9), # T3set - x3 
    (17, 11), # T1set - F1
    (18, 12), # T2set - F2 
    (19, 13), # T3set - F3 
    (17, 14), # T1set - Q1
    (18, 15), # T2set - Q2 
    (19, 16), # T3set - Q3 
]

# Add the edges to the graph
for (u, v) in self_defined_edges
    add_edge!(g, u, v)
end

adj_matrix = LightGraphs.adjacency_matrix(g)

num_nodes = nv(g)

using CSV, DataFrames

df_features_setpoint_T1 = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_10.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
df_features_setpoint_T2 = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_10.0SetChange_T3_0.0.csv", DataFrame)
df_features_setpoint_T3 = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_10.0.csv", DataFrame)

df_features_disturbance_12_same = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_310.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
df_features_disturbance_23_same = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_310.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
df_features_disturbance_13_same = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_310.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

df_features_all_same = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
df_features_combine = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_310.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_10.0.csv", DataFrame)

node_features_setpoint_T1 = Matrix{Float64}(df_features_setpoint_T1[14:20,3:21])
node_features_setpoint_T2 = Matrix{Float64}(df_features_setpoint_T2[14:20,3:21])
node_features_setpoint_T3 = Matrix{Float64}(df_features_setpoint_T3[14:20,3:21])
node_features_disturbance_12 = Matrix{Float64}(df_features_disturbance_12_same[14:20,3:21])
node_features_disturbance_23 = Matrix{Float64}(df_features_disturbance_23_same[14:20,3:21])
node_features_disturbance_13 = Matrix{Float64}(df_features_disturbance_13_same[14:20,3:21])
node_features_all_same = Matrix{Float64}(df_features_all_same[14:20,3:21]) 
node_features_combine = Matrix{Float64}(df_features_combine[14:20,3:21])

# Normalize features
using LinearAlgebra
# Define a function to normalize each column of the matrix
function normalize_columns_between_0_and_1(matrix::Array{Float64, 2})
    normalized_matrix = copy(matrix) # Create a copy to avoid modifying the original matrix
    num_rows = size(matrix, 1) # Number of rows
    num_cols = size(matrix, 2) # Number of columns

    for j in 1:num_cols
        col_min = minimum(matrix[:, j]) # Find the minimum value in the column
        col_max = maximum(matrix[:, j]) # Find the maximum value in the column

        if col_min != col_max
            # Scale values to the range 0 to 1
            normalized_matrix[:, j] .= (matrix[:, j] .- col_min) ./ (col_max - col_min)
        else
            # If col_min == col_max, all values in column are the same. Set them to 0.
            normalized_matrix[:, j] .= 0
        end
    end
    return normalized_matrix
end
node_features_norm_setpoint_T1 = transpose(normalize_columns_between_0_and_1(node_features_setpoint_T1))
node_features_norm_setpoint_T2 = transpose(normalize_columns_between_0_and_1(node_features_setpoint_T2))
node_features_norm_setpoint_T3 = transpose(normalize_columns_between_0_and_1(node_features_setpoint_T3))
node_features_norm_disturbance_12 = transpose(normalize_columns_between_0_and_1(node_features_disturbance_12))
node_features_norm_disturbance_23 = transpose(normalize_columns_between_0_and_1(node_features_disturbance_23))
node_features_norm_disturbance_13 = transpose(normalize_columns_between_0_and_1(node_features_disturbance_13))
node_features_norm_all_same = transpose(normalize_columns_between_0_and_1(node_features_all_same))
node_features_norm_combine = transpose(normalize_columns_between_0_and_1(node_features_combine))

function normalize_rows_between_0_and_1(matrix::Array{Float64, 2})
    normalized_matrix = copy(matrix) # Create a copy to avoid modifying the original matrix
    num_rows = size(matrix, 1) # Number of rows
    num_cols = size(matrix, 2) # Number of columns

    for j in 1:num_rows
        rows_min = minimum(matrix[j, :]) # Find the minimum value in the column
        rows_max = maximum(matrix[j, :]) # Find the maximum value in the column

        if rows_min != rows_max
            # Scale values to the range 0 to 1
            normalized_matrix[j, :] .= (matrix[j, :] .- rows_min) ./ (rows_max - rows_min)
        else
            # If col_min == col_max, all values in column are the same. Set them to 0.
            normalized_matrix[j, :] .= 0
        end
    end
    return normalized_matrix
end

function get_label(node_features_norm)
    node_labels = zeros(num_nodes, num_nodes)
    for i in 1:num_nodes
        for j in 1:num_nodes
            similarity = norm(node_features_norm[i,:]-node_features_norm[j,:])
            # println(i," and ",j)
            # println("Similarity = ", similarity)
        
            if similarity < 10e-3
                node_labels[i, j] = 1
            end
        end
    end

    return node_labels
end

function generate_pairs(similarity_labels)
    pairs = []  # Initialize an empty array to store pairs
    num_nodes = size(similarity_labels, 1)

    for i in 1:num_nodes
        for j in 1:num_nodes
            if i != j  # Ignoring self-loops
                similarity = similarity_labels[i, j]
                push!(pairs, (i, j, similarity))
            end
        end
    end
    
    return pairs
end

pairs_setpoint_T1 = generate_pairs(get_label(node_features_norm_setpoint_T1))
pairs_setpoint_T2 = generate_pairs(get_label(node_features_norm_setpoint_T2))
pairs_setpoint_T3 = generate_pairs(get_label(node_features_norm_setpoint_T3))
pairs_disturbance12 = generate_pairs(get_label(node_features_norm_disturbance_12))
pairs_disturbance23 = generate_pairs(get_label(node_features_norm_disturbance_23))
pairs_disturbance13 = generate_pairs(get_label(node_features_norm_disturbance_13))
pairs_all_same = generate_pairs(get_label(node_features_norm_all_same))
pairs_combine = generate_pairs(get_label(node_features_norm_combine))

# From GPT
function gcn_layer(W, A, X)
    return relu.(A * (X * W))
end

using Flux:params
# Define a GCN with four hidden layers and embedding output
function GCN(A, X, num_hidden, embed_dim)
    W1 = randn(Float32, size(X, 2), num_hidden)
    W2 = randn(Float32, num_hidden, num_hidden)
    W3 = randn(Float32, num_hidden, num_hidden)
    W4 = randn(Float32, num_hidden, num_hidden)
    W_out = randn(Float32, num_hidden, embed_dim)

    # Collect all the parameters for optimization
    params = Flux.params(W1, W2, W3, W4, W_out)
    # params = Flux.params(W1, W4, W_out)

    function forward(x)
        h1 = gcn_layer(W1, A, x)
        # println("h1=", h1)
        h2 = gcn_layer(W2, A, h1)
        # println("h2=",h2)
        h3 = gcn_layer(W3, A, h2)
        # println("h3=",h3)
        h4 = gcn_layer(W4, A, h3)
        # println("h4=",h4)
        return A * (h4 * W_out)  # Output embeddings
    end
    
    return forward, params
end



# Parameters
num_hidden = 8
embed_dim = 5  # Desired embedding dimension

# Build GCN model
forward, parameters = GCN(adj_matrix, node_features_norm_all_same, num_hidden, embed_dim)

# Function to calculate cosine similarity
function cosine_similarity(a, b)
    dot(a, b) / (norm(a) * norm(b))
end

# Helper function to clamp and make sure the similarity is positive
function clamped_positive_similarity(similarity)
    return max(min(similarity, 1 - 1e-7), 1e-7)
end

# Define the similarity loss using cross-entropy
function similarity_loss(embeddings, pairs)
    total_loss = 0.0
    for (i, j, label) in pairs
        node_i = embeddings[i, :]
        node_j = embeddings[j, :]
        similarity = cosine_similarity(node_i, node_j)
        # if similarity > 0.999
        #     y_pred = 1.0
        # else y_pred = 0.0
        # end
        
        # Compute the binary cross-entropy loss
        y_pred = clamped_positive_similarity(similarity)
        y_true = Float32(label)
        # println(y_pred)
        total_loss += -y_true * log(sigmoid(y_pred)) - (1 - y_true) * log(1 - sigmoid(y_pred))
        # total_loss += Flux.logitbinarycrossentropy(y_pred, y_true)
    end
    return total_loss / length(pairs)
end

function multi_task_loss(embedding1, pairs1, embedding2, pairs2, embedding3, pairs3, embedding4, pairs4, embedding5, pairs5, embedding6, pairs6, embedding7, pairs7, embedding8, pairs8)
    loss1 = similarity_loss(embedding1, pairs1)
    loss2 = similarity_loss(embedding2, pairs2)
    loss3 = similarity_loss(embedding3, pairs3)
    loss4 = similarity_loss(embedding4, pairs4)
    loss5 = similarity_loss(embedding5, pairs5)
    loss6 = similarity_loss(embedding6, pairs6)
    loss7 = similarity_loss(embedding7, pairs7)
    loss8 = similarity_loss(embedding8, pairs8)

    # total_loss = (loss1 + loss2 + loss3 + loss4 + loss5 + loss6)/6
    # total_loss = (loss1 + loss2 + loss3 + loss4)/4
    # total_loss = ((loss1 + loss2 + loss3)/3 + (loss4 + loss5 + loss6)/3 + loss7 + loss8)/4
    total_loss = (loss1 + loss2 + loss3 + loss4 + loss5 + loss6+ loss7 + loss8)/8
    return total_loss
end
# Define optimizer
lr = 0.001
opt = Adam(lr)

using Flux: gradient

# Training loop epochs = 20, 50, 100
n_epochs = 200
for epoch in 1:n_epochs
    # Compute gradients and perform optimization step
    loss_value, back = Flux.withgradient(parameters) do
        embeddings_setpoint_T1 = forward(node_features_norm_setpoint_T1)
        embeddings_setpoint_T2 = forward(node_features_norm_setpoint_T2)
        embeddings_setpoint_T3 = forward(node_features_norm_setpoint_T3)
        embeddings_disturbance12 = forward(node_features_norm_disturbance_12)
        embeddings_disturbance23 = forward(node_features_norm_disturbance_23)
        embeddings_disturbance13 = forward(node_features_norm_disturbance_13)
        embeddings_all_same = forward(node_features_norm_all_same)
        embeddings_combine = forward(node_features_norm_combine)

        multi_task_loss(embeddings_setpoint_T1, pairs_setpoint_T1, embeddings_setpoint_T2, pairs_setpoint_T2, embeddings_setpoint_T3, pairs_setpoint_T3, embeddings_disturbance12, pairs_disturbance12, embeddings_disturbance23, pairs_disturbance23, embeddings_disturbance13, pairs_disturbance13, embeddings_all_same, pairs_all_same, embeddings_combine, pairs_combine)
        # multi_task_loss(embeddings_setpoint, pairs_setpoint, embeddings_disturbance12, pairs_disturbance12, embeddings_all_same, pairs_all_same, embeddings_combine, pairs_combine)
        # similarity_loss(embeddings_disturbance, pairs_disturbance)
    end
    
    # Perform the optimization step by updating the parameters
    Flux.Optimise.update!(opt, parameters, back)
    
    if epoch % 10 == 0
        println("Epoch $epoch: Loss = $loss_value")
    end
end

# Get the test data 
# all same
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

# disturbance on R1 
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_310.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

# disturbance on R2 
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_310.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

# disturbance on R3 
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_310.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

# setpoint tracking_T1
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_10.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
# setpoint tracking_T2
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_10.0SetChange_T3_0.0.csv", DataFrame)
# # setpoint tracking_T3
# df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_10.0.csv", DataFrame)

# combine
df_features_test = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_310.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_10.0.csv", DataFrame)

node_features_test = Matrix{Float64}(df_features_test[14:20,3:21])
node_features_norm_test = transpose(normalize_columns_between_0_and_1(node_features_test))

# get the test loss and the symmetry detection result
result_embedding = forward(node_features_norm_test)
pairs_test = generate_pairs(get_label(node_features_norm_test))
similarity_loss_test = similarity_loss(result_embedding, pairs_test)
# Nodes: 1: xBset, 2: Tin1, 3:Tin2, 4:Tin3, 5:T1, 6:T2, 7:T3, 8:x1, 9:x2, 10:x3, 11:xtot, 12:F1, 13:F2, 14:F3, 15:Q1, 16:Q2, 17:Q3
# Figure out identical units (before GNN and after GNN)
# Before GNN: calculatue the cosine_similarity for node_features_norm
# comparison between R1 and R2
# xB, T, F, Q, Tin
xB_12 = cosine_similarity(node_features_norm_test[7,:], node_features_norm_test[8,:])
T_12 = cosine_similarity(node_features_norm_test[4,:], node_features_norm_test[5,:])
F_12 = cosine_similarity(node_features_norm_test[11,:], node_features_norm_test[12,:])
Q_12 = cosine_similarity(node_features_norm_test[14,:], node_features_norm_test[15,:])
Tin_12 = cosine_similarity(node_features_norm_test[1,:], node_features_norm_test[2,:])
# comparison between R1 and R3
xB_13 = cosine_similarity(node_features_norm_test[7,:], node_features_norm_test[9,:])
T_13 = cosine_similarity(node_features_norm_test[4,:], node_features_norm_test[6,:])
F_13 = cosine_similarity(node_features_norm_test[11,:], node_features_norm_test[13,:])
Q_13 = cosine_similarity(node_features_norm_test[14,:], node_features_norm_test[16,:])
Tin_13 = cosine_similarity(node_features_norm_test[1,:], node_features_norm_test[3,:])
# comparison between R2 and R3
xB_23 = cosine_similarity(node_features_norm_test[8,:], node_features_norm_test[9,:])
T_23 = cosine_similarity(node_features_norm_test[5,:], node_features_norm_test[6,:])
F_23 = cosine_similarity(node_features_norm_test[12,:], node_features_norm_test[13,:])
Q_23 = cosine_similarity(node_features_norm_test[15,:], node_features_norm_test[16,:])
Tin_23 = cosine_similarity(node_features_norm_test[2,:], node_features_norm_test[3,:])

# After GNN: calculatue the cosine_similarity for node_features_norm
# xB, T, F, Q, Tin
# NOTE: don't Q and F in the GNN embedded MPC
xB_12 = cosine_similarity(result_embedding[7,:], result_embedding[8,:])
T_12 = cosine_similarity(result_embedding[4,:], result_embedding[5,:])
F_12 = cosine_similarity(result_embedding[11,:], result_embedding[12,:])
Q_12 = cosine_similarity(result_embedding[14,:], result_embedding[15,:])
Tset_12 = cosine_similarity(result_embedding[17,:], result_embedding[18,:])
# comparison between R1 and R3
xB_13 = cosine_similarity(result_embedding[7,:], result_embedding[9,:])
T_13 = cosine_similarity(result_embedding[4,:], result_embedding[6,:])
F_13 = cosine_similarity(result_embedding[11,:], result_embedding[13,:])
Q_13 = cosine_similarity(result_embedding[14,:], result_embedding[16,:])
Tset_13 = cosine_similarity(result_embedding[17,:], result_embedding[19,:])
# comparison between R2 and R3
xB_23 = cosine_similarity(result_embedding[8,:], result_embedding[9,:])
T_23 = cosine_similarity(result_embedding[5,:], result_embedding[6,:])
F_23 = cosine_similarity(result_embedding[12,:], result_embedding[13,:])
Q_23 = cosine_similarity(result_embedding[15,:], result_embedding[16,:])
Tset_23 = cosine_similarity(result_embedding[18,:], result_embedding[19,:])

# Save the model
using BSON: @save, @load

# @save "trained_gnn_8cases_new_200epochs_19nodes.bson" parameters