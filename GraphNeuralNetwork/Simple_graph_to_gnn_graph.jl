
using Graphs
using GraphNeuralNetworks
# using GeometricFlux
using Flux

g = SimpleGraph(17)

# Define the edges based on your description
self_defined_edges = [
    (1, 4),  # x1 - F1
    (1, 16),  # x1 - xtot
    (1, 7),  # x1 - Q1
    (1, 13),  # x1 - Tin1
    (1, 10),  # x1 - T1
    (7, 10),  # Q1 - T1
    (10, 13),   # T1 - Tin
    (4, 16), # F1 - xtot
    (2, 16), # x2 - xtot
    (2, 5),  # x2 - F2
    (2, 14), # x2 - Tin2
    (2, 11), # x2 - T2 
    (2, 8),  # x2 - Q2
    (8, 11), # Q2 - T2 
    (11, 14), # T2 - Tin2 
    (5, 16), # F2 - xtot 
    (3, 16), # x3 - xtot
    (3, 6), # x3 - F3 
    (3, 15), # X3 - Tin3 
    (3, 12), # x3 - T3 
    (3, 9), # x3 - Q2 
    (9, 12), # Q2 - T3 
    (12, 15), # T3 - Tin3 
    (6, 16), # F3 - xtot
    (16, 17), # xtot - xset
    (10, 18), #T1 - T1set
    (11, 19), #T2 - T2set
    (12, 20), #T3 - T3set
]

for (u, v) in self_defined_edges
    add_edge!(g, u, v)
end
# Step 2: Extract Components
edge_list = collect(Graphs.edges(g))
# Extract source and destination nodes from edges
src_nodes = [src(e) for e in edge_list]
dst_nodes = [dst(e) for e in edge_list]
num_nodes = nv(g)

# Step 3: Construct the GNNGraph
# Create adjacency tuple required for GNNGraph
num_edges = length(edge_list)
adjacency_list = (src_nodes, dst_nodes, nothing)

# Step 2: Prepare the Adjacency Matrix for GNN
# Create the adjacency matrix for the graph
adj_matrix = Graphs.adjacency_matrix(g)
using SparseArrays
# Convert the adjacency matrix to a sparse matrix (required by GeometricFlux)
adj_sparse = SparseArrays.sparse(adj_matrix)

# Step 3: Implement and Train a GNN
node_features = rand(Float32, nv(g), 5)  # Example node features (5 features per node)
hidden_dim = 10
target_dim = 4

# Define the GNN model, placing activation functions separately
# din, h_d_1, h_d_2, dout = 5, 10, 5, 2
# chain = GNNChain(GCNConv(din => h_d_1, relu),
#                  GCNConv(h_d_1 => h_d_2, relu),
#                  GCNConv(h_d_2 => dout))

# simple example
# din, d, dout = 3, 4, 2
# chain = GNNChain(GCNConv(din => d, relu),
#                  GCNConv(d => dout))
# g = rand_graph(10, 30)
# model = WithGraph(chain, g)
# X = randn(Float32, din, 10)
# # Pass only X as input, the model already contains the graph.
# y = model(X) 

using CSV, DataFrames

sample = CSV.read("C:\\Users\\10060\\Desktop\\SetChange_xB=0.01.csv", DataFrame)
df1=Matrix(sample[!,2:18]) # nodes are disordered
df_test = [df1[:,8] df1[:,9] df1[:,10] df1[:,12] df1[:,13] df1[:,14] df1[:,15] df1[:,16] df1[:,17] df1[:,5] df1[:,6] df1[:,7] df1[:,2] df1[:,3] df1[:,4] df1[:,11] df1[:,1] ]
# Load the LinearAlgebra package
using LinearAlgebra
# Define a function to normalize each column of the matrix
function normalize_columns_between_neg1_and_1(matrix::Array{Float64, 2})
    normalized_matrix = copy(matrix) # Create a copy to avoid modifying the original matrix
    num_rows = size(matrix, 1) # Number of rows
    num_cols = size(matrix, 2) # Number of columns

    for j in 1:num_cols
        col_min = minimum(matrix[:, j]) # Find the minimum value in the column
        col_max = maximum(matrix[:, j]) # Find the maximum value in the column

        if col_min != col_max
            # Scale values to the range -1 to 1
            normalized_matrix[:, j] .= 2 * ((matrix[:, j] .- col_min) ./ (col_max - col_min)) .- 1
        else
            # If col_min == col_max, all values in column are the same. Set them to 0.
            normalized_matrix[:, j] .= 0
        end
    end
    return normalized_matrix
end

df_norm = normalize_columns_between_neg1_and_1(df_test)
# X = transpose(df_test)
using LightGraphs
using Flux
using GraphNeuralNetworks
using SparseArrays
using GraphPlot # https://github.com/JuliaGraphs/GraphPlot.jl
din, d, dout = 3, 3, 2
chain = GNNChain(GCNConv(din => d, relu),
                 GCNConv(d => dout))
g_test = GNNGraph(adj_matrix )
model = WithGraph(chain, g_test)

# Pass only X as input, the model already contains the graph.
y = model(df_norm) 

# Calculate the similarity using cosine
# node 1 = y[:,1]
# node 1 and 2
cosine_similarity_1_2 = dot(y[:,1], y[:,2]) / (norm(y[:,1]) * norm(y[:,2]))
# node 3 and 4
cosine_similarity_3_4 = dot(y[:,3], y[:,4]) / (norm(y[:,3]) * norm(y[:,4]))
# node 2 and 3
cosine_similarity_2_3 = dot(y[:,3], y[:,2]) / (norm(y[:,3]) * norm(y[:,2]))
# node 1 and 3
cosine_similarity_1_3 = dot(y[:,3], y[:,1]) / (norm(y[:,3]) * norm(y[:,1]))
# node 10 and 11
cosine_similarity_10_11 = dot(y[:,10], y[:,11]) / (norm(y[:,10]) * norm(y[:,11]))
# node 10 and 1
cosine_similarity_10_1 = dot(y[:,10], y[:,1]) / (norm(y[:,10]) * norm(y[:,1]))
# node 16 and 1
cosine_similarity_16_1 = dot(y[:,16], y[:,1]) / (norm(y[:,16]) * norm(y[:,1]))
# node 16 and 1
cosine_similarity_16_1 = dot(y[:,16], y[:,1]) / (norm(y[:,16]) * norm(y[:,1]))
# node 5 and 11
cosine_similarity_5_11 = dot(y[:,5], y[:,11]) / (norm(y[:,5]) * norm(y[:,11]))
# node 5 and 4
cosine_similarity_5_4 = dot(y[:,5], y[:,4]) / (norm(y[:,5]) * norm(y[:,4]))
# ======
# The following is to train a GNN
# The loss function has error when run the training loop
# loss(x, y) = Flux.Losses.mse(model(x, adj_sparse), y)

x = node_features
y = rand(Float32, nv(g), 2)  # Dummy target labels

opt = Descent(0.01)  # Optimizer

# Training loop
for epoch in 1:100
    gs = gradient(Flux.params(model)) do
        l = loss(x, y)
        println("Epoch $epoch, Loss: $l")
        return l
    end
    Flux.Optimise.update!(opt, Flux.params(model), gs)
end

