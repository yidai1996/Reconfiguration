using Graphs
using Flux
using GeometricFlux
using Statistics: mean

# Construct a graph
g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)
add_edge!(g, 5, 1)

# Example node features for the 5 nodes (can be replaced by your data)
node_features = rand(5, 3)  # 5 nodes with 3-dimensional features

# Example placeholders for training and test data
train_data = node_features # In this simple example, we use the same data for training
test_data = rand(5, 3)     # Random test data

# Define a GNN model
gnn = Chain(
    GCNConv(3=>4, relu),  # Convolution layer, 3 input dims to 4 output dims
    GCNConv(4=>4, relu),  # Another layer
    Dense(4=>2),            # Dense layer to get 2-dimensional embeddings
    logsoftmax              # Loss function: log_softmax
)

# Convert graph to adjacency matrix
adj_matrix = Graphs.adjacency_matrix(g)
adj_matrix_sparse = sparse(adj_matrix)

# Helper to get embeddings from the model
function get_embeddings(model, features, adj)
    return model((features, adj))
end

# Define a loss function
loss(x, adj, y) = Flux.Losses.mse(gnn((x, adj)), y)

# Training the model
function train!(model, data, adj, epochs)
    opt = ADAM()
    for epoch in 1:epochs
        gs = gradient(() -> loss(data, adj, data), params(model))
        Flux.Optimise.update!(opt, params(model), gs)
    end
end

# Train the GNN model
train!(gnn, train_data, adj_matrix_sparse, 100)

# Get node embeddings
node_embeddings = get_embeddings(gnn, node_features, adj_matrix_sparse)

# Test the model with test data
test_embeddings = get_embeddings(gnn, test_data, adj_matrix_sparse)

# Compute similarity (cosine similarity in this example)
function cosine_similarity(x, y)
    dot(x, y) / (norm(x) * norm(y))
end

# Example to calculate similarity between embeddings of nodes 1 and 2
similarity = cosine_similarity(test_embeddings[:, 1], test_embeddings[:, 2])

println("Node embeddings: ", node_embeddings)
println("Test embeddings: ", test_embeddings)
println("Similarity between embeddings of nodes 1 and 2: ", similarity)