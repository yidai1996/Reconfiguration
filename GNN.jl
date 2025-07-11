import Graphs, GraphNeuralNetworks
import GraphNeuralNetworks.GNNGraph
source = [1,1,2,2,3,3,3,4];
target = [2,3,1,3,1,2,4,3];
g1 = GNNGraph(source, target)

adjlist = [[2,3], [1,3], [1,2,4], [3]]
g2 = GNNGraph(adjlist)

# Construct a GNNGraph from from a Graphs.jl's graph
# lg = erdos_renyi(10, 30)
# g = GNNGraph(lg)

# use Graphs.jl to construction the undirected graph of the three-reactor system
# assumption: every variable (x_{i,t}) is a node, so does everyone has a relationship with others? Or some nodes may not have direct relationship, e.g. mole fraction at time = 1, and the flow rate at time = 12


# Test of demo code
using LightGraphs
using Flux
using GraphNeuralNetworks
using SparseArrays
using GraphPlot # https://github.com/JuliaGraphs/GraphPlot.jl
using Graphs

# Step 1: Define the graph
g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)

# Step 2: Define node features
input_dim = 3
num_nodes = nv(g)
node_features = rand(Float32, input_dim, num_nodes)

# Step 3: Convert the graph to an adjacency matrix
adjacency_matrix = LightGraphs.adjacency_matrix(g)

# Step 4: Define the GNN model
model = Chain(
    GCNConv(input_dim => 4, relu),
    GCNConv(4 => 2)
)

# Step 5: Forward pass to get initial embeddings with the correct syntax
node_embeddings = model(g, adjacency_matrix)

# (Optional) Step 6: Training the GNN
labels = rand(Bool, num_nodes)
loss_fn = (model, x, g, y) -> Flux.logitcrossentropy(model((x, g)), y)
opt = Flux.Adam()

epochs = 100
for epoch in 1:epochs
    Flux.train!(loss_fn, params(model), [(node_features, adjacency_matrix, labels)], opt)
    println("Epoch $epoch complete")
end

# Step 7: Obtain final node embeddings
final_embeddings = model((node_features, adjacency_matrix))
println("Final node embeddings: ", final_embeddings)


# Demo code from GNN.jl website: https://carlolucibello.github.io/GraphNeuralNetworks.jl/dev/models/
din, d, dout = 3, 4, 2
chain = GNNChain(GCNConv(din => d, relu),
                 GCNConv(d => dout))


g = rand_graph(10, 30)
# g should be replaced by our configurations 
nodelabel = 1:Graphs.nv(g)
gplot(g, nodelabel=nodelabel)
model = WithGraph(chain, g)

X = randn(Float32, din, 10)

# Pass only X as input, the model already contains the graph.
y = model(X) 
