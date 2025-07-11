# Draw the graph for AIChE
using Graphs
using GraphPlot
using Colors
using Plots

out_dir = "G:\\My Drive\\Research\\GNN projects\\Figures\\"
node_count = 19
# Nodes: 1: Tin1, 2:Tin2, 3:Tin3, 4:T1, 5:T2, 6:T3, 7:x1, 8:x2, 9:x3, 10:xtot, 11:F1, 12:F2, 13:F3, 14:Q1, 15:Q2, 16:Q3, 17:Tset1, 18:Tset2, 19:Tset3
node_labels = ["Tin1", "Tin2", "Tin3", "T1", "T2", "T3", "x1", "x2", "x3", "xtot", "F1", "F2", "F3", "Q1", "Q2", "Q3", "Tset1", "Tset2", "Tset3"]
# node_labels = ["1", "2", "3", "4"]
node_kinds = ["R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3", "G", "R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3"]
self_defined_edges = [
    (7, 11),  # x1 - F1
    (7, 10),  # x1 - xtot
    (7, 14),  # x1 - Q1
    (7, 1),  # x1 - Tin1
    (7, 4),  # x1 - T1
    (14, 4),  # Q1 - T1
    (4, 1),   # T1 - Tin1
    (11, 10), # F1 - xtot
    (8, 10), # x2 - xtot
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


g = SimpleGraph(node_count)

# Create a MetaGraph
# mg = MetaGraph(g)
for (u, v) in self_defined_edges
    add_edge!(g, u, v)
end

# Assign colors based on node kinds
# lightblue, light green, lavender, lightpink
kind_to_color = Dict("G" => colorant"#ADD8E6", "R1" => colorant"#90EE90", "R2" => colorant"#E6E6FA", "R3" => colorant"#FFB6C1")

node_colors_original = [kind_to_color[node_kinds[i]] for i in 1:node_count];

# p_original = gplot(g, nodelabel=node_labels, nodefillc=node_colors_original)

# All identical case
node_kinds_all_same = ["R1", "R1", "R1", "R1", "R1", "R1", "R1", "R1", "R1", "G", "R1", "R1", "R1", "R1", "R1", "R1", "R1", "R1", "R1"]
node_colors_all_same = [kind_to_color[node_kinds_all_same[i]] for i in 1:node_count];
p_all_same = gplot(g, nodelabel=node_labels, nodefillc=node_colors_all_same)
save_filename_all_same = out_dir*"all_same_graph.png"

# T1 setpoint tracking case, in which R2=R3
node_kinds_T1_setpoint_tracking = ["R1", "R2", "R2", "R1", "R2", "R2", "R1", "R2", "R2", "G", "R1", "R2", "R2", "R1", "R2", "R2", "R1", "R2", "R2"]
node_colors_T1_setpoint_tracking = [kind_to_color[node_kinds_T1_setpoint_tracking[i]] for i in 1:node_count];
p_T1_setpoint_tracking = gplot(g, nodelabel=node_labels, nodefillc=node_colors_T1_setpoint_tracking)
save_filename_T1_setpoint_tracking = out_dir*"T1_setpoint_tracking_graph.png"

# R2 Disturbed, in which R1=R3
node_kinds_R2_Disturbed = ["R1", "R2", "R1", "R1", "R2", "R1", "R1", "R2", "R1", "G", "R1", "R2", "R1", "R1", "R2", "R1", "R1", "R2", "R1"]
node_colors_R2_Disturbed = [kind_to_color[node_kinds_R2_Disturbed[i]] for i in 1:node_count];
p_R2_Disturbed = gplot(g, nodelabel=node_labels, nodefillc=node_colors_R2_Disturbed)
save_filename_R2_Disturbed = out_dir*"R2_Disturbed_graph.png"

# R3 Disturbed, in which R1=R2
node_kinds_R2_Disturbed = ["R1", "R1", "R3", "R1", "R1", "R3", "R1", "R1", "R3", "G", "R1", "R1", "R3", "R1", "R1", "R3", "R1", "R1", "R3"]
node_colors_R2_Disturbed = [kind_to_color[node_kinds_R2_Disturbed[i]] for i in 1:node_count];
p_R2_Disturbed = gplot(g, nodelabel=node_labels, nodefillc=node_colors_R2_Disturbed)
save_filename_R2_Disturbed = out_dir*"R3_Disturbed_graph.png"

