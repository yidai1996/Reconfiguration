using JuMP, CPLEX

# ------------------------------
# 1. 网络与基准参数
# ------------------------------
nodes = [1,2,3,4]
edges = [(1,2), (2,3), (3,4)]

# 基准
const Sbase = 6e6           # 6 MVA
const Vbase_HV = 12.47e3    # 12.47 kV
const Vbase_LV = 4.16e3     # 4.16 kV

# 线路长度 (mi)
L12 = 2000/5280.0           # Node1–Node2
L34 = 2500/5280.0           # Node3–Node4

# 正序阻抗 Ω/mile（文档值）
Zpos = 0.306 + im*0.6272

# 变压器标幺阻抗 (12.47/4.16, 6 MVA)
Ztr_pu = 0.01 + im*0.06
a = Vbase_HV / Vbase_LV      # 3.0

# 换算到 pu 域
r = Dict{Tuple{Int,Int},Float64}()
x = Dict{Tuple{Int,Int},Float64}()

# Line 1–2 (HV)
z12 = Zpos * L12
r[(1,2)] = real(z12) * Sbase / Vbase_HV^2
x[(1,2)] = imag(z12) * Sbase / Vbase_HV^2

# Transformer (2–3)
r_tr = real(Ztr_pu) / a^2
x_tr = imag(Ztr_pu) / a^2
r[(2,3)] = r_tr
x[(2,3)] = x_tr

# Line 3–4 (LV)
z34 = Zpos * L34
r[(3,4)] = real(z34) * Sbase / Vbase_LV^2
x[(3,4)] = imag(z34) * Sbase / Vbase_LV^2

# 负荷 (Node 4, 平衡 closed-delta)
P4_kW = 1800.0
pf = 0.9
Pload = Dict(4 => P4_kW/1e3/Sbase)             # pu
Qload = Dict(4 => P4_kW/1e3/Sbase * tan(acos(pf)))

# ------------------------------
# 2. 建立 JuMP 模型
# ------------------------------
model = Model(CPLEX.Optimizer)
set_silent(model)   # 取消求解器输出（可选）

# 变量
@variable(model, V[i in nodes] >= 0)                   # 节点电压平方 pu
@variable(model, P[e in edges])                        # 支路有功 pu
@variable(model, Q[e in edges])                        # 支路无功 pu
@variable(model, l[e in edges])                        # convex relaxation
@variable(model, Lvar[e in edges] >= 0)                # 电流平方 pu

# Slack bus
@constraint(model, V[1] == 1.0)

# DistFlow 方程
for (i,j) in edges
    # 下游支路
    downstream = filter(e->e[1]==j, edges)

    # 功率平衡
    @constraint(model,
        P[(i,j)] == sum(P[e] for e in downstream) + get(Pload,j,0.0) + r[(i,j)]*Lvar[(i,j)]
    )
    @constraint(model,
        Q[(i,j)] == sum(Q[e] for e in downstream) + get(Qload,j,0.0) + x[(i,j)]*Lvar[(i,j)]
    )

    # 压降
    @constraint(model,
        V[j] == V[i]
               - 2*(r[(i,j)]*P[(i,j)] + x[(i,j)]*Q[(i,j)])
               + (r[(i,j)]^2 + x[(i,j)]^2)*Lvar[(i,j)]
    )

    # 电流定义（convex relaxation）
    @constraint(model, (P[(i,j)]^2 + Q[(i,j)]^2 <= V[i] * Lvar[(i,j)]) )
end

# 目标：最小化线路损耗 I²R
@objective(model, Min, sum(r[e]*Lvar[e] for e in edges))

# ------------------------------
# 3. 求解并输出
# ------------------------------
optimize!(model)

println("节点电压 (pu)：")
for i in nodes
    println("  Bus $i: ", round(sqrt(value(V[i])), digits=4))
end

println("\n支路潮流 (pu)：")
for (i,j) in edges
    println("  P[$i→$j] = ", round(value(P[(i,j)]), digits=4),
            ", Q[$i→$j] = ", round(value(Q[(i,j)]), digits=4))
end
