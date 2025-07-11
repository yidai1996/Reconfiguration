using Plots, Printf, CSV
using Plots.PlotMeasures, DataFrames
function MyPlots()
    # xf = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Training dataset.xlsx") # Feel free to try any files in your computer
    xf = XLSX.readxlsx("C:\\Users\\sfay\\Documents\\Outputs\\Initial Condition Permutations\\top_initial_conditions.xlsx") # Feel free to try any files in your computer
    # sh1 = xf["Setpoint Tracking"]
    # sh2 = xf["Disturbance Rejection"]
    sh1 = xf["Sheet1"]
    num_rows, num_cols = size(sh1[:])
    println(string(num_rows) * " " * string(num_cols))
    range = "G3:G" * string(num_rows)
    println(range)
    s = convert(Array{Float64,2},sh1["G3:G32"]) # Here is the point: instead of assigning row number (32), let the code search for the number of the last row
    percentage1 = 100*convert(Array{Float64,2},sh1["L3:L32"])
    d = convert(Array{Float64,2},sh2["G6:G42"])
    percentage2 = 100*convert(Array{Float64,2},sh2["L6:L42"])
    
    p1=plot(s,percentage1,markershape = :circle,linestyle = :solid,xlabel="Set point change on xB",ylabel="C ISE percentage change (%)",framestyle=:box,label=false,xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,size=[1100,400],right_margin=0mm,left_margin=10mm,bottom_margin=9mm)
    plot!([0,0.3],[0,0],linestyle=:dot,width=2,label=false)
    display(p1)
    p2=plot(d,percentage2,markershape = :circle,linestyle = :solid,xlabel="Disturbance change on input temperature",ylabel="C ISE percentage change (%)",framestyle=:box,label=false,xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,size=[1100,400],right_margin=0mm,left_margin=10mm,bottom_margin=9mm)
    plot!([3,40],[0,0],linestyle=:dot,width=2,label=false)
    display(p2)
end

# FOPAM Poster
# AIChE PRE
function Poster_setpointtracking()
    ml = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = b[!,"xBtinitial"]

    xB_ml = zeros(3,length(tt))
    xB_p = zeros(3,length(tt))
    xB_h = zeros(3,length(tt))
    xB_m = zeros(3,length(tt))
    xB_s = zeros(3,length(tt))
    xB_b = zeros(3,length(tt))
    xB_ml[1,:] = ml[!,"xB1initial"]
    xB_ml[2,:] = ml[!,"xB2initial"]
    xB_ml[3,:] = ml[!,"xB3initial"]
    xB_p[1,:] = p[!,"xB1initial"]
    xB_p[2,:] = p[!,"xB2initial"]
    xB_p[3,:] = p[!,"xB3initial"]
    xB_h[1,:] = h[!,"xB1initial"]
    xB_h[2,:] = h[!,"xB2initial"]
    xB_h[3,:] = h[!,"xB3initial"]
    xB_m[1,:] = m[!,"xB1initial"]
    xB_m[2,:] = m[!,"xB2initial"]
    xB_m[3,:] = m[!,"xB3initial"]
    xB_s[1,:] = s[!,"xB1initial"]
    xB_s[2,:] = s[!,"xB2initial"]
    xB_s[3,:] = s[!,"xB3initial"]
    xB_b[1,:] = b[!,"xB1initial"]
    xB_b[2,:] = b[!,"xB2initial"]
    xB_b[3,:] = b[!,"xB3initial"]

    T_ml = zeros(3,length(tt))
    T_p = zeros(3,length(tt))
    T_h = zeros(3,length(tt))
    T_m = zeros(3,length(tt))
    T_s = zeros(3,length(tt))
    T_b = zeros(3,length(tt))
    T_ml[1,:] = ml[!,"T1initial"]
    T_ml[2,:] = ml[!,"T2initial"]
    T_ml[3,:] = ml[!,"T3initial"]
    T_p[1,:] = p[!,"T1initial"]
    T_p[2,:] = p[!,"T2initial"]
    T_p[3,:] = p[!,"T3initial"]
    T_h[1,:] = h[!,"T1initial"]
    T_h[2,:] = h[!,"T2initial"]
    T_h[3,:] = h[!,"T3initial"]
    T_m[1,:] = m[!,"T1initial"]
    T_m[2,:] = m[!,"T2initial"]
    T_m[3,:] = m[!,"T3initial"]
    T_s[1,:] = s[!,"T1initial"]
    T_s[2,:] = s[!,"T2initial"]
    T_s[3,:] = s[!,"T3initial"]
    T_b[1,:] = b[!,"T1initial"]
    T_b[2,:] = b[!,"T2initial"]
    T_b[3,:] = b[!,"T3initial"]

    F_ml = zeros(3,length(tt))
    F_p = zeros(3,length(tt))
    F_h = zeros(3,length(tt))
    F_m = zeros(3,length(tt))
    F_s = zeros(3,length(tt))
    F_b = zeros(3,length(tt))
    F_ml[1,:] = ml[!,"flowvt1"]
    F_ml[2,:] = ml[!,"flowvt2"]
    F_ml[3,:] = ml[!,"flowvt3"]
    F_p[1,:] = p[!,"flowvt1"]
    F_p[2,:] = p[!,"flowvt2"]
    F_p[3,:] = p[!,"flowvt3"]
    F_h[1,:] = h[!,"flowvt1"]
    F_h[2,:] = h[!,"flowvt2"]
    F_h[3,:] = h[!,"flowvt3"]
    F_m[1,:] = m[!,"flowvt1"]
    F_m[2,:] = m[!,"flowvt2"]
    F_m[3,:] = m[!,"flowvt3"]
    F_s[1,:] = s[!,"flowvt1"]
    F_s[2,:] = s[!,"flowvt2"]
    F_s[3,:] = s[!,"flowvt3"]
    F_b[1,:] = b[!,"flowvt1"]
    F_b[2,:] = b[!,"flowvt2"]
    F_b[3,:] = b[!,"flowvt3"]

    Q_ml = zeros(3,length(tt))
    Q_p = zeros(3,length(tt))
    Q_h = zeros(3,length(tt))
    Q_m = zeros(3,length(tt))
    Q_s = zeros(3,length(tt))
    Q_b = zeros(3,length(tt))
    Q_ml[1,:] = ml[!,"heatvt1"]
    Q_ml[2,:] = ml[!,"heatvt2"]
    Q_ml[3,:] = ml[!,"heatvt3"]
    Q_p[1,:] = p[!,"heatvt1"]
    Q_p[2,:] = p[!,"heatvt2"]
    Q_p[3,:] = p[!,"heatvt3"]
    Q_h[1,:] = h[!,"heatvt1"]
    Q_h[2,:] = h[!,"heatvt2"]
    Q_h[3,:] = h[!,"heatvt3"]
    Q_m[1,:] = m[!,"heatvt1"]
    Q_m[2,:] = m[!,"heatvt2"]
    Q_m[3,:] = m[!,"heatvt3"]
    Q_s[1,:] = s[!,"heatvt1"]
    Q_s[2,:] = s[!,"heatvt2"]
    Q_s[3,:] = s[!,"heatvt3"]
    Q_b[1,:] = b[!,"heatvt1"]
    Q_b[2,:] = b[!,"heatvt2"]
    Q_b[3,:] = b[!,"heatvt3"]

    label = reshape(["R$i" for i in 1:3],1,3)
    
    # plot(tt[1:9],[xBt_m[1:9],xBt_ml[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["mixing" "ML configuration switch" "series"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.pdf")
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.png")

    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid" "mixing" "sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # plot!([tt[begin],90],[0.2,0.2],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    # plot!([90,90],[0.2,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)
    # plot!([90,tt[9]],[0.16,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)

    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.pdf")
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.png")
    p_xB_ml = plot(tt, transpose(xB_ml), label=label, title="Config. Switched\n\n", titlefont=30)
    p_xB_p = plot(tt, transpose(xB_p), label=label, title="Config. Parallel\n\n", titlefont=30)
    p_xB_h = plot(tt, transpose(xB_h), label=label, title="Config. Hybrid\n\n", titlefont=30)
    p_xB_m = plot(tt, transpose(xB_m), label=label, title="Config. Mixing\n\n", titlefont=30)
    p_xB_s = plot(tt, transpose(xB_s), label=label, title="Config. Series\n\n", titlefont=30)
    p_xB_b = plot(tt, transpose(xB_b), label=label,ylabel="Mol Frac. xB (Local)", title="Continous Flow Vars\n\n", titlefont=30)
    p_T_ml = plot(tt,transpose(T_ml),label=label)
    p_T_p = plot(tt,transpose(T_p),label=label)
    p_T_h = plot(tt,transpose(T_h),label=label)
    p_T_m = plot(tt,transpose(T_m),label=label)
    p_T_s = plot(tt,transpose(T_s),label=label)
    p_T_b = plot(tt,transpose(T_b),label=label, ylabel="Temp (K) (Local)")
    p_F_ml = plot(tt,transpose(F_ml)*1e3, label=label)
    p_F_p = plot(tt,transpose(F_p)*1e3, label=label)
    p_F_h = plot(tt,transpose(F_h)*1e3, label=label)
    p_F_m = plot(tt,transpose(F_m)*1e3, label=label)
    p_F_s = plot(tt,transpose(F_s)*1e3, label=label)
    p_F_b = plot(tt,transpose(F_b)*1e3, label=label, ylabel="Flow rate (L/s)")
    p_Q_ml = plot(tt,transpose(Q_ml), label=label,xlabel="Time (s) " )
    p_Q_p = plot(tt,transpose(Q_p), label=label, xlabel="Time (s) ") 
    p_Q_h = plot(tt,transpose(Q_h), label=label, xlabel="Time (s) ")
    p_Q_m = plot(tt,transpose(Q_m), label=label, xlabel="Time (s) ")
    p_Q_s = plot(tt,transpose(Q_s), label=label, xlabel="Time (s) ")
    p_Q_b = plot(tt,transpose(Q_b), label=label, xlabel="Time (s) ", ylabel="Heating Rate (kW)")
    p_all = plot(p_xB_b, p_xB_ml, p_xB_p, p_xB_h, p_xB_m, p_xB_s, p_T_b, p_T_ml, p_T_p, p_T_h, p_T_m, p_T_s, p_F_b, p_F_ml, p_F_p, p_F_h, p_F_m, p_F_s, p_Q_b, p_Q_ml, p_Q_p, p_Q_h, p_Q_m, p_Q_s,layout=(4,6),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=26,yguidefontsize=26,framestyle=:box,size=[3000,1600],right_margin=10mm,left_margin=18mm,linewidth=3, bottom_margin=10mm)
    display(p_all)
    # p_all=plot(p_xB_ml, p_T_ml, p_F_ml, p_Q_ml, p_xB_p, p_T_p, p_F_p, p_Q_p, p_xB_h, p_T_h, p_F_h, p_Q_h, p_xB_m, p_T_m, p_F_m, p_Q_m, p_xB_s, p_T_s, p_F_h, p_Q_s,layout=(5,4),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[1500,2500],right_margin=10mm,left_margin=18mm,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_setpointtracking.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_setpointtracking.png")
    label2 = ["Continous Flow Vars" "Config. Switched" "Config. Parallel" "Config. Hybrid" "Config. Mixing" "Config. Series"]
    plot(tt[2:end],[PI_b[2:end],PI_ml[2:end],PI_p[2:end], PI_h[2:end], PI_m[2:end], PI_s[2:end]], xlabel="Times (s)", ylabel = "Performance Index of xBt", label=label2,framestyle=:box,legend=:topright,legendfont=12,xtickfontsize=8,ytickfontsize=8,xguidefontsize=12,yguidefontsize=12,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\PIcomparison_setpointtracking.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\PIcomparison_setpoingtrackiing.png")
end

function Poster_DisturbanceRejection()
    ml = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = b[!,"xBtinitial"]

    xB_ml = zeros(3,length(tt))
    xB_p = zeros(3,length(tt))
    xB_h = zeros(3,length(tt))
    xB_m = zeros(3,length(tt))
    xB_s = zeros(3,length(tt))
    xB_b = zeros(3,length(tt))
    xB_ml[1,:] = ml[!,"xB1initial"]
    xB_ml[2,:] = ml[!,"xB2initial"]
    xB_ml[3,:] = ml[!,"xB3initial"]
    xB_p[1,:] = p[!,"xB1initial"]
    xB_p[2,:] = p[!,"xB2initial"]
    xB_p[3,:] = p[!,"xB3initial"]
    xB_h[1,:] = h[!,"xB1initial"]
    xB_h[2,:] = h[!,"xB2initial"]
    xB_h[3,:] = h[!,"xB3initial"]
    xB_m[1,:] = m[!,"xB1initial"]
    xB_m[2,:] = m[!,"xB2initial"]
    xB_m[3,:] = m[!,"xB3initial"]
    xB_s[1,:] = s[!,"xB1initial"]
    xB_s[2,:] = s[!,"xB2initial"]
    xB_s[3,:] = s[!,"xB3initial"]
    xB_b[1,:] = b[!,"xB1initial"]
    xB_b[2,:] = b[!,"xB2initial"]
    xB_b[3,:] = b[!,"xB3initial"]

    T_ml = zeros(3,length(tt))
    T_p = zeros(3,length(tt))
    T_h = zeros(3,length(tt))
    T_m = zeros(3,length(tt))
    T_s = zeros(3,length(tt))
    T_b = zeros(3,length(tt))
    T_ml[1,:] = ml[!,"T1initial"]
    T_ml[2,:] = ml[!,"T2initial"]
    T_ml[3,:] = ml[!,"T3initial"]
    T_p[1,:] = p[!,"T1initial"]
    T_p[2,:] = p[!,"T2initial"]
    T_p[3,:] = p[!,"T3initial"]
    T_h[1,:] = h[!,"T1initial"]
    T_h[2,:] = h[!,"T2initial"]
    T_h[3,:] = h[!,"T3initial"]
    T_m[1,:] = m[!,"T1initial"]
    T_m[2,:] = m[!,"T2initial"]
    T_m[3,:] = m[!,"T3initial"]
    T_s[1,:] = s[!,"T1initial"]
    T_s[2,:] = s[!,"T2initial"]
    T_s[3,:] = s[!,"T3initial"]
    T_b[1,:] = b[!,"T1initial"]
    T_b[2,:] = b[!,"T2initial"]
    T_b[3,:] = b[!,"T3initial"]

    F_ml = zeros(3,length(tt))
    F_p = zeros(3,length(tt))
    F_h = zeros(3,length(tt))
    F_m = zeros(3,length(tt))
    F_s = zeros(3,length(tt))
    F_b = zeros(3,length(tt))
    F_ml[1,:] = ml[!,"flowvt1"]
    F_ml[2,:] = ml[!,"flowvt2"]
    F_ml[3,:] = ml[!,"flowvt3"]
    F_p[1,:] = p[!,"flowvt1"]
    F_p[2,:] = p[!,"flowvt2"]
    F_p[3,:] = p[!,"flowvt3"]
    F_h[1,:] = h[!,"flowvt1"]
    F_h[2,:] = h[!,"flowvt2"]
    F_h[3,:] = h[!,"flowvt3"]
    F_m[1,:] = m[!,"flowvt1"]
    F_m[2,:] = m[!,"flowvt2"]
    F_m[3,:] = m[!,"flowvt3"]
    F_s[1,:] = s[!,"flowvt1"]
    F_s[2,:] = s[!,"flowvt2"]
    F_s[3,:] = s[!,"flowvt3"]
    F_b[1,:] = b[!,"flowvt1"]
    F_b[2,:] = b[!,"flowvt2"]
    F_b[3,:] = b[!,"flowvt3"]

    Q_ml = zeros(3,length(tt))
    Q_p = zeros(3,length(tt))
    Q_h = zeros(3,length(tt))
    Q_m = zeros(3,length(tt))
    Q_s = zeros(3,length(tt))
    Q_b = zeros(3,length(tt))
    Q_ml[1,:] = ml[!,"heatvt1"]
    Q_ml[2,:] = ml[!,"heatvt2"]
    Q_ml[3,:] = ml[!,"heatvt3"]
    Q_p[1,:] = p[!,"heatvt1"]
    Q_p[2,:] = p[!,"heatvt2"]
    Q_p[3,:] = p[!,"heatvt3"]
    Q_h[1,:] = h[!,"heatvt1"]
    Q_h[2,:] = h[!,"heatvt2"]
    Q_h[3,:] = h[!,"heatvt3"]
    Q_m[1,:] = m[!,"heatvt1"]
    Q_m[2,:] = m[!,"heatvt2"]
    Q_m[3,:] = m[!,"heatvt3"]
    Q_s[1,:] = s[!,"heatvt1"]
    Q_s[2,:] = s[!,"heatvt2"]
    Q_s[3,:] = s[!,"heatvt3"]
    Q_b[1,:] = s[!,"heatvt1"]
    Q_b[2,:] = s[!,"heatvt2"]
    Q_b[3,:] = s[!,"heatvt3"]

    label = reshape(["R$i" for i in 1:3],1,3)
    
    # plot(tt[1:9],[xBt_m[1:9],xBt_ml[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["mixing" "ML configuration switch" "series"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.pdf")
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.png")

    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid" "mixing" "sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # plot!([tt[begin],90],[0.2,0.2],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    # plot!([90,90],[0.2,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)
    # plot!([90,tt[9]],[0.16,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)

    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.pdf")
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.png")
    p_xB_ml = plot(tt, transpose(xB_ml), label=label, title="Config. Switched\n\n", titlefont=30)
    p_xB_p = plot(tt, transpose(xB_p), label=label, title="Config. Parallel\n\n", titlefont=30)
    p_xB_h = plot(tt, transpose(xB_h), label=label, title="Config. Hybrid\n\n", titlefont=30)
    p_xB_m = plot(tt, transpose(xB_m), label=label, title="Config. Mixing\n\n", titlefont=30)
    p_xB_s = plot(tt, transpose(xB_s), label=label, title="Config. Series\n\n", titlefont=30)
    p_xB_b = plot(tt, transpose(xB_b), label=label,ylabel="Mol Frac. xB (Local)", title="Continous Flow Vars\n\n", titlefont=30)
    p_T_ml = plot(tt,transpose(T_ml),label=label)
    p_T_p = plot(tt,transpose(T_p),label=label)
    p_T_h = plot(tt,transpose(T_h),label=label)
    p_T_m = plot(tt,transpose(T_m),label=label)
    p_T_s = plot(tt,transpose(T_s),label=label)
    p_T_b = plot(tt,transpose(T_b),label=label, ylabel="Temp (K) (Local)")
    p_F_ml = plot(tt,transpose(F_ml)*1e3, label=label)
    p_F_p = plot(tt,transpose(F_p)*1e3, label=label)
    p_F_h = plot(tt,transpose(F_h)*1e3, label=label)
    p_F_m = plot(tt,transpose(F_m)*1e3, label=label)
    p_F_s = plot(tt,transpose(F_s)*1e3, label=label)
    p_F_b = plot(tt,transpose(F_b)*1e3, label=label, ylabel="Flow rate (L/s)")
    p_Q_ml = plot(tt,transpose(Q_ml), label=label,xlabel="Time (s) " )
    p_Q_p = plot(tt,transpose(Q_p), label=label, xlabel="Time (s) ") 
    p_Q_h = plot(tt,transpose(Q_h), label=label, xlabel="Time (s) ")
    p_Q_m = plot(tt,transpose(Q_m), label=label, xlabel="Time (s) ")
    p_Q_s = plot(tt,transpose(Q_s), label=label, xlabel="Time (s) ")
    p_Q_b = plot(tt,transpose(Q_b), label=label, xlabel="Time (s) ", ylabel="Heating Rate (kW)")
    p_all = plot(p_xB_b, p_xB_ml, p_xB_p, p_xB_h, p_xB_m, p_xB_s, p_T_b, p_T_ml, p_T_p, p_T_h, p_T_m, p_T_s, p_F_b, p_F_ml, p_F_p, p_F_h, p_F_m, p_F_s, p_Q_b, p_Q_ml, p_Q_p, p_Q_h, p_Q_m, p_Q_s,layout=(4,6),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=26,yguidefontsize=26,framestyle=:box,size=[3000,1600],right_margin=10mm,left_margin=18mm,linewidth=3, bottom_margin=10mm)
    display(p_all)
    # p_all=plot(p_xB_ml, p_T_ml, p_F_ml, p_Q_ml, p_xB_p, p_T_p, p_F_p, p_Q_p, p_xB_h, p_T_h, p_F_h, p_Q_h, p_xB_m, p_T_m, p_F_m, p_Q_m, p_xB_s, p_T_s, p_F_h, p_Q_s,layout=(5,4),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[1500,2500],right_margin=10mm,left_margin=18mm,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_disturbrejection.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_disturbrejection.png")
    label2 = ["Continous Flow Vars" "Config. Switched" "Config. Parallel" "Config. Hybrid" "Config. Mixing" "Config. Series"]
    plot(tt[2:end],[PI_b[2:end],PI_ml[2:end],PI_p[2:end], PI_h[2:end], PI_m[2:end], PI_s[2:end]], xlabel="Times (s)", ylabel = "Performance Index of xBt", label=label2,framestyle=:box,legend=:topright,legendfont=12,xtickfontsize=8,ytickfontsize=8,xguidefontsize=12,yguidefontsize=12,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\PIcomparison_disturbrejection.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\PIcomparison_disturbrejection.png")
end

# Poster_DisturbanceRejection()
# Poster_setpointtracking()
function Plot_for_MLvsNonML_setpointtracking()
    ml = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    # PI_ml_ab = ml_ab[!, "xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    # xBt_ml_ab = ml_ab[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = s[!,"xBtinitial"]

    plot([tt[1],tt[2]],[0.20,0.20],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot([tt[2],tt[2]],[0.20,0.16],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot!([tt[2],tt[9]],[0.16,0.16],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    # plot!(tt[1:9],[xBt_ml_ab[1:9], xBt_ml[1:9]],title="Mole Frac. of Product B",label=["ML-AdaBoost MPC" "ML-KNN MPC"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    plot!(tt[1:9],[xBt_b[1:9], xBt_ml[1:9],xBt_p[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)


    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_xBtComparison.png")
    plot(tt[2:9],[PI_b[2:9], PI_ml[2:9],PI_p[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    # plot(tt[2:9],[PI_ml_ab[2:9], PI_ml[2:9]],title="PI of xBt",label=["ML-AdaBoost MPC" "ML-KNN MPC" ],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_PI_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_PI_xBtComparison.png")
end

function Plot_for_MLvsNonML_disturbance()
    ml = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    # PI_ml_ab = ml_ab[!, "xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    # xBt_ml_ab = ml_ab[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = b[!,"xBtinitial"]
    plot([tt[1],tt[9]],[0.11,0.11],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    # plot!(tt[1:9],[xBt_ml_ab[1:9], xBt_ml[1:9]],title="Mole Frac. of Product B",label=["ML-AdaBoost MPC" "ML-KNN MPC"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    plot!(tt[1:9],[xBt_b[1:9], xBt_ml[1:9],xBt_p[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)


    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\DisturbRejection_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\DisturbRejection_xBtComparison.png")
    # plot(tt[2:9],[PI_b[2:9], PI_ml[2:9],PI_p[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    plot(tt[2:9],[PI_b[2:9], PI_ml[2:9],PI_p[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    # plot(tt[2:9],[PI_ml_ab[2:9], PI_ml[2:9]],title="PI of xBt",label=["ML-AdaBoost MPC" "ML-KNN MPC" ],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\DisturbRejection_PI_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\DisturbRejection_PI_xBtComparison.png")
end
# Plot_for_MLvsNonML_setpointtracking()
function AIChE_TALK()
    ml = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    tt = b[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = b[!,"xBtinitial"]

    xB_ml = zeros(3,length(tt))
    xB_p = zeros(3,length(tt))
    xB_h = zeros(3,length(tt))
    xB_m = zeros(3,length(tt))
    xB_s = zeros(3,length(tt))
    xB_b = zeros(3,length(tt))
    xB_ml[1,:] = ml[!,"xB1initial"]
    xB_ml[2,:] = ml[!,"xB2initial"]
    xB_ml[3,:] = ml[!,"xB3initial"]
    xB_p[1,:] = p[!,"xB1initial"]
    xB_p[2,:] = p[!,"xB2initial"]
    xB_p[3,:] = p[!,"xB3initial"]
    xB_h[1,:] = h[!,"xB1initial"]
    xB_h[2,:] = h[!,"xB2initial"]
    xB_h[3,:] = h[!,"xB3initial"]
    xB_m[1,:] = m[!,"xB1initial"]
    xB_m[2,:] = m[!,"xB2initial"]
    xB_m[3,:] = m[!,"xB3initial"]
    xB_s[1,:] = s[!,"xB1initial"]
    xB_s[2,:] = s[!,"xB2initial"]
    xB_s[3,:] = s[!,"xB3initial"]
    xB_b[1,:] = b[!,"xB1initial"]
    xB_b[2,:] = b[!,"xB2initial"]
    xB_b[3,:] = b[!,"xB3initial"]

    T_ml = zeros(3,length(tt))
    T_p = zeros(3,length(tt))
    T_h = zeros(3,length(tt))
    T_m = zeros(3,length(tt))
    T_s = zeros(3,length(tt))
    T_b = zeros(3,length(tt))
    T_ml[1,:] = ml[!,"T1initial"]
    T_ml[2,:] = ml[!,"T2initial"]
    T_ml[3,:] = ml[!,"T3initial"]
    T_p[1,:] = p[!,"T1initial"]
    T_p[2,:] = p[!,"T2initial"]
    T_p[3,:] = p[!,"T3initial"]
    T_h[1,:] = h[!,"T1initial"]
    T_h[2,:] = h[!,"T2initial"]
    T_h[3,:] = h[!,"T3initial"]
    T_m[1,:] = m[!,"T1initial"]
    T_m[2,:] = m[!,"T2initial"]
    T_m[3,:] = m[!,"T3initial"]
    T_s[1,:] = s[!,"T1initial"]
    T_s[2,:] = s[!,"T2initial"]
    T_s[3,:] = s[!,"T3initial"]
    T_b[1,:] = b[!,"T1initial"]
    T_b[2,:] = b[!,"T2initial"]
    T_b[3,:] = b[!,"T3initial"]

    F_ml = zeros(3,length(tt))
    F_p = zeros(3,length(tt))
    F_h = zeros(3,length(tt))
    F_m = zeros(3,length(tt))
    F_s = zeros(3,length(tt))
    F_b = zeros(3,length(tt))
    F_ml[1,:] = ml[!,"flowvt1"]
    F_ml[2,:] = ml[!,"flowvt2"]
    F_ml[3,:] = ml[!,"flowvt3"]
    F_p[1,:] = p[!,"flowvt1"]
    F_p[2,:] = p[!,"flowvt2"]
    F_p[3,:] = p[!,"flowvt3"]
    F_h[1,:] = h[!,"flowvt1"]
    F_h[2,:] = h[!,"flowvt2"]
    F_h[3,:] = h[!,"flowvt3"]
    F_m[1,:] = m[!,"flowvt1"]
    F_m[2,:] = m[!,"flowvt2"]
    F_m[3,:] = m[!,"flowvt3"]
    F_s[1,:] = s[!,"flowvt1"]
    F_s[2,:] = s[!,"flowvt2"]
    F_s[3,:] = s[!,"flowvt3"]
    F_b[1,:] = b[!,"flowvt1"]
    F_b[2,:] = b[!,"flowvt2"]
    F_b[3,:] = b[!,"flowvt3"]

    Q_ml = zeros(3,length(tt))
    Q_p = zeros(3,length(tt))
    Q_h = zeros(3,length(tt))
    Q_m = zeros(3,length(tt))
    Q_s = zeros(3,length(tt))
    Q_b = zeros(3,length(tt))
    Q_ml[1,:] = ml[!,"heatvt1"]
    Q_ml[2,:] = ml[!,"heatvt2"]
    Q_ml[3,:] = ml[!,"heatvt3"]
    Q_p[1,:] = p[!,"heatvt1"]
    Q_p[2,:] = p[!,"heatvt2"]
    Q_p[3,:] = p[!,"heatvt3"]
    Q_h[1,:] = h[!,"heatvt1"]
    Q_h[2,:] = h[!,"heatvt2"]
    Q_h[3,:] = h[!,"heatvt3"]
    Q_m[1,:] = m[!,"heatvt1"]
    Q_m[2,:] = m[!,"heatvt2"]
    Q_m[3,:] = m[!,"heatvt3"]
    Q_s[1,:] = s[!,"heatvt1"]
    Q_s[2,:] = s[!,"heatvt2"]
    Q_s[3,:] = s[!,"heatvt3"]
    Q_b[1,:] = s[!,"heatvt1"]
    Q_b[2,:] = s[!,"heatvt2"]
    Q_b[3,:] = s[!,"heatvt3"]

    label = reshape(["R$i" for i in 1:3],1,3)
    
    # plot(tt[1:9],[xBt_m[1:9],xBt_ml[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["mixing" "ML configuration switch" "series"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.pdf")
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.png")

    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid" "mixing" "sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    # plot!([tt[begin],90],[0.2,0.2],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    # plot!([90,90],[0.2,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)
    # plot!([90,tt[9]],[0.16,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)

    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.pdf")
    # savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.png")
    p_xB_ml = plot(tt, transpose(xB_ml), label=label, title="Config. Switched\n\n", titlefont=30)
    p_xB_p = plot(tt, transpose(xB_p), label=label, title="Config. Parallel\n\n", titlefont=30)
    p_xB_h = plot(tt, transpose(xB_h), label=label, title="Config. Hybrid\n\n", titlefont=30)
    p_xB_m = plot(tt, transpose(xB_m), label=label, title="Config. Mixing\n\n", titlefont=30)
    p_xB_s = plot(tt, transpose(xB_s), label=label, title="Config. Series\n\n", titlefont=30)
    p_xB_b = plot(tt, transpose(xB_b), label=label,ylabel="Mol Frac. xB (Local)", title="Continous Flow Vars\n\n", titlefont=30)
    p_T_ml = plot(tt,transpose(T_ml),label=label)
    p_T_p = plot(tt,transpose(T_p),label=label)
    p_T_h = plot(tt,transpose(T_h),label=label)
    p_T_m = plot(tt,transpose(T_m),label=label)
    p_T_s = plot(tt,transpose(T_s),label=label)
    p_T_b = plot(tt,transpose(T_b),label=label, ylabel="Temp (K) (Local)")
    p_F_ml = plot(tt,transpose(F_ml)*1e3, label=label)
    p_F_p = plot(tt,transpose(F_p)*1e3, label=label)
    p_F_h = plot(tt,transpose(F_h)*1e3, label=label)
    p_F_m = plot(tt,transpose(F_m)*1e3, label=label)
    p_F_s = plot(tt,transpose(F_s)*1e3, label=label)
    p_F_b = plot(tt,transpose(F_b)*1e3, label=label, ylabel="Flow rate (L/s)")
    p_Q_ml = plot(tt,transpose(Q_ml), label=label,xlabel="Time (s) " )
    p_Q_p = plot(tt,transpose(Q_p), label=label, xlabel="Time (s) ") 
    p_Q_h = plot(tt,transpose(Q_h), label=label, xlabel="Time (s) ")
    p_Q_m = plot(tt,transpose(Q_m), label=label, xlabel="Time (s) ")
    p_Q_s = plot(tt,transpose(Q_s), label=label, xlabel="Time (s) ")
    p_Q_b = plot(tt,transpose(Q_b), label=label, xlabel="Time (s) ", ylabel="Heating Rate (kW)")
    p_all = plot(p_xB_b, p_xB_ml, p_xB_p, p_xB_h, p_xB_m, p_xB_s, p_T_b, p_T_ml, p_T_p, p_T_h, p_T_m, p_T_s, p_F_b, p_F_ml, p_F_p, p_F_h, p_F_m, p_F_s, p_Q_b, p_Q_ml, p_Q_p, p_Q_h, p_Q_m, p_Q_s,layout=(4,6),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=26,yguidefontsize=26,framestyle=:box,size=[3000,1600],right_margin=10mm,left_margin=18mm,linewidth=3, bottom_margin=10mm)
    display(p_all)
    # p_all=plot(p_xB_ml, p_T_ml, p_F_ml, p_Q_ml, p_xB_p, p_T_p, p_F_p, p_Q_p, p_xB_h, p_T_h, p_F_h, p_Q_h, p_xB_m, p_T_m, p_F_m, p_Q_m, p_xB_s, p_T_s, p_F_h, p_Q_s,layout=(5,4),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[1500,2500],right_margin=10mm,left_margin=18mm,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_SPandDR.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_SPandDR.png")
    label2 = ["Continous Flow Vars" "Config. Switched" "Config. Parallel" "Config. Hybrid" "Config. Mixing" "Config. Series"]
    plot(tt[2:end],[PI_b[2:end],PI_ml[2:end],PI_p[2:end], PI_h[2:end], PI_m[2:end], PI_s[2:end]], xlabel="Times (s)", ylabel = "Performance Index of xBt", label=label2,framestyle=:box,legend=:topright,legendfont=12,xtickfontsize=8,ytickfontsize=8,xguidefontsize=12,yguidefontsize=12,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\PIcomparison_SPandDR.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\PIcomparison_SPandDR.png")
end
function Plot_for_MLvsNonML_setpointtracking()
    ml = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    # PI_ml_ab = ml_ab[!, "xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    # xBt_ml_ab = ml_ab[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = s[!,"xBtinitial"]

    plot([tt[1],tt[2]],[0.20,0.20],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot([tt[2],tt[2]],[0.20,0.16],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot!([tt[2],tt[9]],[0.16,0.16],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    # plot!(tt[1:9],[xBt_ml_ab[1:9], xBt_ml[1:9]],title="Mole Frac. of Product B",label=["ML-AdaBoost MPC" "ML-KNN MPC"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    plot!(tt[1:9],[xBt_b[1:9], xBt_ml[1:9],xBt_p[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)


    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_xBtComparison.png")
    plot(tt[2:9],[PI_b[2:9], PI_ml[2:9],PI_p[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    # plot(tt[2:9],[PI_ml_ab[2:9], PI_ml[2:9]],title="PI of xBt",label=["ML-AdaBoost MPC" "ML-KNN MPC" ],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_PI_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\Setpoint_PI_xBtComparison.png")
end

function aichetalk()
    ml = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    p = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    h = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    m = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    s = CSV.read("G:\\My Drive\\Research\\SummerResearch\\PermutationWithoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    b = CSV.read("G:\\My Drive\\Research\\MLReconfiguration\\nobinaryflow_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.19_xB2_0.19_xB3_0.19_T0_340.0SetChange_xB_-0.09.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    # PI_ml_ab = ml_ab[!, "xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]
    PI_b = b[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    # xBt_ml_ab = ml_ab[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]
    xBt_b = b[!,"xBtinitial"]
    plot([tt[1],tt[2]],[0.19,0.19],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot!([tt[2],tt[2]],[0.19,0.10],linestyle=:dot,width=1,linecolor="Crimson")
    plot!([tt[2],tt[9]],[0.1,0.1],linestyle=:dot,width=1,linecolor="Crimson")
    # plot!(tt[1:9],[xBt_ml_ab[1:9], xBt_ml[1:9]],title="Mole Frac. of Product B",label=["ML-AdaBoost MPC" "ML-KNN MPC"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    plot!(tt[1:9],[xBt_b[1:9], xBt_ml[1:9],xBt_p[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)


    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\SPandDR_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\SPandDR_xBtComparison.png")
    # plot(tt[2:9],[PI_b[2:9], PI_ml[2:9],PI_p[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    plot(tt[2:9],[PI_b[2:9], PI_ml[2:9],PI_p[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["Continous F Vars" "Switched" "Parallel"  "Hybrid" "Mixing" "Sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    # plot(tt[2:9],[PI_ml_ab[2:9], PI_ml[2:9]],title="PI of xBt",label=["ML-AdaBoost MPC" "ML-KNN MPC" ],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\SPandDR_PI_xBtComparison.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\SPandDR_PI_xBtComparison.png")


    xB_ml = zeros(3,length(tt))
    xB_p = zeros(3,length(tt))
    xB_h = zeros(3,length(tt))
    xB_m = zeros(3,length(tt))
    xB_s = zeros(3,length(tt))
    xB_b = zeros(3,length(tt))
    xB_ml[1,:] = ml[!,"xB1initial"]
    xB_ml[2,:] = ml[!,"xB2initial"]
    xB_ml[3,:] = ml[!,"xB3initial"]
    xB_p[1,:] = p[!,"xB1initial"]
    xB_p[2,:] = p[!,"xB2initial"]
    xB_p[3,:] = p[!,"xB3initial"]
    xB_h[1,:] = h[!,"xB1initial"]
    xB_h[2,:] = h[!,"xB2initial"]
    xB_h[3,:] = h[!,"xB3initial"]
    xB_m[1,:] = m[!,"xB1initial"]
    xB_m[2,:] = m[!,"xB2initial"]
    xB_m[3,:] = m[!,"xB3initial"]
    xB_s[1,:] = s[!,"xB1initial"]
    xB_s[2,:] = s[!,"xB2initial"]
    xB_s[3,:] = s[!,"xB3initial"]
    xB_b[1,:] = b[!,"xB1initial"]
    xB_b[2,:] = b[!,"xB2initial"]
    xB_b[3,:] = b[!,"xB3initial"]

    T_ml = zeros(3,length(tt))
    T_p = zeros(3,length(tt))
    T_h = zeros(3,length(tt))
    T_m = zeros(3,length(tt))
    T_s = zeros(3,length(tt))
    T_b = zeros(3,length(tt))
    T_ml[1,:] = ml[!,"T1initial"]
    T_ml[2,:] = ml[!,"T2initial"]
    T_ml[3,:] = ml[!,"T3initial"]
    T_p[1,:] = p[!,"T1initial"]
    T_p[2,:] = p[!,"T2initial"]
    T_p[3,:] = p[!,"T3initial"]
    T_h[1,:] = h[!,"T1initial"]
    T_h[2,:] = h[!,"T2initial"]
    T_h[3,:] = h[!,"T3initial"]
    T_m[1,:] = m[!,"T1initial"]
    T_m[2,:] = m[!,"T2initial"]
    T_m[3,:] = m[!,"T3initial"]
    T_s[1,:] = s[!,"T1initial"]
    T_s[2,:] = s[!,"T2initial"]
    T_s[3,:] = s[!,"T3initial"]
    T_b[1,:] = b[!,"T1initial"]
    T_b[2,:] = b[!,"T2initial"]
    T_b[3,:] = b[!,"T3initial"]

    F_ml = zeros(3,length(tt))
    F_p = zeros(3,length(tt))
    F_h = zeros(3,length(tt))
    F_m = zeros(3,length(tt))
    F_s = zeros(3,length(tt))
    F_b = zeros(3,length(tt))
    F_ml[1,:] = ml[!,"flowvt1"]
    F_ml[2,:] = ml[!,"flowvt2"]
    F_ml[3,:] = ml[!,"flowvt3"]
    F_p[1,:] = p[!,"flowvt1"]
    F_p[2,:] = p[!,"flowvt2"]
    F_p[3,:] = p[!,"flowvt3"]
    F_h[1,:] = h[!,"flowvt1"]
    F_h[2,:] = h[!,"flowvt2"]
    F_h[3,:] = h[!,"flowvt3"]
    F_m[1,:] = m[!,"flowvt1"]
    F_m[2,:] = m[!,"flowvt2"]
    F_m[3,:] = m[!,"flowvt3"]
    F_s[1,:] = s[!,"flowvt1"]
    F_s[2,:] = s[!,"flowvt2"]
    F_s[3,:] = s[!,"flowvt3"]
    F_b[1,:] = b[!,"flowvt1"]
    F_b[2,:] = b[!,"flowvt2"]
    F_b[3,:] = b[!,"flowvt3"]

    Q_ml = zeros(3,length(tt))
    Q_p = zeros(3,length(tt))
    Q_h = zeros(3,length(tt))
    Q_m = zeros(3,length(tt))
    Q_s = zeros(3,length(tt))
    Q_b = zeros(3,length(tt))
    Q_ml[1,:] = ml[!,"heatvt1"]
    Q_ml[2,:] = ml[!,"heatvt2"]
    Q_ml[3,:] = ml[!,"heatvt3"]
    Q_p[1,:] = p[!,"heatvt1"]
    Q_p[2,:] = p[!,"heatvt2"]
    Q_p[3,:] = p[!,"heatvt3"]
    Q_h[1,:] = h[!,"heatvt1"]
    Q_h[2,:] = h[!,"heatvt2"]
    Q_h[3,:] = h[!,"heatvt3"]
    Q_m[1,:] = m[!,"heatvt1"]
    Q_m[2,:] = m[!,"heatvt2"]
    Q_m[3,:] = m[!,"heatvt3"]
    Q_s[1,:] = s[!,"heatvt1"]
    Q_s[2,:] = s[!,"heatvt2"]
    Q_s[3,:] = s[!,"heatvt3"]
    Q_b[1,:] = s[!,"heatvt1"]
    Q_b[2,:] = s[!,"heatvt2"]
    Q_b[3,:] = s[!,"heatvt3"]

    label = reshape(["R$i" for i in 1:3],1,3)

    p_xB_ml = plot(tt, transpose(xB_ml), label=label, title="Config. Switched\n\n", titlefont=30)
    p_xB_p = plot(tt, transpose(xB_p), label=label, title="Config. Parallel\n\n", titlefont=30)
    p_xB_h = plot(tt, transpose(xB_h), label=label, title="Config. Hybrid\n\n", titlefont=30)
    p_xB_m = plot(tt, transpose(xB_m), label=label, title="Config. Mixing\n\n", titlefont=30)
    p_xB_s = plot(tt, transpose(xB_s), label=label, title="Config. Series\n\n", titlefont=30)
    p_xB_b = plot(tt, transpose(xB_b), label=label,ylabel="Mol Frac. xB (Local)", title="Continous Flow Vars\n\n", titlefont=30)
    p_T_ml = plot(tt,transpose(T_ml),label=label)
    p_T_p = plot(tt,transpose(T_p),label=label)
    p_T_h = plot(tt,transpose(T_h),label=label)
    p_T_m = plot(tt,transpose(T_m),label=label)
    p_T_s = plot(tt,transpose(T_s),label=label)
    p_T_b = plot(tt,transpose(T_b),label=label, ylabel="Temp (K) (Local)")
    p_F_ml = plot(tt,transpose(F_ml)*1e3, label=label)
    p_F_p = plot(tt,transpose(F_p)*1e3, label=label)
    p_F_h = plot(tt,transpose(F_h)*1e3, label=label)
    p_F_m = plot(tt,transpose(F_m)*1e3, label=label)
    p_F_s = plot(tt,transpose(F_s)*1e3, label=label)
    p_F_b = plot(tt,transpose(F_b)*1e3, label=label, ylabel="Flow rate (L/s)")
    p_Q_ml = plot(tt,transpose(Q_ml), label=label,xlabel="Time (s) " )
    p_Q_p = plot(tt,transpose(Q_p), label=label, xlabel="Time (s) ") 
    p_Q_h = plot(tt,transpose(Q_h), label=label, xlabel="Time (s) ")
    p_Q_m = plot(tt,transpose(Q_m), label=label, xlabel="Time (s) ")
    p_Q_s = plot(tt,transpose(Q_s), label=label, xlabel="Time (s) ")
    p_Q_b = plot(tt,transpose(Q_b), label=label, xlabel="Time (s) ", ylabel="Heating Rate (kW)")
    p_all = plot(p_xB_b, p_xB_ml, p_xB_p, p_xB_h, p_xB_m, p_xB_s, p_T_b, p_T_ml, p_T_p, p_T_h, p_T_m, p_T_s, p_F_b, p_F_ml, p_F_p, p_F_h, p_F_m, p_F_s, p_Q_b, p_Q_ml, p_Q_p, p_Q_h, p_Q_m, p_Q_s,layout=(4,6),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=26,yguidefontsize=26,framestyle=:box,size=[3000,1600],right_margin=10mm,left_margin=18mm,linewidth=3, bottom_margin=10mm)
    display(p_all)
    # p_all=plot(p_xB_ml, p_T_ml, p_F_ml, p_Q_ml, p_xB_p, p_T_p, p_F_p, p_Q_p, p_xB_h, p_T_h, p_F_h, p_Q_h, p_xB_m, p_T_m, p_F_m, p_Q_m, p_xB_s, p_T_s, p_F_h, p_Q_s,layout=(5,4),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[1500,2500],right_margin=10mm,left_margin=18mm,linewidth=3)
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_SPandDR.pdf")
    savefig("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\otherPerformanceComparison_SPandDR.png")
end

function aichetalk2024()
    # all same
    df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
    # combined
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_310.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_10.0.csv", DataFrame)
    # T1 setpoint
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_10.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
    # R2 disturbance
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_310.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
    # R3 disturbance
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_310.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

    tt = df[!,"times"]
    # PI_ml = ml[!,"xBt_PI"]

    # xBt = round.(df[!,"xBtinitial"],digits=3)
    xB1 = round.(df[!, "xB1initial"], digits=3)
    xB2 = round.(df[!, "xB2initial"], digits=3)
    xB3 = round.(df[!, "xB3initial"], digits=3)
    
    T1 = round.(df[!, "T1initial"], digits=3)
    T2 = round.(df[!, "T2initial"], digits=3)
    T3 = round.(df[!, "T3initial"], digits=3)

    colors = ["orange", "mediumpurple1", "coral1"]
    
    plot(tt[14:23],
        # alldifferent
        T1[14:23],
        label="R1=R2=R3",
        # all same
        # [xB1[14:23]],
        # label="R1=R2=R3" ,
        # R3 disturbance
        # xB1[14:23],
        # label="R1=R2=R3",
        xlabel="Times (s)",
        # ylabel="Individual xB",
        ylabel="Individual T (K)",
        framestyle=:box,
        linewidth=2,
        legend=:topright,
        size=(400, 300),
        ylim=(385, 395),
        lc="orange")
    # plot!(tt[14:23],
    #     # alldifferent
    #     T2[14:23], 
    #     label= "R2",
    #     # all same
    #     # [xB1[14:23]],
    #     # label="R1=R2=R3" ,
    #     # R3 disturbance
    #     # T3[14:23],
    #     # label="R3",
    #     linewidth=2,
    #     lc=:mediumpurple1)
    # plot!(tt[14:23],
    #     # alldifferent
    #     T3[14:23], 
    #     label= "R3",
    #     # all same
    #     # [xB1[14:23]],
    #     # label="R1=R2=R3" ,
    #     # R3 disturbance
    #     # T3[14:23],
    #     # label="R3",
    #     linewidth=2,
    #     lc=:coral1)
    # plot(tt[2:9],[PI_ml_ab[2:9], PI_ml[2:9]],title="PI of xBt",label=["ML-AdaBoost MPC" "ML-KNN MPC" ],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\allsame_T.pdf")
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\combine_T.pdf")
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\T1setpoint_T.pdf")
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R2disturbance_T.pdf")
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R3disturbance_T.pdf")
end
aichetalk2024()

function aichetalk2024_2()
    # all same
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
    # combined
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_310.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_10.0.csv", DataFrame)
    # T1 setpoint
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_300.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_10.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
    # R2 disturbance
    # df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\19nodes\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_310.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)
    df = CSV.read("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_Tin1_300.0_Tin2_320.0_Tin3_300.0SetChange_xB_0.0SetChange_T1_0.0SetChange_T2_0.0SetChange_T3_0.0.csv", DataFrame)

    tt = df[!,"times"]
    # PI_ml = ml[!,"xBt_PI"]

    xBt = round.(df[!,"xBtinitial"],digits=3)
    xB1 = round.(df[!, "xB1initial"], digits=3)
    xB2 = round.(df[!, "xB2initial"], digits=3)
    xB3 = round.(df[!, "xB3initial"], digits=3)
    
    T1 = round.(df[!, "T1initial"],digits=3)
    T2 = round.(df[!, "T2initial"],digits=3)
    T3 = round.(df[!, "T3initial"],digits=3)

    F1 = round.(df[!, "flowvt1"]*10e3,digits=3)
    F2 = round.(df[!, "flowvt2"]*10e3,digits=3)
    F3 = round.(df[!, "flowvt3"]*10e3,digits=3)

    Q1 = round.(df[!, "heatvt1"],digits=3)
    Q2 = round.(df[!, "heatvt2"],digits=3)
    Q3 = round.(df[!, "heatvt3"],digits=3)
    
        
    # plot(tt[14:23],
    #     # alldifferent
    #     xB1[14:23],
    #     label="R1=R3",
    #     # all same
    #     # [xB1[14:23]],
    #     # label="R1=R2=R3" ,
    #     # R3 disturbance
    #     # xB1[14:23],
    #     # label="R1=R2=R3",
    #     xlabel="Times (s)",
    #     ylabel="Individual xB",
    #     # ylabel="Flow rate (L/s)",
    #     # ylabel="Heating Rate (kW)",
    #     # ylabel="Individual T (K)",
    #     framestyle=:box,
    #     linewidth=2,
    #     legend=:topright,
    #     size=(400, 300),
    #     # ylim=(385, 395),
    #     lc="orange")
    # plot!(tt[14:23],
    #     # alldifferent
    #     xB2[14:23], 
    #     label= "R2",
    #     # all same
    #     # [xB1[14:23]],
    #     # label="R1=R2=R3" ,
    #     # R3 disturbance
    #     # T3[14:23],
    #     # label="R3",
    #     linewidth=2,
    #     lc=:mediumpurple1)

    #     savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R2_DIST_x.pdf")
    # plot!(tt[14:23],
    #     # alldifferent
    #     T3[14:23], 
    #     label= "R3",
    #     # all same
    #     # [xB1[14:23]],
    #     # label="R1=R2=R3" ,
    #     # R3 disturbance
    #     # T3[14:23],
    #     # label="R3",
    #     linewidth=2,
    #     lc=:coral1)
    
    # plot(tt[14:23],
    #     [xB1[14:23], xB2[14:23], xB3[14:23]],
    #     ylabel = "Individual xB",
    #     label=["xB1" "xB2"  "xB3" ],
    #     xlabel="Times (s)",
    #     framestyle=:box,
    #     linewidth=2,
    #     legend=:topright,
    #     size=(400, 300))
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R2_x.pdf")

    # plot(tt[14:23],
    # [T1[14:23], T2[14:23], T3[14:23]],
    # ylabel="Temperature",
    # label=["T1" "T2"  "T3" ],
    # xlabel="Times (s)",
    # framestyle=:box,
    # linewidth=2,
    # legend=:topright,
    # size=(400, 300))
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R2_T.pdf")

    # plot(tt[14:23],
    #     xBt[14:23],
    #     ylabel="xBtot",
    #     xlabel="Times (s)",
    #     framestyle=:box,
    #     linewidth=2,
    #     legend=:topright,
    #     size=(400, 300))
    # savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R2_xtot.pdf")
    # p_all=plot(p_xB_ml, p_T_ml, p_F_ml, p_Q_ml, p_xB_p, p_T_p, p_F_p, p_Q_p, p_xB_h, p_T_h, p_F_h, p_Q_h, p_xB_m, p_T_m, p_F_m, p_Q_m, p_xB_s, p_T_s, p_F_h, p_Q_s,layout=(5,4),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[1500,2500],right_margin=10mm,left_margin=18mm,linewidth=3)
    plot([tt[14],tt[16]],[1,1],linestyle=:dot,width=2,yaxis=true, label="R1=R2=R3",framestyle=:box,lc=:red, xlabel="Times (s)",legend=:topright, size=(300, 100),ylim=(0.5, 1.5))
    plot!([tt[16],tt[23]],[1,1],linestyle=:dot,width=2,label="R1=R3",lc=:blue)
    savefig("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel\\AICHE\\R2_DIST_MPC.pdf")
end
aichetalk2024_2()
