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
function Poster()
    ml = CSV.read("C:\\Users\\10060\\Downloads\\Research\\posterdata\\withML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    p = CSV.read("C:\\Users\\10060\\Downloads\\Research\\posterdata\\withoutML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    h = CSV.read("C:\\Users\\10060\\Downloads\\Research\\posterdata\\withoutML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    m = CSV.read("C:\\Users\\10060\\Downloads\\Research\\posterdata\\withoutML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    s = CSV.read("C:\\Users\\10060\\Downloads\\Research\\posterdata\\withoutML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.2_xB2_0.2_xB3_0.2_T0_300.0SetChange_xB_-0.04.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]

    xB_ml = zeros(3,length(tt))
    xB_p = zeros(3,length(tt))
    xB_h = zeros(3,length(tt))
    xB_m = zeros(3,length(tt))
    xB_s = zeros(3,length(tt))
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

    T_ml = zeros(3,length(tt))
    T_p = zeros(3,length(tt))
    T_h = zeros(3,length(tt))
    T_m = zeros(3,length(tt))
    T_s = zeros(3,length(tt))
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

    F_ml = zeros(3,length(tt))
    F_p = zeros(3,length(tt))
    F_h = zeros(3,length(tt))
    F_m = zeros(3,length(tt))
    F_s = zeros(3,length(tt))
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

    Q_ml = zeros(3,length(tt))
    Q_p = zeros(3,length(tt))
    Q_h = zeros(3,length(tt))
    Q_m = zeros(3,length(tt))
    Q_s = zeros(3,length(tt))
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

    label = reshape(["R$i" for i in 1:3],1,3)
    
    plot(tt[1:9],[xBt_m[1:9],xBt_ml[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["mixing" "ML configuration switch" "series"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.pdf")
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison_mixingseries.png")

    # plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    plot(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid" "mixing" "sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)
    plot!([tt[begin],90],[0.2,0.2],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot!([90,90],[0.2,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)
    plot!([90,tt[9]],[0.16,0.16],linestyle=:dot,width=1,linecolor="Crimson",label=false)

    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.pdf")
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\xBtcomparison.png")
    p_xB_ml = plot(tt, transpose(xB_ml), label=label,ylabel="Mol Frac. xB (Local)")
    p_xB_p = plot(tt, transpose(xB_p), label=label)
    p_xB_h = plot(tt, transpose(xB_h), label=label)
    p_xB_m = plot(tt, transpose(xB_m), label=label)
    p_xB_s = plot(tt, transpose(xB_s), label=label)
    p_T_ml = plot(tt,transpose(T_ml),label=label, ylabel="Temp (K) (Local)")
    p_T_p = plot(tt,transpose(T_p),label=label)
    p_T_h = plot(tt,transpose(T_h),label=label)
    p_T_m = plot(tt,transpose(T_m),label=label)
    p_T_s = plot(tt,transpose(T_s),label=label)
    p_F_ml = plot(tt,transpose(F_ml)*1e3, label=label, ylabel="Flow rate (L/s)")
    p_F_p = plot(tt,transpose(F_p)*1e3, label=label)
    p_F_h = plot(tt,transpose(F_h)*1e3, label=label)
    p_F_m = plot(tt,transpose(F_m)*1e3, label=label)
    p_F_s = plot(tt,transpose(F_s)*1e3, label=label)
    p_Q_ml = plot(tt,transpose(Q_ml), label=label, ylabel="Heating Rate (kW)",xlabel="Time (s) (Conf. switched)" )
    p_Q_p = plot(tt,transpose(Q_p), label=label, xlabel="Time (s) (a)")
    p_Q_h = plot(tt,transpose(Q_h), label=label, xlabel="Time (s) (b)")
    p_Q_m = plot(tt,transpose(Q_m), label=label, xlabel="Time (s) (c)")
    p_Q_s = plot(tt,transpose(Q_s), label=label, xlabel="Time (s) (d)")
    p_all = plot(p_xB_ml, p_xB_p, p_xB_h, p_xB_m, p_xB_s, p_T_ml, p_T_p, p_T_h, p_T_m, p_T_s, p_F_ml, p_F_p, p_F_h, p_F_m, p_F_s, p_Q_ml, p_Q_p, p_Q_h, p_Q_m, p_Q_s,layout=(4,5),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[2500,1500],right_margin=10mm,left_margin=18mm,linewidth=3, bottom_margin=10mm)
    # p_all=plot(p_xB_ml, p_T_ml, p_F_ml, p_Q_ml, p_xB_p, p_T_p, p_F_p, p_Q_p, p_xB_h, p_T_h, p_F_h, p_Q_h, p_xB_m, p_T_m, p_F_m, p_Q_m, p_xB_s, p_T_s, p_F_h, p_Q_s,layout=(5,4),legend=:topright,legendfont=16,xtickfontsize=16,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,framestyle=:box,size=[1500,2500],right_margin=10mm,left_margin=18mm,linewidth=3)
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\otherPerformanceComparison.pdf")
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\otherPerformanceComparison.png")
    label2 = ["Conf. Switched" "Conf.(a)" "Conf.(b)" "Conf.(c)" "Conf.(d)"]
    plot(tt[2:end],[PI_ml[2:end],PI_p[2:end], PI_h[2:end], PI_m[2:end], PI_s[2:end]], xlabel="Times (s)", ylabel = "Performance Index of xBt", label=label2,framestyle=:box,legend=:topright,legendfont=12,xtickfontsize=8,ytickfontsize=8,xguidefontsize=12,yguidefontsize=12,linewidth=3)
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\PIcomparison.pdf")
    savefig("C:\\Users\\10060\\Downloads\\Research\\posterdata\\PIcomparison.png")
end

function Plot_for_MLvsNonML()
    ml = CSV.read("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\ML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    p = CSV.read("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\parallel\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    h = CSV.read("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\hybrid\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    m = CSV.read("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\mixing\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    s = CSV.read("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\series\\noML_initial_T1_388.7_T2_388.7_T3_388.7_xB1_0.11_xB2_0.11_xB3_0.11_T0_340.0SetChange_xB_0.0.csv",DataFrame)
    tt = ml[!,"times"]
    PI_ml = ml[!,"xBt_PI"]
    PI_p = p[!,"xBt_PI"]
    PI_h = h[!,"xBt_PI"]
    PI_m = m[!,"xBt_PI"]
    PI_s = s[!,"xBt_PI"]

    xBt_ml = ml[!,"xBtinitial"]
    xBt_p = p[!,"xBtinitial"]
    xBt_h = h[!,"xBtinitial"]
    xBt_m = m[!,"xBtinitial"]
    xBt_s = s[!,"xBtinitial"]

    plot([tt[1],tt[9]],[0.11,0.11],linestyle=:dot,width=1,linecolor="Crimson",label="Set Point")
    plot!(tt[1:9],[xBt_p[1:9],xBt_ml[1:9],xBt_h[1:9],xBt_m[1:9],xBt_s[1:9]],title="Mole Frac. of Product B",label=["parallel" "ML configuration switch" "hybrid" "mixing" "sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3)

    savefig("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\xBtComparison.pdf")
    savefig("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\xBtComparison.png")
    plot(tt[2:9],[PI_p[2:9],PI_ml[2:9],PI_h[2:9],PI_m[2:9],PI_s[2:9]],title="PI of xBt",label=["parallel" "ML configuration switch" "hybrid" "mixing" "sereis"],xlabel="Times (s)",framestyle=:box,linewidth=3,legend=:bottomright)
    savefig("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\PI_xBtComparison.pdf")
    savefig("C:\\Users\\10060\\Downloads\\SummerResearch\\posterdata\\withML\\PI_xBtComparison.png")
end