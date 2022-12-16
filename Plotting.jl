using Plots, Printf, XLSX
using Plots.PlotMeasures
function MyPlots()
    xf = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Training dataset.xlsx")
    sh1 = xf["Setpoint Tracking"]
    sh2 = xf["Disturbance Rejection"]
    s = convert(Array{Float64,2},sh1["G3:G32"])
    percentage1 = 100*convert(Array{Float64,2},sh1["L3:L32"])
    d = convert(Array{Float64,2},sh2["G6:G42"])
    percentage2 = 100*convert(Array{Float64,2},sh2["L6:L42"])
    
    # p1=plot(s,percentage1,markershape = :circle,linestyle = :solid,xlabel="Set point change on xB",ylabel="C ISE percentage change (%)",framestyle=:box,label=false,xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,size=[1100,400],right_margin=0mm,left_margin=10mm,bottom_margin=9mm)
    # plot!([0,0.3],[0,0],linestyle=:dot,width=2,label=false)
    # display(p1)
    p2=plot(d,percentage2,markershape = :circle,linestyle = :solid,xlabel="Disturbance change on input temperature",ylabel="C ISE percentage change (%)",framestyle=:box,label=false,xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,size=[1100,400],right_margin=0mm,left_margin=10mm,bottom_margin=9mm)
    plot!([3,40],[0,0],linestyle=:dot,width=2,label=false)
    display(p2)
    # p2=plot(tt,X, label=false, ylabel="Final Mole Frac. of B")
    # plot!([tt[begin],1350],[0.11,0.11],linestyle=:dot,width=2,linecolor="Crimson",label="Setpoint")
    # plot!([1350,1350],[0.11,0.11+0.01],linestyle=:dot,width=2,linecolor="Crimson",label=false)
    # plot!([1350,tt[end]],[0.11+0.01,0.11+0.01],linestyle=:dot,width=2,linecolor="Crimson",label=false)
    # p3=plot(tt,transpose(Q),xlabel="Time (s)", label=label,ylabel="Heating Rate (kW)")
    # p4=plot(tt,1000*transpose(F),xlabel="Time (s)", label=label,ylabel="Flow rate (L/s)")
    # # p5=plot(tt,transpose(Xr),xlabel="Time (s)", label=label,ylabel="Mole fraction of B in each reactor")

    # # p_all=plot(p2,layout=(1,1),legend=:bottomright,xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,framestyle=:box)
    # # p_all=plot(p2,p3,layout=(1,2),legend=:bottomleft,xtickfontsize=10,ytickfontsize=10,xguidefontsize=12,yguidefontsize=12,framestyle=:box,size=[1100,400],legendfont=10,right_margin=0mm,left_margin=10mm,bottom_margin=9mm)
    # p_all=plot(p1,p2,p3,p4,layout=(2,2),legend=:topright,xtickfontsize=6,ytickfontsize=6,xguidefontsize=8,yguidefontsize=8,framestyle=:box)

    # # p1=plot(times,transpose(heatvt),label=["R1" "R2" "R3"],ylabel="Heating Rate (kW)")
    # # p2=plot(times,1000*transpose(flowvt),xlabel="Time (s)", label=["R1" "R2" "R3"],ylabel="Flow rate (L/s)")
    # # p3=plot(times,transpose(Tvt),label=["R1" "R2" "R3"],ylabel="Temp (K) (Local)")
    # # p4=plot(times,Ttvt,label=false,ylabel="Temp (K) (Global)")
    # # p5=plot(times,transpose(xBvt),label=["R1" "R2" "R3"],xlabel="Time (s)", ylabel="Mol Frac. xB (Local)")
    # # p6=plot(times,xBtvt,xlabel="Time (s)", label=false,ylabel="Mol Frac. xB (Global)")
    # # p_all=plot(p1,p3,p4,p2,p5,p6,layout=(2,3),legend=:topright,size=[800,500],margin=5mm,framestyle=:box)

    # # savefig("proposal comparison parallel to parallel +7K and +8K.pdf")
end
