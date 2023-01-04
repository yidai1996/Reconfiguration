using Plots, Printf, XLSX
using Plots.PlotMeasures
function MyPlots()
    xf = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Training dataset.xlsx") # Feel free to try any files in your computer
    sh1 = xf["Setpoint Tracking"]
    sh2 = xf["Disturbance Rejection"]
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
