using CSV, DataFrames,XLSX

df1 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, parallel=0.0)
df2 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, hybrid=0.0)
df3 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, mixing=0.0)
df4 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, series=0.0)

# Permutation on all files in the directory of Disturbance.
LoopInputList = collect(-1:-1:-20)
for i in LoopInputList
    xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\parallel\\Dist=$i.xlsx")
    xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\2and1parallel\\Dist=$i.xlsx")
    xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\mixing\\Dist=$i.xlsx")
    xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\series\\Dist=$i.xlsx")
    sh1 = xf1["Sheet1"]
    sh2 = xf2["Sheet1"]
    sh3 = xf3["Sheet1"]
    sh4 = xf4["Sheet1"]
    NewTin = sh1["C16"]
    NewxBset = sh1["B16"]
    NewT1 = sh1["F16"]
    NewT2 = sh1["G16"]
    NewT3 = sh1["H16"]
    NewxB1 = sh1["I16"]
    NewxB2 = sh1["J16"]
    NewxB3 = sh1["K16"]
    NewxBt = sh1["L16"]
    ObjValueParallel = sh1["T18"]
    ObjValueHybrid = sh2["T18"]
    ObjValueHybridMixing = sh3["T18"]
    ObjValueHybridSeries = sh4["T18"]
    push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
    push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
    push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
    push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
    println("i=: ",i)
end

LoopInputList_PositiveDisturbance = collect(1:1:40)
for i in LoopInputList_PositiveDisturbance
    xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\parallel\\positive\\Dist=$i.xlsx")
    xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\2and1parallel\\positive\\Dist=$i.xlsx")
    xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\mixing\\positive\\Dist=$i.xlsx")
    xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\series\\positive\\Dist=$i.xlsx")
    sh1 = xf1["Sheet1"]
    sh2 = xf2["Sheet1"]
    sh3 = xf3["Sheet1"]
    sh4 = xf4["Sheet1"]
    NewTin = sh1["C16"]
    NewxBset = sh1["B16"]
    NewT1 = sh1["F16"]
    NewT2 = sh1["G16"]
    NewT3 = sh1["H16"]
    NewxB1 = sh1["I16"]
    NewxB2 = sh1["J16"]
    NewxB3 = sh1["K16"]
    NewxBt = sh1["L16"]
    ObjValueParallel = sh1["T18"]
    ObjValueHybrid = sh2["T18"]
    ObjValueHybridMixing = sh3["T18"]
    ObjValueHybridSeries = sh4["T18"]
    push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
    push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
    push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
    push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
    println("i=: ",i)
end


# For Setpoint Tracking
dir_setpoint_tracking1 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\parallel"
dir_setpoint_tracking2 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\2and1 parallel"
dir_setpoint_tracking3 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\mixing"
dir_setpoint_tracking4 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\series"

x1 = readdir(dir_setpoint_tracking1)
x2 = readdir(dir_setpoint_tracking2)
x3 = readdir(dir_setpoint_tracking3)
x4 = readdir(dir_setpoint_tracking4)

for i in 1:length(x1)
    if occursin("xlsx",x1[i])
        println(x1[i])
        xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\parallel\\"*x1[i])
        xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\2and1 parallel\\"*x1[i])
        xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\mixing\\"*x1[i])
        xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\series\\"*x1[i])

        sh1 = xf1["Sheet1"]
        sh2 = xf2["Sheet1"]
        sh3 = xf3["Sheet1"]
        sh4 = xf4["Sheet1"]

        NewTin = sh1["C3"]
        NewxBset = sh1["B3"]
        NewT1 = sh1["F3"]
        NewT2 = sh1["G3"]
        NewT3 = sh1["H3"]
        NewxB1 = sh1["I3"]
        NewxB2 = sh1["J3"]
        NewxB3 = sh1["K3"]
        NewxBt = sh1["L3"]

        ObjValueParallel = sh1["T4"]
        ObjValueHybrid = sh2["T4"]
        ObjValueHybridMixing = sh3["T4"]
        ObjValueHybridSeries = sh4["T4"]
        push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
        push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
        push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
        push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
    else
        continue
    end 
    println("i=: ",i)
end

# For initial conditions Permutation
dir_initalconditions1 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\parallel"
dir_initalconditions2 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\2and1parallel"
dir_initalconditions3 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\mixing"
dir_initalconditions4 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\series"

x1 = readdir(dir_initalconditions1)
x2 = readdir(dir_initalconditions2)
x3 = readdir(dir_initalconditions3)
x4 = readdir(dir_initalconditions4)

for i in 1:length(x1)
    if occursin("xlsx",x1[i])
        println(x1[i])
        xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\parallel\\"*x1[i])
        xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\2and1parallel\\"*x1[i])
        xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\mixing\\"*x1[i])
        xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Initial conditions permutation\\series\\"*x1[i])

        sh1 = xf1["Sheet1"]
        sh2 = xf2["Sheet1"]
        sh3 = xf3["Sheet1"]
        sh4 = xf4["Sheet1"]

        NewTin = sh1["C2"]
        NewxBset = sh1["B2"]
        NewT1 = sh1["F2"]
        NewT2 = sh1["G2"]
        NewT3 = sh1["H2"]
        NewxB1 = sh1["I2"]
        NewxB2 = sh1["J2"]
        NewxB3 = sh1["K2"]
        NewxBt = sh1["L2"]

        ObjValueParallel = sh1["T4"]
        ObjValueHybrid = sh2["T4"]
        ObjValueHybridMixing = sh3["T4"]
        ObjValueHybridSeries = sh4["T4"]
        push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
        push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
        push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
        push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
    else
        continue
    end 
    println("i=: ",i)
end

CSV.write("G:\\My Drive\\Research\\SVM\\Training dataset\\Training set of parallel.csv",df1)
CSV.write("G:\\My Drive\\Research\\SVM\\Training dataset\\Training set of hybrid.csv",df2)
CSV.write("G:\\My Drive\\Research\\SVM\\Training dataset\\Training set of mixing.csv",df3)
CSV.write("G:\\My Drive\\Research\\SVM\\Training dataset\\Training set of series.csv",df4)
# file = CSV.File(IOBuffer(data); header=["times", "xBset", "T01", "T02", "T03", "Tvt1", "Tvt2", "Tvt3", "xBvt1",	"xBvt2", "xBvt3", "xBtvt", "flowvt1", "flowvt2", "flowvt3", "heatvt1", "heatvt2", "heatvt3", "PerformanceIndex", "xBtPI"],select=["xBtPI"])
# df = CSV.write("tempcsvfile.csv")
# i=0;
# # TODO When write results from reactor_reconfigure_simulator.jl delete all spaces in the column name!
# temp=0;
# for row in file
#     i+=1
#     println(i)
#     if i==18
#         println("xBset=",row.xBtPI)
#         temp=row.xBtPI
#         numbertemp=parse(Float64,temp)
#         CSV.write("G:\\My Drive\\Research\\Reconfiguration\\tempcsvfile.csv",numbertemp)
#     end
# end
# i=0;