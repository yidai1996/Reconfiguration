# RawDataset
using CSV, DataFrames,XLSX

df1 = DataFrame(xB1=0.11, xB2=0.11, xB3=0.11, F1=0.07, F2=0.07, F3=0.07, Q1=200, Q2=200, Q3=200, T1=388.7, T2=388.7, T3=388.7, Tin1=300, Tin2=300, Tin3=300, xtot=0.11, xBset=0.11, T1set=388.7, T2set=388.7, T3set=388.7)
df2 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, hybrid=0.0)
df3 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, mixing=0.0)
df4 = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11, series=0.0)

# # Permutation on all files in the directory of Disturbance.
# LoopInputList = collect(-1:-1:-20)
# for i in LoopInputList
#     xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\parallel\\Dist=$i.xlsx")
#     xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\2and1parallel\\Dist=$i.xlsx")
#     xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\mixing\\Dist=$i.xlsx")
#     xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\series\\Dist=$i.xlsx")
#     sh1 = xf1["Sheet1"]
#     sh2 = xf2["Sheet1"]
#     sh3 = xf3["Sheet1"]
#     sh4 = xf4["Sheet1"]
#     NewTin = sh1["C16"]
#     NewxBset = sh1["B16"]
#     NewT1 = sh1["F16"]
#     NewT2 = sh1["G16"]
#     NewT3 = sh1["H16"]
#     NewxB1 = sh1["I16"]
#     NewxB2 = sh1["J16"]
#     NewxB3 = sh1["K16"]
#     NewxBt = sh1["L16"]
#     ObjValueParallel = sh1["T18"]
#     ObjValueHybrid = sh2["T18"]
#     ObjValueHybridMixing = sh3["T18"]
#     ObjValueHybridSeries = sh4["T18"]
#     push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
#     push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
#     push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
#     push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
#     println("i=: ",i)
# end

# LoopInputList_PositiveDisturbance = collect(1:1:40)
# for i in LoopInputList_PositiveDisturbance
#     xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\parallel\\positive\\Dist=$i.xlsx")
#     xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\2and1parallel\\positive\\Dist=$i.xlsx")
#     xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\mixing\\positive\\Dist=$i.xlsx")
#     xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Disturbance rejection\\series\\positive\\Dist=$i.xlsx")
#     sh1 = xf1["Sheet1"]
#     sh2 = xf2["Sheet1"]
#     sh3 = xf3["Sheet1"]
#     sh4 = xf4["Sheet1"]
#     NewTin = sh1["C16"]
#     NewxBset = sh1["B16"]
#     NewT1 = sh1["F16"]
#     NewT2 = sh1["G16"]
#     NewT3 = sh1["H16"]
#     NewxB1 = sh1["I16"]
#     NewxB2 = sh1["J16"]
#     NewxB3 = sh1["K16"]
#     NewxBt = sh1["L16"]
#     ObjValueParallel = sh1["T18"]
#     ObjValueHybrid = sh2["T18"]
#     ObjValueHybridMixing = sh3["T18"]
#     ObjValueHybridSeries = sh4["T18"]
#     push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
#     push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
#     push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
#     push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
#     println("i=: ",i)
# end


# # For Setpoint Tracking
# dir_setpoint_tracking1 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\parallel"
# dir_setpoint_tracking2 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\2and1 parallel"
# dir_setpoint_tracking3 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\mixing"
# dir_setpoint_tracking4 = "G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\series"

# x1 = readdir(dir_setpoint_tracking1)
# x2 = readdir(dir_setpoint_tracking2)
# x3 = readdir(dir_setpoint_tracking3)
# x4 = readdir(dir_setpoint_tracking4)

# for i in 1:length(x1)
#     if occursin("xlsx",x1[i])
#         println(x1[i])
#         xf1 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\parallel\\"*x1[i])
#         xf2 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\2and1 parallel\\"*x1[i])
#         xf3 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\mixing\\"*x1[i])
#         xf4 = XLSX.readxlsx("G:\\My Drive\\Research\\SVM\\Training dataset\\Raw data\\Setpoint tracking\\series\\"*x1[i])

#         sh1 = xf1["Sheet1"]
#         sh2 = xf2["Sheet1"]
#         sh3 = xf3["Sheet1"]
#         sh4 = xf4["Sheet1"]

#         NewTin = sh1["C3"]
#         NewxBset = sh1["B3"]
#         NewT1 = sh1["F3"]
#         NewT2 = sh1["G3"]
#         NewT3 = sh1["H3"]
#         NewxB1 = sh1["I3"]
#         NewxB2 = sh1["J3"]
#         NewxB3 = sh1["K3"]
#         NewxBt = sh1["L3"]

#         ObjValueParallel = sh1["T4"]
#         ObjValueHybrid = sh2["T4"]
#         ObjValueHybridMixing = sh3["T4"]
#         ObjValueHybridSeries = sh4["T4"]
#         push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))
#         push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
#         push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
#         push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
#     else
#         continue
#     end 
#     println("i=: ",i)
# end

# For initial conditions Permutation
dir_initalconditions1 = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\parallel"
dir_initalconditions2 = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\hybrid"
dir_initalconditions3 = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\mixing"
dir_initalconditions4 = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\series"

x1 = readdir(dir_initalconditions1)
x2 = readdir(dir_initalconditions2)
x3 = readdir(dir_initalconditions3)
x4 = readdir(dir_initalconditions4)

for i in 1:length(x1)
    xf1 = XLSX.readxlsx("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\parallel\\"*x1[i])
    sh1 = xf1["Sheet1"]
    NewTin = sh1["C2"]
    NewxBset = sh1["B3"]
    NewT1 = sh1["F2"]
    NewT2 = sh1["G2"]
    NewT3 = sh1["H2"]
    NewxB1 = sh1["I2"]
    NewxB2 = sh1["J2"]
    NewxB3 = sh1["K2"]
    NewxBt = sh1["L2"]
    ObjValueParallel = sh1["T4"]
    push!(df1,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueParallel))

end

for i in eachindex(x2)
    xf2 = XLSX.readxlsx("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\hybrid\\"*x2[i])
    sh = xf2["Sheet1"]
    NewTin = sh["C2"]
    NewxBset = sh["B3"]
    NewT1 = sh["F2"]
    NewT2 = sh["G2"]
    NewT3 = sh["H2"]
    NewxB1 = sh["I2"]
    NewxB2 = sh["J2"]
    NewxB3 = sh["K2"]
    NewxBt = sh["L2"]
    ObjValueHybrid = sh["T4"]
    push!(df2,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybrid))
end

for i in eachindex(x3)
    xf3 = XLSX.readxlsx("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\mixing\\"*x3[i])
    sh = xf3["Sheet1"]
    NewTin = sh["C2"]
    NewxBset = sh["B3"]
    NewT1 = sh["F2"]
    NewT2 = sh["G2"]
    NewT3 = sh["H2"]
    NewxB1 = sh["I2"]
    NewxB2 = sh["J2"]
    NewxB3 = sh["K2"]
    NewxBt = sh["L2"]
    ObjValueHybridMixing = sh["T4"]
    push!(df3,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridMixing))
end

for i in eachindex(x4)
    xf4 = XLSX.readxlsx("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\series\\"*x4[i])
    sh = xf4["Sheet1"]
    NewTin = sh["C3"]
    NewxBset = sh["B3"]
    NewT1 = sh["F3"]
    NewT2 = sh["G3"]
    NewT3 = sh["H3"]
    NewxB1 = sh["I3"]
    NewxB2 = sh["J3"]
    NewxB3 = sh["K3"]
    NewxBt = sh["L3"]
    ObjValueHybridSeries = sh["T4"]
    push!(df4,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueHybridSeries))
end

CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\dataset\\Training set of parallel.csv",df1)
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\dataset\\Training set of hybrid.csv",df2)
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\dataset\\Training set of mixing.csv",df3)
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\dataset\\Training set of series.csv",df4)
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
#         temp=row.xBtPIreturn the index of sorted element of a given array julia
#         numbertemp=parse(Float64,temp)
#         CSV.write("G:\\My Drive\\Research\\Reconfiguration\\tempcsvfile.csv",numbertemp)
#     end
# end
# i=0;

# Find the best configuration and create a new CSV file for it 
dff = CSV.read("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\dataset\\Training set of all configurations.csv", DataFrame)
row_number = nrow(dff)
x = zeros(row_number,4)
x[:,1] = dff.parallel
x[:,2] = dff.hybrid
x[:,3] = dff.mixing
x[:,4] = dff.series
index = zeros(row_number, 4)
for i in 1:row_number
    index[i,:] = sortperm(x[i,:])
end

# println(index[4502,:])
# min_index = argmin(x; dims = 2)
# best_configuration_column = String[]
# second_best_configuration_column = String[]
# third_best_configuration_column = String[]
# worst_configuration_column = String[]

Astring = Matrix{String}(undef, row_number, 4)

for i in 1:row_number
    c_row = index[i,:]
    for j in 1:4
        c = c_row[j]
        if c == 1
            Astring[i,j] = "parallel"
        elseif c==2
            Astring[i,j] = "hybrid"
        elseif c==3
            Astring[i,j] = "mixing"
        else
            Astring[i,j] = "series"
        end
    end
end
print(Astring[:,1])
# df_bestconfiguration = select!(dff, Not([:parallel, :hybrid, :mixing, :series]))
# df_bestconfiguration[!,:BestConfiguration] = best_configuration_column
# CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\dataset\\Training set of best configurations.csv",df_bestconfiguration)

df_bestconfiguration = dff
df_bestconfiguration[!,:BestConfiguration] = Astring[:,1]
df_bestconfiguration[!,:SecondBestConfiguration] = Astring[:,2]
df_bestconfiguration[!,:ThirdBestConfiguration] = Astring[:,3]
df_bestconfiguration[!,:WorstBestConfiguration] = Astring[:,4]
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration3\\dataset\\Training set of best configurations with sorted.csv",df_bestconfiguration)


 