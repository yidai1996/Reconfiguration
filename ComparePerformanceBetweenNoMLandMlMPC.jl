using DataFrames, CSV, XLSX

df = DataFrame(Tin=300, xBset=0.11, T1initial=388.7, T2initial=388.7, T3initial=388.7, xB1initial=0.11, xB2initial=0.11, xB3initial=0.11, xBtinitial=0.11,  PI_ML=0.0)
initial_condition_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration_ML_MPC\\didML\\"
x1 = readdir(initial_condition_dir)
# df[!,:PI_ML] = map( (x) -> string(round(x; digits=8)) , df[!,:PI_ML])
for i in eachindex(x1)
    # xf1 = XLSX.readxlsx(initial_condition_dir*x1[i])
    # # println(xf1)
    # sh1 = xf1["Sheet1"]
    # NewTin = sh1["C2"]
    # NewxBset = sh1["B3"]
    # NewT1 = sh1["F2"]
    # NewT2 = sh1["G2"]
    # NewT3 = sh1["H2"]
    # NewxB1 = sh1["I2"]
    # NewxB2 = sh1["J2"]
    # NewxB3 = sh1["K2"]
    # NewxBt = sh1["L2"]
    # ObjValueML = sh1["T4"]
    xf1 = CSV.read(initial_condition_dir*x1[i],DataFrame)
    NewTin = xf1[1,"T01"]
    NewxBset = xf1[2,"xBset"]
    NewT1 = xf1[1,"T1initial"]
    NewT2 = xf1[1,"T2initial"]
    NewT3 = xf1[1,"T3initial"]
    NewxB1 = xf1[1,"xB1initial"]
    NewxB2 = xf1[1,"xB2initial"]
    NewxB3 = xf1[1,"xB3initial"]
    NewxBt = xf1[1,"xBtinitial"]
    ObjValueML = xf1[3,"xBt_PI"]
    # println(ObjValueML)
    if ObjValueML>=100000
        println(ObjValueML)
        println(i)
    end
    push!(df,(NewTin,NewxBset,NewT1,NewT2,NewT3,NewxB1,NewxB2,NewxB3,NewxBt,ObjValueML))
end
check = df[!,:PI_ML]
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration_ML_MPC\\PerformanceIndexOfMLMPC.csv",df)
dff = CSV.read("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\Space_filling_sampling\\dataset\\Training set of best configurations with sorted.csv", DataFrame)
PI_NO_ML_BEST = zeros(nrow(dff))
for i in 1:nrow(dff)
    println(i)
    configuration = dff[i,"BestConfiguration"]
    if configuration == "parallel"
        PI_NO_ML_BEST[i] = dff[i,"parallel"]
    elseif configuration == "hybrid"
        PI_NO_ML_BEST[i] = dff[i,"hybrid"]
    elseif configuration == "series"
        PI_NO_ML_BEST[i] = dff[i,"series"]
    elseif configuration == "mixing"
        PI_NO_ML_BEST[i] = dff[i,"mixing"]
    else println("ERROR!")
        break
    end
    if configuration==NaN
        println("NAN ERROR!")
        break
    end
end
dfff = CSV.read("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration_ML_MPC\\PerformanceIndexOfMLMPC.csv", DataFrame)

PI_ML = dfff[!,:PI_ML]

discrepancy = PI_ML-PI_NO_ML_BEST
percentage = discrepancy./PI_NO_ML_BEST
index_max = argmax(percentage)
println(discrepancy[index_max])
dfff[!,:PI_NO_ML_BEST] = PI_NO_ML_BEST
dfff[!,:percentage] = percentage
CSV.write("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration_ML_MPC\\ComparedPerformanceIndexOfMLMPC.csv",dfff)

