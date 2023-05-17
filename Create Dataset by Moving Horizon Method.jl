using XLSX, CSV, DataFrames
# For additional data (generate from moving horizon method)
include("reactor_reconfigure_simulator.jl")
# read the initial condition from the moving horizon file
initial_condition_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\initial_condition\\parallel\\"
# initial_condition_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\initial_condition\\hybrid\\"
# initial_condition_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\initial_condition\\mixing\\"
# initial_condition_dir = "C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\initial_condition\\series\\"
x1 = readdir(initial_condition_dir)

for i in eachindex(x1)
    xf1 = XLSX.readxlsx(initial_condition_dir*x1[i])
    println(xf1)
    sh1 = xf1["Sheet1"]
    
    for j in 3:7
        NewTin = sh1["C"*string(j)] - 300
        NewxBset = sh1["B"*string(j)] - 0.11
        NewT1 = sh1["F"*string(j)]
        NewT2 = sh1["G"*string(j)]
        NewT3 = sh1["H"*string(j)]
        NewxB1 = sh1["I"*string(j)]
        NewxB2 = sh1["J"*string(j)]
        NewxB3 = sh1["K"*string(j)]
        NewxBt = sh1["L"*string(j)]
        new_initial_condition = zeros(3,3)
        new_initial_condition[:,1] .= 300
        new_initial_condition[:,2] = [NewT1 NewT2 NewT3]
        new_initial_condition[:,3] = [NewxB1 NewxB2 NewxB3]
        println(new_initial_condition," ", NewTin)
        MPC_tracking("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration\\parallel", [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+NewTin 0+NewTin;0+NewTin 0+NewTin;0+NewTin 0+NewTin],[0;0;0+NewxBset],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],0,new_initial_condition;tmax=1500,print=false,save_plots=false,plot_name="all_plots.pdf")
        MPC_tracking("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration\\hybrid", [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+NewTin 0+NewTin;0+NewTin 0+NewTin;0+NewTin 0+NewTin],[0;0;0+NewxBset],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],0,new_initial_condition;tmax=1500,print=false,save_plots=false,plot_name="all_plots.pdf")
        MPC_tracking("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration\\mixing", [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0],[0+NewTin 0+NewTin;0+NewTin 0+NewTin;0+NewTin 0+NewTin],[0;0;0+NewxBset],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],0,new_initial_condition;tmax=1500,print=false,save_plots=false,plot_name="all_plots.pdf")
        MPC_tracking("C:\\Users\\yid\\TemporaryResearchDataStorage\\Reconfiguration\\additional_data\\PreDataSetForReconfiguration\\series", [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0],[0+NewTin 0+NewTin;0+NewTin 0+NewTin;0+NewTin 0+NewTin],[0;0;0+NewxBset],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],0,new_initial_condition;tmax=1500,print=false,save_plots=false,plot_name="all_plots.pdf")
        # break
    end
end

# input the initial condition into MPC system for 4 configurations and generate new data
