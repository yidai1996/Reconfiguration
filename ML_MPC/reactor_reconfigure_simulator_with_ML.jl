# For any given number of reactors and potential configurations
# https://aiche.onlinelibrary.wiley.com/doi/full/10.1002/aic.11801
using Plots, JuMP, DifferentialEquations, NLsolve, BenchmarkTools, Ipopt
using MathOptInterface, Printf, ProgressBars, DelimitedFiles, Profile, XLSX, CSV
using DataFrames
using MLJ
# include("G:\\My Drive\\Research\\SVM\\Git\\Reclone\\Adaboost.jl")

# include("permutation.jl")

function loadProcessData(N::Int,n,initial_values;print=true)
    # global F0=9/3600/N #m^3/s
    global Vlittle=0.5/3
    # Parallel
    global V=fill(Vlittle, N) #m^3
    global d_H1=-6e4 #KJ/kmol
    global d_H2=-7e4 #KJ/kmol
    global k1=2.77e3 #s^-1
    global k2=2.5e3 #s^-1
    global E1=5e4 #KJ/kmol
    global E2=6e4 #KJ/kmol
    global c_p=4.2 #KJ/kg/K
    global mass=.00279 #kmol/kg
    global rho=1000 #kg/m^3
    global R_gas=8.314 #KJ/kmol/K
    global xA0=1
    global xB0=0
    # 4R parallel
    # global T0=[300 300 300 300]
    # global Ts=[388.7 388.7 388.7 388.7]
    # global xBs=[0.11 0.11 0.11 0.11]
    # global xAs=[1-xBs[1] 1-xBs[2] 1-xBs[3] 1-xBs[4]]
    # 3R parallel
    # global T0=[300 300 300]
    # global Ts=[388.7 388.7 388.7]
    # global xBs=[0.11 0.11 0.11]
    # global xAs=[1-xBs[1] 1-xBs[2] 1-xBs[3]]
    # global T0=[300 300]
    # global Ts=[388.7 388.7]
    # global xBs=[0.11 0.11]
    # global xAs=[1-xBs[1] 1-xBs[2]]
    # 3R 2&1parallel
    # global T0=[300 300 300]
    # global Ts=[370 388.7 388.7]
    # global xBs=[0.055 0.11 0.11]
    # global xAs=[1-xBs[1] 1-xBs[2] 1-xBs[3]]
    # 3R series
    # global T0=[300 300 300]
    # global Ts=[370 380 388.7]
    # global xBs=[0.055 0.08 0.11]
    # global xAs=[1-xBs[1] 1-xBs[2] 1-xBs[3]]
    # 3R mixing
    # global T0=[300 300 300]
    # global Ts=[370 370 388.7]
    # global xBs=[0.055 0.055 0.11]
    # global xAs=[1-xBs[1] 1-xBs[2] 1-xBs[3]]

    # 4R mixing
    # global T0=[300 300 300 300]
    # global Ts=[370 375 380 388.7]
    # global xBs=[0.055 0.07 0.085 0.11]
    # global xBs=[0.0515 0.0752 0.1038 0.115]
    # global xAs=[1-xBs[1] 1-xBs[2] 1-xBs[3] 1-xBs[4]]
    # global T0=fill(initial_values[1],N) #K
    # global Ts=fill(initial_values[2],N) # will change with different input n and other initial conditions
    # global xBs=fill(initial_values[3],N) # will change with different input n and other initial conditions
    # global xAs=fill(1-initial_values[3],N) # will change with different input n and other initial conditions
    # global Ftest=0.000709
    global Ftest=0.000709 # For 4R
    # 2R intial condition
    # global T0=[300 300] #K
    # global Ts=[388.7;388.7] # will change with different input n and other initial conditions
    # global xBs=[0.11;0.11] # will change with different input n and other initial conditions
    # global xAs=[1-xBs[1];1-xBs[2]] # will change with different input n and other initial conditions

    global T0=initial_values[:,1] #K
    global Ts=initial_values[:,2] # will change with different input n and other initial conditions
    # println("Ts=",Ts)
    global xBs=initial_values[:,3] # will change with different input n and other initial conditions
    global xAs=1 .- xBs # will change with different input n and other initial conditions

    global F0=(-k1*exp(-E1/R_gas/Ts[1])*(1-xBs[1])+(k2*exp(-E2/R_gas/Ts[1])*xBs[1]))*V/(xB0-xBs[1])
    global Flow0=zeros(N+1,N+1)
    global Q_nom=zeros(N)
    global F_nom
    Q_nom,F_nom,ini_lookup=findSS_all(T0,Ts,xBs,n,print=print)
    for k=1:length(ini_lookup)
        for i=1:N+1
            for j=1:N+1
                if ini_lookup[k][1]==i&&ini_lookup[k][2]==j
                    Flow0[i,j]=F_nom[k]
                end
            end
        end
    end

    if print
        println("Parameters Loaded!")
    end
end

function MPC_solve(xBset,Tset,n,Flow,T0_inreal,T_0real,xA_0real,xB_0real,q_T,q_xA,q_xB,r_heat,r_flow,dt,P,N;heat_init=0,flow_init=0,print=true)
    # println("n=",n)
    global count

    K=round(Int,P/dt)

    MPC=JuMP.Model(Ipopt.Optimizer)
    # MPC=Model(() -> BARON.Optimizer(MaxTime=10000))
    MOI.set(MPC, MOI.RawOptimizerAttribute("print_level"), 1)

    T0_in=T0_inreal
    T_0=T_0real
    xA_0=xA_0real
    xB_0=xB_0real
    # println("T0_in=",T0_in)
    # println("T0_0=",T_0)
    # println("xB=",xB_0)

    # Only steady states for input streams from outside instead of other reactors
    # heat_ss,flow_ss=findSS_all(T0_in,Ts,xBs,n,Flow)
    # println(Tset,xBset)
    heat_ss,flow_ss,mpclook=findSS_all(T0_in,Tset,xBset,n,print=print)


    if print
        println("Heat_ss=",heat_ss)
        println("Flow_ss=",flow_ss)
        println(mpclook)
    end

    for k=1:length(mpclook)
        for i=1:N+1
            for j=1:N+1
                if mpclook[k][1]==i&&mpclook[k][2]==j
                    if flow_ss[k]<0
                        Flow[i,j]=1e-7
                    else  Flow[i,j]=flow_ss[k]
                    end
                end
            end
        end
    end
    for i=1:N
        if heat_ss[i]<0
            heat_ss[i]=1e-7
        end
    end

    if print
        println("Flow=",Flow)
    end

    xA_guess=zeros(N,K+1)
    xB_guess=zeros(N,K+1)
    T_guess=zeros(N,K+1)
    for i=1:N
        xA_guess[i,1]=xA_0[i]
        xB_guess[i,1]=xB_0[i]
        T_guess[i,1]=T_0[i]
    end
    xB_tot_guess=zeros(K+1)
    xB_tot_guess[1]=sum(n[i,N+1]*Flow[i,N+1]*xB_guess[i,1] for i=1:N)/sum(n[i,N+1]*Flow[i,N+1] for i=1:N)

    for k=1:K
        for i=1:N
            T_guess[i,k+1] = (1/V[i]*(sum(n[j,i]*Flow[j,i]*T_guess[j,k] for j=1:N) + n[N+1,i]*Flow[N+1,i]*T0_in[i]- sum(n[i,j]*Flow[i,j]*T_guess[i,k] for j=1:N+1)) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_guess[i,k])*xA_guess[i,k])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_guess[i,k])*xB_guess[i,k]) + heat_ss[i]/rho/c_p/V[i])*dt + T_guess[i,k]
            xA_guess[i,k+1] = (1/V[i]*(sum(n[j,i]*Flow[j,i]*xA_guess[j,k] for j=1:N) + n[N+1,i]*Flow[N+1,i]*xA0 - sum(n[i,j]*Flow[i,j]*xA_guess[i,k] for j=1:N+1)) + (-k1*exp(-E1/R_gas/T_guess[i,k])*xA_guess[i,k]))*dt + xA_guess[i,k]
            xB_guess[i,k+1] = (1/V[i]*(sum(n[j,i]*Flow[j,i]*xB_guess[j,k] for j=1:N) - sum(n[i,j]*Flow[i,j]*xB_guess[i,k] for j=1:N+1)) + k1*exp(-E1/R_gas/T_guess[i,k])*xA_guess[i,k] + (-k2*exp(-E2/R_gas/T_guess[i,k])*xB_guess[i,k]))*dt + xB_guess[i,k]
        end
        xB_tot_guess[k+1] = sum(n[i,N+1]*Flow[i,N+1]*xB_guess[i,k] for i=1:N)/sum(n[i,N+1]*Flow[i,N+1] for i=1:N)
        # Tt_guess[k+1]=sum((n[i,N+1])*Flow[i,N+1]*T_guess[i,k] for i=1:N)/sum(n[i,N+1]*Flow[i,N+1] for i=1:N)
    end

    if print
        println("xB_guess=",xB_guess)
        println("xBt_guess=",xB_tot_guess)
    end

    JuMP.@variables MPC begin
        # Q[i=1:N,k=0:K-1], (lower_bound=0.2*heat_ss[i], upper_bound=1.8*heat_ss[i],start=heat_ss[i])# Q of the reactors
        Q[i=1:N,k=0:K-1], (lower_bound=0, upper_bound=1.8*heat_ss[i],start=heat_ss[i])# Q of the reactors
        # F[i=1:N+1,j=1:N+1,k=0:K-1], (lower_bound=n[i,j]*0.2*Flow[i,j], upper_bound=n[i,j]*1.8*Flow[i,j],start=Flow[i,j])# Flowrate between reactors
        F[i=1:N+1,j=1:N+1,k=0:K-1], (lower_bound=n[i,j]*0.2*Flow[i,j], upper_bound=n[i,j]*(1+Flow[i,j]),start=Flow[i,j])# Flowrate between reactors
        # F[i=1:N+1,j=1:N+1,k=0:K-1], (lower_bound=n[i,j]*0.2*Flow[i,j],start=Flow[i,j])# Flowrate between reactors
        T[i=1:N,k=0:K], (lower_bound=T0[i],upper_bound=2000,start=T_guess[i,k+1])
        xA[i=1:N,k=0:K], (lower_bound=0, upper_bound=1,start=xA_guess[i,k+1])
        xB[i=1:N,k=0:K], (lower_bound=0, upper_bound=1,start=xB_guess[i,k+1])
        xBt[k=0:K], (lower_bound=0, upper_bound=1,start=xB_tot_guess[k+1])

        # m[k=0:K], (lower_bound=0, upper_bound=1,start=m_init)
    end

    JuMP.@constraints MPC begin
        T_init[i=1:N], T[i,0]==T_0[i]
        xA_init[i=1:N], xA[i,0]==xA_0[i]
        xB_init[i=1:N], xB[i,0]==xB_0[i]
        xBt_init, xBt[0]==sum(n[i,N+1]*Flow[i,N+1]*xB_0[i] for i=1:N)/sum(n[i,N+1]*Flow[i,N+1] for i=1:N)
        MassB[i=1:N,k=0:K-1], sum(n[j,i]*F[j,i,k] for j=1:N+1)==sum(n[i,j]*F[i,j,k] for j=1:N+1)
    end

    JuMP.@NLconstraints MPC begin
        Temp[i=1:N,k=0:K-1], T[i,k+1] == (1/V[i]*(sum(n[j,i]*F[j,i,k]*T[j,k] for j=1:N) + n[N+1,i]*F[N+1,i,k]*T0_in[i]- sum(n[i,j]*F[i,j,k]*T[i,k] for j=1:N+1)) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T[i,k])*xA[i,k])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T[i,k])*xB[i,k]) + Q[i,k]/rho/c_p/V[i])*dt + T[i,k]
        MoleFractionxA[i=1:N,k=0:K-1], xA[i,k+1] == (1/V[i]*(sum(n[j,i]*F[j,i,k]*xA[j,k] for j=1:N) + n[N+1,i]*F[N+1,i,k]*xA0 - sum(n[i,j]*F[i,j,k]*xA[i,k] for j=1:N+1)) + (-k1*exp(-E1/R_gas/T[i,k])*xA[i,k]))*dt + xA[i,k]
        MoleFractionxB[i=1:N,k=0:K-1], xB[i,k+1] == (1/V[i]*(sum(n[j,i]*F[j,i,k]*xB[j,k] for j=1:N) - sum(n[i,j]*F[i,j,k]*xB[i,k] for j=1:N+1)) + k1*exp(-E1/R_gas/T[i,k])*xA[i,k] + (-k2*exp(-E2/R_gas/T[i,k])*xB[i,k]))*dt + xB[i,k]
        OutputMoleFraction[k=0:K-1], xBt[k+1] == sum(n[i,N+1]*F[i,N+1,k]*xB[i,k+1] for i=1:N)/sum(n[i,N+1]*F[i,N+1,k] for i=1:N)
    end

    JuMP.@objective(MPC,Min,sum(q_T*(T[i,k]-Tset[i])^2 for i=1:N for k=0:K)+sum(q_xB*(xBt[k]-xBset[end])^2 for k=0:K)+sum(r_heat*(Q[i,k]-Q[i,k-1])^2 for i=1:N for k=1:K-1) + sum(r_flow*(n[i,j]*F[i,j,k]-n[i,j]*F[i,j,k-1])^2 for i=1:N+1 for j=1:N+1 for k=1:K-1) + sum(r_heat*(Q[i,0]-Q[i,K-1])^2 for i=1:N) + sum(r_flow*(n[i,j]*F[i,j,0]-n[i,j]*F[i,j,K-1])^2 for i=1:N+1 for j=1:N+1))

    JuMP.optimize!(MPC)


    st=MathOptInterface.RawStatusString()
    if st=="INFEASIBLE_POINT"
        println(xA1_guess,xA2_guess)
        error("Solver infeasible, problem stopping")
    end
    # obj=getobjectivevalue(MPC) # works for Julia 1.15.3
    obj=JuMP.objective_value(MPC) # works for Julia 1.17.2
    if print
        println("Obj in MPC=",obj)
    end

    results_T=JuMP.value.(T)
    # println("results_T=",results_T)
    results_xB=JuMP.value.(xB)
    # println("results_xB=",results_xB)
    results_xBt=JuMP.value.(xBt)
    # println("results_xBt=",results_xBt)
    results_heat = JuMP.value.(Q)
    # println("results_heat=",results_heat)
    results_flow = JuMP.value.(F)
    # println("results_flow=",results_flow)
    results_heat0 = JuMP.value.(Q[:,0])
    results_flow0 = JuMP.value.(F[:,:,0])

    obj_T=sum(q_T*(results_T[i,k]-Tset[i])^2 for i=1:N for k=0:K)
    obj_xBt=sum(q_xB*(results_xBt[k]-xBset[end])^2 for k=0:K)
    obj_Q=sum(r_heat*(results_heat[i,k]-results_heat[i,k-1])^2 for i=1:N for k=1:K-1)+sum(r_heat*(results_heat[i,0]-results_heat[i,K-1])^2 for i=1:N)
    obj_F=sum(r_flow*(n[i,j]*results_flow[i,j,k]-n[i,j]*results_flow[i,j,k-1])^2 for i=1:N for j=1:N+1 for k=1:K-1)+sum(r_flow*(n[i,j]*results_flow[i,j,0]-n[i,j]*results_flow[i,j,K-1])^2 for i=1:N for j=1:N+1)

    # println("results_heat=",results_heat)
    # println("results_flow=",results_flow)
    # println("Obj_T= ",obj_T)
    # println("Obj_xBt= ",obj_xBt)
    # println("Obj_Q= ",obj_Q)
    # println("Obj_F= ",obj_F)

    if print
        println("soln_heat=",results_heat0)
        println("soln_flow=",results_flow0)
    end
    return results_heat0, results_flow0, obj_xBt,obj_T,obj_Q,obj_F,obj

end

# SetChange_xB = [1xN]

function MPC_tracking(out_dir, n1::Array{Int,2},n2,Dist_T0,SetChange_xB,SetChange_T,q_T,q_xA,q_xB,r_heat,r_flow,dt,P,
    dist_time,setpoint_time,initial_values; tmax=200,print=true,save_plots=false,plot_name="all_plots.png",MLcheck=false) # This is for continous disturbance on the (unstable) input temperature
    # (runs the moving horizon loop for set point tracking)
    # N=length(Dist_T0)
    # When testing continous disturbance system, the Dist_T0 contains the beginning point
    println(tmax)
    global N=size(n1)[1]-1
    if print
        println("N=",N)
    end
    l=length(dist_time)
    ll=length(setpoint_time)
    # Check the length of disturbance vectors and dist_time vector are the same
    if print
        if size(Dist_T0)[2] == l
            println("The length of disturbance variables == the one of dist_time vector=",l)
        else println("The length of disturbance variables are not equal to the one of dist_time vector")
            return
        end
    end
    loadProcessData(N,n1,initial_values,print=print)
    # loadProcessData(N,n1,initial_setpoints,print=print)

    time_steps=round(Int,tmax/dt)

    global times=zeros(time_steps+1)
    pt=round(Int,P/dt)

    global count=1

    global ObjValue=zeros(time_steps) # To storage the optimal objective value from MPC for each iteration
    global ObjValue=zeros(time_steps) # To storage the optimal objective value from MPC for each iteration
    obj_output_total=zeros(time_steps+1)
    obj_output_xBt=zeros(time_steps+1)
    obj_output_T=zeros(time_steps+1)
    obj_output_Q=zeros(time_steps+1)
    obj_output_F=zeros(time_steps+1)

    global Tvt=zeros(N,time_steps+1)
    global xAvt=zeros(N,time_steps+1)
    global xBvt=zeros(N,time_steps+1)
    global heatvt=zeros(N,time_steps+1)
    global flowvt=zeros(N+1,N+1,time_steps+1)
    global adjacentM=zeros(N+1,N+1,time_steps+1)
    global record_configuration=zeros(time_steps+1)
    global xBtvt=zeros(1,time_steps+1)
    newstate=zeros(3*N)
    # Y=zeros(time_steps+1)
    global T0_invt=zeros(N,time_steps+1)
    global xBsetpoint=zeros(N,time_steps+1)
    global Tsetpoint=zeros(N,time_steps+1)

    # global recordFindSS=zeros()
    # global recordStepAll=zeros()
    nn1 = [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0];
    nn2 = [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0];
    nn3 = [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0];
    nn4 = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0];
    # xBsetpoint[3,1]=xBs[3]+SetChange_xB[3]
    xBsetpoint[:,1]=xBs
    Tsetpoint[:,1]=Ts
    T0_invt[:,1]=initial_values[:,1]
    Tvt[:,1]=initial_values[:,2]
    xBvt[:,1]=initial_values[:,3]
    xAvt[:,1]=1 .- xBs
    heatvt[:,1]=Q_nom
    flowvt[1:N+1,1:N+1,1]=Flow0
    xBtvt[1]=sum(n1[i,N+1]*flowvt[i,N+1,1]*xBvt[i,1] for i=1:N)/sum(n1[i,N+1]*flowvt[i,N+1,1] for i=1:N)
    # mach = machine("KNN_Zavreal_not_bigdata.jl")
    # mach = machine("KNN_t.jl")
    mach_3 = []
    My_Machines = ["KNN_t.jl", "Tree_max_depth_9.jl", "SVM.jl"]
    if MLcheck==true
        # Inplement reconfiguration based on ML mathematical guideline.
        features = DataFrame(Tin=T0_invt[1,1], xBset=xBsetpoint[3,1], T1initial=Tvt[1,1], T2initial=Tvt[2,1], T3initial=Tvt[3,1], xB1initial=xBvt[1,1], xB2initial=xBvt[2,1], xB3initial=xBvt[3,1], xBtinitial=xBtvt[1])
        println(features)
        # configuration = predict_mode(mach, features) # For member ML 
        
        for i in eachindex(My_Machines)
            push!(mach_3, machine(My_Machines[i]))
        end
        configuration = finalH([4.683, 2.592, 2.619], features, mach_3) # For Adaboost

        println(configuration)
       if configuration == ["parallel"]
           adjacentM[:,:,1]=nn1
            record_configuration[1]=1
       elseif configuration == ["hybrid"]
           adjacentM[:,:,1]=nn2
           record_configuration[1]=2
       elseif configuration == ["mixing"]
           adjacentM[:,:,1]=nn3
           record_configuration[1]=3
       elseif configuration == ["series"]
           adjacentM[:,:,1]=nn4
           record_configuration[1]=4
       else pringln("ERROR IN RECONFIGURATION MACHINE")
       end
    else
        adjacentM[1:N+1,1:N+1,1]=n1
        record_configuration[1]=1
   end
    
    
    times[1]=0
    tt=1
    

    for tt=1:time_steps
        if MLcheck==true
             # Inplement reconfiguration based on ML mathematical guideline.
             features = DataFrame(Tin=T0_invt[1,tt], xBset=xBsetpoint[3, tt], T1initial=Tvt[1,tt], T2initial=Tvt[2,tt], T3initial=Tvt[3,tt], xB1initial=xBvt[1,tt], xB2initial=xBvt[2,tt], xB3initial=xBvt[3,tt], xBtinitial=xBtvt[tt])
         
             println(features)
            #  configuration = predict_mode(mach, features)

            for i in eachindex(My_Machines)
                push!(mach_3, machine(My_Machines[i]))
            end
            configuration = finalH([4.683, 2.592, 2.619], features, mach_3) # For Adaboost

            println(configuration)
            if configuration == ["parallel"]
                adjacentM[:,:,tt+1]=nn1
                 record_configuration[tt+1]=1
            elseif configuration == ["hybrid"]
                adjacentM[:,:,tt+1]=nn2
                record_configuration[tt+1]=2
            elseif configuration == ["mixing"]
                adjacentM[:,:,tt+1]=nn3
                 record_configuration[tt+1]=3
            elseif configuration == ["series"]
                adjacentM[:,:,tt+1]=nn4
                record_configuration[tt+1]=4
            else pringln("ERROR IN RECONFIGURATION MACHINE")
                break
            end
        end
 
        resultsheatvt,resultsflowvt,obj_output_xBt[tt+1],obj_output_T[tt+1],obj_output_Q[tt+1],obj_output_F[tt+1],obj_output_total[tt+1]=MPC_solve(xBsetpoint[:,tt],Tsetpoint[:,tt],adjacentM[:,:,tt],flowvt[:,:,tt],T0_invt[:,tt],Tvt[:,tt],xAvt[:,tt],xBvt[:,tt],q_T,q_xA,q_xB,r_heat,r_flow,dt,P,N;
            heat_init=heatvt[1,tt],flow_init=flowvt[1,1,tt],print=print)

        for i=1:N
            heatvt[i,tt+1]=resultsheatvt[i]
            flowvt[:,:,tt+1]=resultsflowvt
        end
        if print
            println("count(tt)=",count)
        end

        for i=1:N
            # println("The disturb has been added for calculating system reactions")
            j=Int(l)
            # println("j=",j)
            while j>=0
                if j==0
                    T0_invt[i,tt]=T0_invt[i,tt]
                    break
                end
                if tt>=dist_time[j]
                    T0_invt[i,tt]=T0[i]+Dist_T0[i,j]
                    # println("i=",i," j=",j," T0_invt[i,tt]=",T0_invt[i,tt]," dist_time[j]=",dist_time[j])
                    break
                else j=j-1
                end

            end
        end

        newstate=MPC_step_all(T0_invt[:,tt],Tvt[:,tt],xAvt[:,tt],xBvt[:,tt],heatvt[:,tt+1],flowvt[:,:,tt+1],adjacentM[:,:,tt],dt,print=print)
        for aa=0:N-1
            Tvt[aa+1,tt+1]=newstate[3*aa+1]
            xAvt[aa+1,tt+1]=newstate[3*aa+2]
            xBvt[aa+1,tt+1]=newstate[3*aa+3]
        end

        for i=1:N
            j=Int(l)
            # println("j=",j)
            while j>=0
                if j==0
                    T0_invt[i,tt+1]=T0_invt[i,tt]
                    break
                end
                if tt>=dist_time[j]
                    T0_invt[i,tt+1]=T0[i]+Dist_T0[i,j]
                    # println("i=",i," j=",j," T0_invt[i,tt+1]=",T0_invt[i,tt+1]," dist_time[j]=",dist_time[j])
                    break
                else j=j-1
                end

            end
        end

        

        for i=1:N
            j=Int(ll) # j is number of times we change setpoint
            while j>=0 #
                if j==0
                    xBsetpoint[i,tt+1]=xBs[i]
                    Tsetpoint[i,tt+1]=Ts[i]
                    if MLcheck==false
                        adjacentM[:,:,tt+1]=n1
                    end
                    break
                end
                if tt>=setpoint_time[j]
                    xBsetpoint[i,tt+1]=xBs[i]+SetChange_xB[i]
                    Tsetpoint[i,tt+1]=Ts[i]+SetChange_T[i]
                    if MLcheck==false
                        adjacentM[:,:,tt+1]=n2
                    end
                    break
                else j=j-1
                end
            end
        end

       
        # println("For ",tt+1," iteration, the xBsetpoint is:", xBsetpoint[:,tt+1])
        xBtvt[tt+1]=sum(adjacentM[i,N+1,tt+1]*flowvt[i,N+1,tt+1]*xBvt[i,tt+1] for i=1:N)/sum(adjacentM[i,N+1,tt+1]*flowvt[i,N+1,tt+1] for i=1:N)
        times[tt+1]=times[tt]+dt
        count=count+1
    end

     # have to reshape because plot accepts a matrix not a vector, also must be 1xN not Nx1
    label = reshape(["R$i" for i in 1:N],1,N)
    if print || save_plots
        p1=plot(times,transpose(T0_invt),xlabel="Time (s)",label=label,ylabel="Input Temperature")
        p2=plot(times,transpose(xBvt),xlabel="Time (s)", label=label,ylabel="Individual xB")
        p3=plot(times,transpose(xBtvt),xlabel="Time (s)", label=false,ylabel="Final Output xB(xB3)")
        p4=plot(times,transpose(heatvt),xlabel="Time (s)", label=label,ylabel="Q (kW)")
        flow_plot=zeros(N,time_steps+1)
        for i=1:N
            flow_plot[i,:]=flowvt[N+1,i,:]
        end
        p5=plot(times,transpose(flow_plot),xlabel="Time (s)", label=label,ylabel="F (m^3/s)")
        p6=plot(times,transpose(Tvt),xlabel="Time (s)",label=label,ylabel="Reactor Temperature")
        p_all=plot(p1,p2,p3,p4,p5,p6,layout=(2,3),legend=:bottomright,xtickfontsize=6,ytickfontsize=6,xguidefontsize=8,yguidefontsize=8)
        if print
            display(p_all)
        end
        
    end

    s = zeros(6)
    b = zeros(count)
    b1 = zeros(count)
    b2 = zeros(count)
    b3 = zeros(count)
    b4 = zeros(count)
    for t = 2:count
        s += [sum(q_xB*(xBtvt[t] - xBsetpoint[end,t])^2), q_T*sum((Tvt[i,t]-Ts[i])^2 for i=1:N), r_flow*sum((flowvt[i,j,t] - flowvt[i,j,t-1])^2 for i=1:N+1 for j=1:N+1),
                r_heat*sum((heatvt[i,t] - heatvt[i,t-1])^2 for i=1:N), 0,0]
        b[t] = b[t-1] + q_xB*sum((xBtvt[t] - xBsetpoint[end,t])^2) + q_T*sum((Tvt[i,t]-Ts[i])^2 for i=1:N) + r_flow*sum((flowvt[i,j,t] - flowvt[i,j,t-1])^2 for i=1:N+1 for j=1:N+1) + r_heat*sum((heatvt[i,t] - heatvt[i,t-1])^2 for i=1:N)
        b1[t] = b1[t-1] + q_xB*sum((xBtvt[t] - xBsetpoint[end,t])^2)
        b2[t] = b2[t-1] + q_T*sum((Tvt[i,t]-Ts[i])^2 for i=1:N)
        b3[t] = b3[t-1] + r_flow*sum((flowvt[i,j,t] - flowvt[i,j,t-1])^2 for i=1:N+1 for j=1:N+1)
        b4[t] = b4[t-1] + r_heat*sum((heatvt[i,t] - heatvt[i,t-1])^2 for i=1:N)
    end
    s[5] = maximum(Tvt[1,:])
    epsilon = 0.01 * xBs[end]
    for i in 1:length(xBtvt)
        if i > dist_time[1] && xBtvt[i] < xBs[1] + epsilon
            s[6] = i
            break
        end

    end

    
    println("writing performance to file")
    # txt file
    # top_file = out_dir * "\\initial_T1_" * string(initial_values[1,2]) *"_T2_" * string(initial_values[2,2]) * "_T3_" * string(initial_values[3,2]) * "_xB1_" *string(initial_values[1,3]) * "_xB2_" *string(initial_values[2,3]) * "_xB3_" *string(initial_values[3,3]) * "_T0_" *string(initial_values[1,1]) * "SetChange_xB_" * string(SetChange_xB[end]) * ".txt"
    # touch(top_file)
    # file = open(top_file, "w")
    
    column_names = ["times","xBset","T01","T02", "T03", "Tvt1","Tvt2","Tvt3", "xBvt1","xBvt2","xBvt3", "xBtvt", "flowvt1", "flowvt2","flowvt3","heatvt1","heatvt2","heatvt3", "Tset1", "Tset2", "Tset3", "Performance index", "xBt PI","Tvt PI","Fvt PI","Qvt PI","tt_stable","Configuration_record"]
    # data=[times,xBsetpoint[end,:],T0_invt[1,:],T0_invt[2,:],T0_invt[3,:],Tvt[1,:],Tvt[2,:],Tvt[3,:],xBvt[1,:],xBvt[2,:],xBvt[3,:],xBtvt,flowvt[N+1,1,:],flowvt[N+1,2,:],flowvt[N+1,3,:],heatvt[1,:],heatvt[2,:],heatvt[3,:],obj_output_total,obj_output_xBt,obj_output_T,obj_output_F,obj_output_Q,fill(s[6],length(times)),record_configuration]
    data=[times,xBsetpoint[end,:],T0_invt[1,:],T0_invt[2,:],T0_invt[3,:],Tvt[1,:],Tvt[2,:],Tvt[3,:],xBvt[1,:],xBvt[2,:],xBvt[3,:],xBtvt,flowvt[N+1,1,:],flowvt[N+1,2,:],flowvt[N+1,3,:],heatvt[1,:],heatvt[2,:],heatvt[3,:],b,b1,b2,b3,b4,fill(s[6],length(times)),record_configuration]

    # write to txt file
    # write(file, join(column_names, "\t") * "\n")
    # writedlm(file, data)
    
    # write to excel file
    # top_excel_file = out_dir * "\\ML_initial_T1_" * string(round(initial_values[1,2];digits=4)) *"_T2_" * string(round(initial_values[2,2];digits=4)) * "_T3_" * string(round(initial_values[3,2];digits=4)) * "_xB1_" *string(round(initial_values[1,3];digits=4)) * "_xB2_" *string(round(initial_values[2,3];digits=4)) * "_xB3_" *string(round(initial_values[3,3];digits=4)) * "_T0_" *string(round(initial_values[1,1]+Dist_T0[1,1];digits=4))* "SetChange_xB_" * string(round(SetChange_xB[end];digits = 4)) * ".xlsx"
    if MLcheck == false
        # top_excel_file = out_dir * "\\noML_initial_T1_" * string(round(initial_values[1,2];digits=4)) *"_T2_" * string(round(initial_values[2,2];digits=4)) * "_T3_" * string(round(initial_values[3,2];digits=4)) * "_xB1_" *string(round(initial_values[1,3];digits=4)) * "_xB2_" *string(round(initial_values[2,3];digits=4)) * "_xB3_" *string(round(initial_values[3,3];digits=4)) * "_Tin1_" *string(round(initial_values[1,1]+Dist_T0[1,2];digits=4))* "_Tin2_" *string(round(initial_values[2,1]+Dist_T0[2,2];digits=4))*  "_Tin3_" *string(round(initial_values[3,1]+Dist_T0[3,2];digits=4))* "SetChange_xB_" * string(round(SetChange_xB[end];digits = 4))  

        top_excel_file = out_dir * "\\noML_initial_T1_" * string(round(initial_values[1,2];digits=4)) *"_T2_" * string(round(initial_values[2,2];digits=4)) * "_T3_" * string(round(initial_values[3,2];digits=4)) * "_xB1_" *string(round(initial_values[1,3];digits=4)) * "_xB2_" *string(round(initial_values[2,3];digits=4)) * "_xB3_" *string(round(initial_values[3,3];digits=4)) * "_Tin1_" *string(round(initial_values[1,1]+Dist_T0[1,2];digits=4))* "_Tin2_" *string(round(initial_values[2,1]+Dist_T0[2,2];digits=4))*  "_Tin3_" *string(round(initial_values[3,1]+Dist_T0[3,2];digits=4))* "SetChange_xB_" * string(round(SetChange_xB[end];digits = 4)) * "SetChange_T1_" * string(round(SetChange_T[1];digits = 4)) * "SetChange_T2_" * string(round(SetChange_T[2];digits = 4)) * "SetChange_T3_" * string(round(SetChange_T[3];digits = 4)) 
    else 
        top_excel_file = out_dir * "\\ML_initial_T1_" * plot_name* string(round(initial_values[1,2];digits=4)) *"_T2_" * string(round(initial_values[2,2];digits=4)) * "_T3_" * string(round(initial_values[3,2];digits=4)) * "_xB1_" *string(round(initial_values[1,3];digits=4)) * "_xB2_" *string(round(initial_values[2,3];digits=4)) * "_xB3_" *string(round(initial_values[3,3];digits=4)) * "_T0_" *string(round(initial_values[1,1]+Dist_T0[1,2];digits=4))* "SetChange_xB_" * string(round(SetChange_xB[end];digits = 4))
    end
    # XLSX.writetable(top_excel_file, data, column_names)
    # close(file)
    
    # DataFrame
    # println(convert(Vector,times))
    df_MPC = DataFrame(times=vec(times), xBset=vec(xBsetpoint[end,:]), T01=vec(T0_invt[1,:]), T02=vec(T0_invt[2,:]), T03=vec(T0_invt[3,:]), T1initial=vec(Tvt[1,:]), T2initial=vec(Tvt[2,:]), T3initial=convert(Vector,Tvt[3,:]), xB1initial=vec(xBvt[1,:]), xB2initial=vec(xBvt[2,:]), xB3initial=vec(xBvt[3,:]), xBtinitial=vec(xBtvt), flowvt1=vec(flowvt[N+1,1,:]), flowvt2=vec(flowvt[N+1,2,:]), flowvt3=vec(flowvt[N+1,3,:]), heatvt1=vec(heatvt[1,:]), heatvt2=vec(heatvt[2,:]), heatvt3=vec(heatvt[3,:]), Tset1=vec(Tsetpoint[1,:]), Tset2=vec(Tsetpoint[2,:]), Tset3=vec(Tsetpoint[3,:]), Performance_index=vec(b), xBt_PI=vec(b1), Tvt_PI=vec(b2), Fvt_PI=vec(b3), Qvt_PI=vec(b4), tt_stable=vec(fill(s[6],length(times))), Configuration_record=vec(record_configuration))
    # df_MPC = DataFrame(data)
    # df_MPC = convert(DataFrame, data)
    CSV.write(top_excel_file* ".csv", df_MPC)
    # println(obj_output_xBt)
    # column_names = ["times","xBset","T01","T02", "T03", "T04","Tvt1","Tvt2","Tvt3","Tvt4", "xBvt1","xBvt2","xBvt3","xBvt4", "xBtvt", "flowvt1", "flowvt2","flowvt3","flowvt4","heatvt1","heatvt2","heatvt3","heatvt4", "ObjValue from MPC", "xBt ISE from MPC","Tvt PI","Fvt PI","Qvt PI","tt_stable"]

    # excel file
    if save_plots
        println("saving fig to $plot_name")
        savefig(top_excel_file*plot_name)
    end
    # data=[times,xBsetpoint[end,:],T0_invt[1,:],T0_invt[2,:],T0_invt[3,:],Tvt[1,:],Tvt[2,:],Tvt[3,:],xBvt[1,:],xBvt[2,:],xBvt[3,:],xBtvt,flowvt[1,N+1,:],flowvt[2,N+1,:],flowvt[3,N+1,:],heatvt[1,:],heatvt[2,:],heatvt[3,:],b,fill(s[6],length(times))]
    # data=[times,xBsetpoint[end,:],T0_invt[1,:],T0_invt[2,:],T0_invt[3,:],Tvt[1,:],Tvt[2,:],Tvt[3,:],xBvt[1,:],xBvt[2,:],xBvt[3,:],xBtvt,flowvt[N+1,1,:],flowvt[N+1,2,:],flowvt[N+1,3,:],heatvt[1,:],heatvt[2,:],heatvt[3,:],b,b1,b2,b3,b4,fill(s[6],length(times))]
    # data=[times,xBsetpoint[end,:],T0_invt[1,:],T0_invt[2,:],T0_invt[3,:],T0_invt[4,:],Tvt[1,:],Tvt[2,:],Tvt[3,:],Tvt[4,:],xBvt[1,:],xBvt[2,:],xBvt[3,:],xBvt[4,:],xBtvt,flowvt[N+1,1,:],flowvt[N+1,2,:],flowvt[N+1,3,:],flowvt[N+1,4,:],heatvt[1,:],heatvt[2,:],heatvt[3,:],heatvt[4,:],b,b1,b2,b3,b4,fill(s[6],length(times))]
    
    

    return s

end


function MPC_step_all(T0_in,T_0,xA_0,xB_0,heat,Flow,n,dt;print=true) # Use one ODE solver to solve the whole system
    
    if print
        println("These are the inputs for MPC_step_all")
        println("T=",T_0)
        println("Tin=",T0_in)
        println("xB=",xB_0)
        println("heat=",heat)
        println("flow=",Flow)
    end
    function odeodes!(du,u,p,t)
        for i=1:N # N reactors in total
            du[3*(i-1)+1] = 1/V[i]*(sum(n[j,i]*Flow[j,i]*u[3*(j-1)+1] for j=1:N) + n[N+1,i]*Flow[N+1,i]*T0_in[i] - sum(n[i,j]*Flow[i,j]*u[3*(i-1)+1] for j=1:N+1)) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/u[3*(i-1)+1])*u[3*(i-1)+2])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/u[3*(i-1)+1])*u[3*(i-1)+3]) + heat[i]/rho/c_p/V[i]  # Temperature of the i th reactor
            du[3*(i-1)+2] = 1/V[i]*(sum(n[j,i]*Flow[j,i]*u[3*(j-1)+2] for j=1:N) + n[N+1,i]*Flow[N+1,i]*xA0 - sum(n[i,j]*Flow[i,j]*u[3*(i-1)+2] for j=1:N+1)) + (-k1*exp(-E1/R_gas/u[3*(i-1)+1])*u[3*(i-1)+2])
            du[3*(i-1)+3] = 1/V[i]*(sum(n[j,i]*Flow[j,i]*u[3*(j-1)+3] for j=1:N) - sum(n[i,j]*Flow[i,j]*u[3*(i-1)+3] for j=1:N+1)) + k1*exp(-E1/R_gas/u[3*(i-1)+1])*u[3*(i-1)+2] + (-k2*exp(-E2/R_gas/u[3*(i-1)+1])*u[3*(i-1)+3])
        end
    end
    # u[3*i+1] T
    # u[3*i+2] xA
    # u[3*i+3] xB

    initial_vec=zeros(3*N)

    for i=0:N-1
        initial_vec[3*i+1]=T_0[i+1]
        initial_vec[3*i+2]=xA_0[i+1]
        initial_vec[3*i+3]=xB_0[i+1]
    end
    if print
        println(initial_vec)
    end

    prob=ODEProblem(odeodes!,initial_vec,(0.0,dt))
    # prob=ODEProblem(odeodes!,initial_vec,(0.0,20*dt))
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    if print
        println("Next measurement is: ", last(soln.u))
    end
    return last(soln.u)
    # a=soln.t
    # A=Array(soln)
    # return a,A[1,:],A[2,:],A[3,:]
    # return a,A[4,:],A[5,:],A[6,:]
    # return the time array and simulated T, xA, xB arrays

end


function findSS_all(T0_in,T_0,xB_0,n;print=true)
    # assume there is no spliting
    # TODO negative flowrate occurs for the mixing reactor with n=[0 0 0 1 0; 0 0 0 1 0; 0 0 0 1 0; 0 0 0 0 1; 1 1 1 1 0]
    # TODO BoundErrors occur if n=[0 0 0 1 0; 0 0 0 1 0; 0 0 0 1 0; 0 0 0 0 1; 1 1 1 0 0]
    # println(T_0)
    Lookup=findall(isone,n) # find all index of open streams
    L=length(Lookup)
    if print
        println("L=",L)
    end
    flow_start=zeros(L)
    flow_start[:].=Ftest
    Ttot=zeros(N+1,N)
    xBtot=zeros(N+1,N)
    for i=1:N+1
        if i!=N+1
            Ttot[i,:].=T_0[i]
            xBtot[i,:].=xB_0[i]
        else
            for j=1:N
                Ttot[i,j]=T0_in[j]
                xBtot[i,j]=0
            end
        end
    end
    # println("FinalT=",Ttot," Length=",length(Ttot))
    heat_start=zeros(N)
    for i=1:N
        heat_start[i] = -rho*c_p*V[i]*(1/V[i]*(sum(flow_start[k]*Ttot[Lookup[k][1],i] for k=1:L if Lookup[k][2]==i) - sum(flow_start[k]*T_0[i] for k=1:L if Lookup[k][1]==i)) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_0[i])*(1-xB_0[i])) + (-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_0[i])*xB_0[i]))
        # println("i=",i)
    end

    # initial_vec2=zeros(4*N*(N+4))
    initial_vec2=zeros(L + N) # flow+heat+xA

    initial_vec2[1:L] = flow_start
    initial_vec2[L+1:end] = heat_start

    function f!(du,u)
        for i=1:N # N reactors in total
            du[3*(i-1)+1] = 1/V[i]*(sum(u[k]*Ttot[Lookup[k][1],i] for k=1:L if Lookup[k][2]==i) - sum(u[k]*T_0[i] for k=1:L if Lookup[k][1]==i)) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_0[i])*(1-xB_0[i]))+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_0[i])*xB_0[i]) + u[L+i]/rho/c_p/V[i]
            du[3*(i-1)+2] = 1/V[i]*(sum(u[k]*xBtot[Lookup[k][1],i] for k=1:L if Lookup[k][2]==i) - sum(u[k]*xB_0[i] for k=1:L if Lookup[k][1]==i)) + k1*exp(-E1/R_gas/T_0[i])*(1-xB_0[i]) + (-k2*exp(-E2/R_gas/T_0[i])*xB_0[i])
            du[3*(i-1)+3] = sum(u[k] for k=1:L if Lookup[k][2]==i) - sum(u[k] for k=1:L if Lookup[k][1]==i)
        end
    end
    # u[1:L] are flow rates
    # u[L+1:L+N] is heating rate
    # u[L+N+1:end] is xA

    soln=nlsolve(f!,initial_vec2)
    heat_ss=zeros(N)
    flow_allconnected=zeros(L)
    # flow_ss=zeros(N)
    flow_ss=soln.zero[1:L]
    heat_ss=soln.zero[L+1:end]
    # for i=1:N
    #     for k=1:L
    #         if Lookup[k][1]==N+1&&Lookup[k][2]==i
    #             # println("i=",i," k=",k)
    #             flow_ss[i]=soln.zero[k]
    #         end
    #     end
    # end
    # println("Lookup=",Lookup)
    # println("flow_ss=",flow_ss," and the length =",length(flow_ss))
    # for i=1:L
    #     if flow_ss[i]<0
    #         println("flow_ss=",flow_ss)
    #         error("Negative flowrate occurs")
    #     end
    # end
    return heat_ss,flow_ss,Lookup
end


# parallel_3R = [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0]
# series_3R = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0]
# parallel_2_and_1_3R = [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0]
# mixing_3R = [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0]
# parallel_4R = [0 0 0 0 1; 0 0 0 0 1; 0 0 0 0 1; 0 0 0 0 1; 1 1 1 1 0]
# non_parallel_4R = [0 0 1 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 1 1 1 0] # just an example, this is 1 and 2 mix into 3 and 4 is in series after 3

# initial_conditions = repeat([300 388.7 0.11],size(parallel_3R)[1] - 1)
# initial_conditions_3R_series = [300 370 0.055;300 380 0.08; 300 388.7 0.11] # 3R series
# initial_conditions_3R_2_and_1 = [300 370 0.055;300 388.7 0.11; 300 388.7 0.11] # 3R 2and1 parallel
# initial_conditions_3R_mixing = [300 370 0.055;300 370 0.055; 300 388.7 0.11] # 3R mixing
# initial_conditions_4R_parallel = repeat([300 388.7 0.11],size(parallel_4R)[1] - 1)
# initial_conditions_4R_non_parallel = repeat([300 388.7 0.11],size(parallel_4R)[1] - 1) # same as above, just an example

# disturbances = [0 0; 0 0; 0 0]

# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\parallel", [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+0 0+40;0+0 0+40;0+0 0+40],[0;0;0],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=1500,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\hybrid", [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+0 0+40;0+0 0+40;0+0 0+40],[0;0;0],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=1500,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\mixing", [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0],[0+0 0+40;0+0 0+40;0+0 0+40],[0;0;0],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=1500,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for AIChE pre\\withoutML\\series", [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0],[0+0 0+40;0+0 0+40;0+0 0+40],[0;0;0],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=1500,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for ESCAPE", [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+0 0+0;0+0 0+0;0+0 0+0],[0;0;0-0.09],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],1,[335 438.7 0.11;335 438.7 0.11;335 438.7 0.11];tmax=3000,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for ESCAPE", [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+0 0+0;0+0 0+0;0+0 0+0],[0;0;0-0.08],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[335 438.7 0.11;335 438.7 0.11;335 438.7 0.11];tmax=3000,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=true)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for ESCAPE", [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0 0 1 0; 0 0 1 0; 0 0 0 1; 1 1 1 0],[0+0 0+0;0+0 0+0;0+0 0+0],[0;0;0-0.09],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[335 438.7 0.11;335 438.7 0.11;335 438.7 0.11];tmax=3000,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for ESCAPE", [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0],[0+0 0+0;0+0 0+0;0+0 0+0],[0;0;0-0.09],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[335 438.7 0.11;335 438.7 0.11;335 438.7 0.11];tmax=3000,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
# MPC_tracking("G:\\My Drive\\Research\\MLReconfiguration\\data for ESCAPE", [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+0 0+0;0+0 0+0;0+0 0+0],[0;0;0-0.09],[0 ;0 ;0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,2],0,[340 388.7 0.19;340 388.7 0.19;340 388.7 0.19];tmax=3000,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=true)

# MPC_tracking("G:\\My Drive\\Class slices\\2024Fall\\IOE618\\project", [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0], [0+0 0+0;0+0 0+20;0+0 0+0],[0;0;0],[0+0 ;0+0 ;0+0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],15,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=3000,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)

MPC_tracking("G:\\My Drive\\Research\\", [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0+40 0+40;0+40 0+40;0+40 0+40],[0;0;0],[0+0 ;0+0 ;0+0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],15,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=1500,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)

# @btime begin 
#     MPC_tracking("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel", [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0], [0 0 0 1; 0 0 0 1; 0 0 0 1; 1 1 1 0],[0+0 0+0;0+0 0+0;0+0 0+20],[0;0;0],[0+0 ;0+0 ;0+0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],15,[300 388.7 0.11;300 388.7 0.11;300 388.7 0.11];tmax=3000,print=false,save_plots=false,plot_name="all_plots.pdf",MLcheck=false)
# end
# MPC_tracking(out_dir, n1::Array{Int,2},n2,Dist_T0,SetChange_xB,SetChange_T,q_T,q_xA,q_xB,r_heat,r_flow,dt,P,
#     dist_time,setpoint_time,initial_values; tmax=200,print=true,save_plots=false,plot_name="all_plots.png",MLcheck=false) # This is for continous disturbance on the (unstable) input temperature
# MPC_tracking("G:\\My Drive\\Research\\GNN projects\\Data\\Parallel", [0 0 1; 0 0 1; 1 1 0], [0 0 1; 0 0 1; 1 1 0],[0+0 0+0; 0 0],[0; 0],[0+0 ;0 + 0 ],1,1e7,1e7,1e-5,1e7,90,1000,[0,15],15,[300 388.7 0.11];tmax=1500,print=false,save_plots=true,plot_name="all_plots.pdf",MLcheck=false)
