#Single Reactor Disturbance Simulation
#https://www.sciencedirect.com/science/article/pii/S000925090800403X?casa_token=aY6Jl0CMNX4AAAAA:JUSu3a5swBkQP8394S3Tfvg0XHZKA5THcWVmWFVhob7QOhQIER3YlNL0F7cW2IbdYC5hzNqg#fig6

#Functions for solving MPC optimization problem and implementing results in receding horizon framework

using Plots, JuMP, Ipopt, DifferentialEquations, NLsolve

function loadProcessData()
    # global F0=9/3600/2 #m^3/s
    global V=0.5/3 #m^3
    global T0=300 #K
    global d_H1=-6e4 #KJ/kmol
    global d_H2=-7e4 #KJ/kmol
    global k1=2.77e3 #s^-1
    global k2=2.5e3 #s^-1
    global E1=5e4 #KJ/kmol
    global E2=6e4 #KJ/kmol
    global c_p=4.2 #KJ/kg/K
    global m=.00279 #kmol/kg
    global rho=1000 #kg/m^3
    global R_gas=8.314 #KJ/kmol/K
    global xA0=1
    global xB0=0
    global T10=388.7
    global T20=388.7
    global xA10=0.89
    global xB10=0.11
    global xA20=xA10
    global xB20=xB10
    global F0=(-k1*exp(-E1/R_gas/T10)*(1-xB10)+(k2*exp(-E2/R_gas/T10)*xB10))*V/(xB0-xB10)
    # global Q1_nom=1.26e6
    # global Q1_nom=rho*c_p*V*(d_H1*m/c_p*k1*exp(-E1/R_gas/388.7)*0.89+d_H2*m/c_p*k2*exp(-E2/R_gas/388.7)*0.11-F0/V*(T0-388.7))
    # global Q2_nom=rho*c_p*V*(d_H1*m/c_p*k1*exp(-E1/R_gas/388.7)*0.89+d_H2*m/c_p*k2*exp(-E2/R_gas/388.7)*0.11-F0/V*(T0-388.7))
    global Q1_nom,ff1=findSS(T0,T10,xB10) #ff1 doesn't matter
    global Q2_nom,ff2=findSS(T0,T10,xB10) #ff2 doesn't matter
    println("Parameters Loaded!")
end

# #Ipopt is a nonlinear programming solver
function MPC_solve(T0_1in,T0_2in,T1_0,xA1_0,xB1_0,T2_0,xA2_0,xB2_0,q_T,q_xA,q_xB,r_heat,r_flow,dt,P;heat_init=0,flow_init=0)
    K=round(Int,P/dt)
    #round(Int,Number)
    MPC=Model(Ipopt.Optimizer)
    MOI.set(MPC, MOI.RawOptimizerAttribute("print_level"), 1)

    heat_ss1,flow_ss1=findSS(T0_1in,T10,xB10)
    heat_ss2,flow_ss2=findSS(T0_2in,T20,xB20)

    println("Q_ss1=",heat_ss1,"flow_ss1=",flow_ss1)

    # xA_guess=zeros(2,K+1)
    # xB_guess=zeros(2,K+1)
    # T_guess=zeros(2,K+1)
    # xA_guess[1,1]=xA1_0
    # xB_guess[1,1]=xB1_0
    # T_guess[1,1]=T1_0
    # xA_guess[2,1]=xA2_0
    # xB_guess[1]=xB2_0
    # T_guess[1]=T2_0
    # xB_tot_guess=zeros(K+1)
    # xA_guess[2,1]=xA2_0
    # xB_guess[2,1]=xB2_0
    # T_guess[2,1]=T2_0

    xA_guess1=zeros(K+1)
    xB_guess1=zeros(K+1)
    T_guess1=zeros(K+1)
    xA_guess2=zeros(K+1)
    xB_guess2=zeros(K+1)
    T_guess2=zeros(K+1)
    xA_guess1[1]=xA1_0
    xB_guess1[1]=xB1_0
    T_guess1[1]=T1_0
    xA_guess2[1]=xA2_0
    xB_guess2[1]=xB2_0
    T_guess2[1]=T2_0
    xB_tot_guess=zeros(K+1)
    xB_tot_guess[1]=(flow_ss1*xB_guess1[1]+flow_ss2*xB_guess2[1])/(flow_ss1+flow_ss2)


    # for k=1:K
    #     xA_guess[1,k+1] = (flow_ss1/V*(xA0-xA_guess[1,k])+(-k1*exp(-E1/R_gas/T_guess[1,k])*xA_guess[1,k]))*dt + xA_guess[1,k]
    #     xB_guess[1,k+1] = (flow_ss1/V*(xB0-xB_guess[1,k])+k1*exp(-E1/R_gas/T_guess[1,k])*xA_guess[1,k]+(-k2*exp(-E2/R_gas/T_guess[1,k])*xB_guess[1,k]))*dt + xB_guess[1,k]
    #     T_guess[1,k+1] = (flow_ss1/V*(T0_1in-T_guess[1,k])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T_guess[1,k])*xA_guess[1,k])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T_guess[1,k])*xB_guess[1,k])+heat_ss1/rho/c_p/V)*dt + T_guess[1,k]
    #     xA_guess[2,k+1] = (flow_ss2/V*(xA0-xA_guess[2,k])+(-k1*exp(-E1/R_gas/T_guess[2,k])*xA_guess[2,k]))*dt + xA_guess[2,k]
    #     xB_guess[2,k+1] = (flow_ss2/V*(xB0-xB_guess[2,k])+k1*exp(-E1/R_gas/T_guess[2,k])*xA_guess[2,k]+(-k2*exp(-E2/R_gas/T_guess[2,k])*xB_guess[2,k]))*dt + xB_guess[2,k]
    #     T_guess[2,k+1] = (flow_ss2/V*(T0_2in-T_guess[2,k])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T_guess[2,k])*xA_guess[2,k])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T_guess[2,k])*xB_guess[2,k])+heat_ss2/rho/c_p/V)*dt + T_guess[2,k]
    #     xB_tot_guess[k+1]=(flow_ss1*xB_guess[1,k+1]+flow_ss2*xB_guess[2,k+1])/(flow_ss1+flow_ss2)
    # end

    for k=1:K
        xA_guess1[k+1] = (flow_ss1/V*(xA0-xA_guess1[k])+(-k1*exp(-E1/R_gas/T_guess1[k])*xA_guess1[k]))*dt + xA_guess1[k]
        xB_guess1[k+1] = (flow_ss1/V*(xB0-xB_guess1[k])+k1*exp(-E1/R_gas/T_guess1[k])*xA_guess1[k]+(-k2*exp(-E2/R_gas/T_guess1[k])*xB_guess1[k]))*dt + xB_guess1[k]
        T_guess1[k+1] = (flow_ss1/V*(T0_1in-T_guess1[k])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T_guess1[k])*xA_guess1[k])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T_guess1[k])*xB_guess1[k])+heat_ss1/rho/c_p/V)*dt + T_guess1[k]
        xA_guess2[k+1] = (flow_ss2/V*(xA0-xA_guess2[k])+(-k1*exp(-E1/R_gas/T_guess2[k])*xA_guess2[k]))*dt + xA_guess2[k]
        xB_guess2[k+1] = (flow_ss2/V*(xB0-xB_guess2[k])+k1*exp(-E1/R_gas/T_guess2[k])*xA_guess2[k]+(-k2*exp(-E2/R_gas/T_guess2[k])*xB_guess2[k]))*dt + xB_guess2[k]
        T_guess2[k+1] = (flow_ss2/V*(T0_2in-T_guess2[k])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T_guess2[k])*xA_guess2[k])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T_guess2[k])*xB_guess2[k])+heat_ss2/rho/c_p/V)*dt + T_guess2[k]
        xB_tot_guess[k+1]=(flow_ss1*xB_guess1[k+1]+flow_ss2*xB_guess2[k+1])/(flow_ss1+flow_ss2)
    end

    JuMP.@variables MPC begin
        Q1[k=0:K-1], (lower_bound=0.2*heat_ss1, upper_bound=1.8*heat_ss1,start=heat_ss1)#Q of the reactor 1
        Q2[k=0:K-1], (lower_bound=0.2*heat_ss2, upper_bound=1.8*heat_ss2,start=heat_ss2)#Q of the reactor 2
        F1[k=0:K-1], (lower_bound=0.2*flow_ss1, upper_bound=1.8*flow_ss1,start=flow_ss1)#flowrate of the feed 1
        F2[k=0:K-1], (lower_bound=0.2*flow_ss2, upper_bound=1.8*flow_ss2,start=flow_ss2)#flowrate of the feed 2

        T1[k=0:K], (lower_bound=T0,upper_bound=2000,start=T_guess1[k+1])
        T2[k=0:K], (lower_bound=T0,upper_bound=2000,start=T_guess2[k+1])
        xA1[k=0:K], (lower_bound=0, upper_bound=1,start=xA_guess1[k+1])
        xA2[k=0:K], (lower_bound=0, upper_bound=1,start=xA_guess2[k+1])
        xB1[k=0:K], (lower_bound=0, upper_bound=1,start=xB_guess1[k+1])
        xB2[k=0:K], (lower_bound=0, upper_bound=1,start=xB_guess2[k+1])
        xB[k=0:K], (lower_bound=0, upper_bound=1,start=xB_tot_guess[k+1])
    end

    @constraints MPC begin
        T1_init, T1[0]==T1_0
        T2_init, T2[0]==T2_0
        xA1_init, xA1[0]==xA1_0
        xA2_init, xA2[0]==xA2_0
        xB1_init, xB1[0]==xB1_0
        xB2_init, xB2[0]==xB2_0
        xB_init, xB[0]==(flow_ss1*xB1_0+flow_ss2*xB2_0)/(flow_ss1+flow_ss2)
    end

#     #Discrete Time for finite differeces
    @NLconstraints MPC begin

        dT1dt[k=0:K-1], T1[k+1]== (F1[k]/V*(T0_1in-T1[k])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T1[k])*xA1[k])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T1[k])*xB1[k])+Q1[k]/rho/c_p/V)*dt + T1[k]
        dT2dt[k=0:K-1], T2[k+1]== (F2[k]/V*(T0_2in-T2[k])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T2[k])*xA2[k])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T2[k])*xB2[k])+Q2[k]/rho/c_p/V)*dt + T2[k]
        dxA1dt[k=0:K-1], xA1[k+1]== (F1[k]/V*(xA0-xA1[k])+(-k1*exp(-E1/R_gas/T1[k])*xA1[k]))*dt + xA1[k]
        dxA2dt[k=0:K-1], xA2[k+1]== (F2[k]/V*(xA0-xA2[k])+(-k1*exp(-E1/R_gas/T2[k])*xA2[k]))*dt + xA2[k]
        dxB1dt[k=0:K-1], xB1[k+1]== (F1[k]/V*(xB0-xB1[k])+k1*exp(-E1/R_gas/T1[k])*xA1[k]+(-k2*exp(-E2/R_gas/T1[k])*xB1[k]))*dt + xB1[k]
        dxB2dt[k=0:K-1], xB2[k+1]== (F2[k]/V*(xB0-xB2[k])+k1*exp(-E1/R_gas/T2[k])*xA2[k]+(-k2*exp(-E2/R_gas/T2[k])*xB2[k]))*dt + xB2[k]
        dxBdt[k=0:K-1], xB[k+1]==(F1[k]*xB1[k+1]+F2[k]*xB2[k+1])/(F1[k]+F2[k])
    end

    # @objective(MPC,Min,q_T*(T1[1]-T10)^2+q_T*(T2[1]-T20)^2+q_xB*(xB[1]-xB10)^2+sum(q_T*(T1[k+1]-T10)^2+q_T*(T2[k+1]-T20)^2+q_xB*(xB[k+1]-xB10)^2+r_heat*(Q1[k]-Q1[k-1])^2+r_heat*(Q2[k]-Q2[k-1])^2+r_flow*(F1[k]-F1[k-1])^2+r_flow*(F2[k]-F2[k-1])^2 for k=1:K-1))
    @objective(MPC,Min,sum(q_T*(T1[k]-T10)^2+q_T*(T2[k]-T20)^2+q_xB*(xB[k]-xB10)^2 for k=0:K)+sum(r_heat*(Q1[k]-Q1[k-1])^2+r_heat*(Q2[k]-Q2[k-1])^2+r_flow*(F1[k]-F1[k-1])^2+r_flow*(F2[k]-F2[k-1])^2 for k=1:K-1)+sum(r_heat*((Q1[0]-Q1[K-1])^2+(Q2[0]-Q2[K-1])^2)+r_flow*((F1[0]-F1[K-1])^2+(F2[0]-F2[K-1])^2)))


    # @objective(MPC,Min,sum(q_T*(T1[k]-T10)^2+q_T*(T2[k]-T20)^2+q_xB*(xB[k]-xB10)^2+r_heat*(Q1[k-1]-heat_ss1)^2+r_heat*(Q2[k-1]-heat_ss2)^2
    #     +r_flow*(F1[k-1]-flow_ss1)^2+r_flow*(F2[k-1]-flow_ss2)^2 for k=1:K))

    JuMP.optimize!(MPC)

    heat_soln1=JuMP.value.(Q1)
    heat_soln2=JuMP.value.(Q2)
    flow_soln1=JuMP.value.(F1)
    flow_soln2=JuMP.value.(F2)
    xA_soln1=JuMP.value.(xA1)
    xA_soln2=JuMP.value.(xA2)
    xB_soln1=JuMP.value.(xB1)
    xB_soln2=JuMP.value.(xB2)
    println("xB_soln1=",xB_soln1)
    println("xB_soln2=",xB_soln2)
    xB_soln=JuMP.value.(xB)
    println("xB_soln=",xB_soln)
    T_soln1=JuMP.value.(T1)
    T_soln2=JuMP.value.(T2)

    obj_T=sum(q_T*((T_soln1[k]-T10)^2+(T_soln2[k]-T20)^2) for k=0:K)
    obj_xBt=sum(q_xB*(xB_soln[k]-xB10)^2 for k=0:K)
    obj_Q=sum(r_heat*((heat_soln1[k]-heat_soln1[k-1])^2+(heat_soln2[k]-heat_soln2[k-1])^2) for k=1:K-1)+sum(r_heat*((heat_soln1[0]-heat_soln1[K-1])^2+(heat_soln2[0]-heat_soln2[K-1])^2))
    obj_F=sum(r_flow*((flow_soln1[k]-flow_soln1[k-1])^2+(flow_soln2[k]-flow_soln2[k-1])^2) for k=1:K-1)+sum(r_flow*((flow_soln1[0]-flow_soln1[K-1])^2+(flow_soln2[0]-flow_soln2[K-1])^2))
    println("Obj_T= ",obj_T)
    println("Obj_xBt= ",obj_xBt)
    println("Obj_Q= ",obj_Q)
    println("Obj_F= ",obj_F)
    st=raw_status(MPC)
    if st=="Infeasible_Problem_Detected"
        println(xA1_guess,xA2_guess)
        error("Solver infeasible, problem stopping")
    end
    # The error function throws a generic ErrorException. This will interrupt execution of the function or block immediately.
    # receive a solver specific string explaning why the optimization stopped
    obj=objective_value(MPC)
    println("Q1soln=",heat_soln1[0],"Q2soln=",heat_soln2[0], " fsoln1=",flow_soln1[0], " fsoln2=",flow_soln2[0]," Tsoln1=",T_soln1[0], " Tsoln1=",T_soln2[0]," xA1=",xA_soln1[0]," xA2=",xA_soln2[0]," obj=",obj)

    return heat_soln1[0], flow_soln1[0], heat_soln2[0], flow_soln2[0], xA_soln1, xA_soln2, T_soln1, T_soln2

end

function MPC_tracking(Dist_T10,Dist_T20,q_T,q_xA,q_xB,r_heat,r_flow,dt,P,dist_time;tmax=200)
    #(runs the moving horizon loop for set point tracking)
    loadProcessData()
    time_steps=round(Int,tmax/dt)
    times=zeros(time_steps+1)

    T0_1invt=zeros(time_steps+1)
    T1vt=zeros(time_steps+1)
    xA1vt=zeros(time_steps+1)
    xB1vt=zeros(time_steps+1)
    heat1vt=zeros(time_steps+1)
    flow1vt=zeros(time_steps+1)
    # For reactor 1
    T0_2invt=zeros(time_steps+1)
    T2vt=zeros(time_steps+1)
    xA2vt=zeros(time_steps+1)
    xB2vt=zeros(time_steps+1)
    heat2vt=zeros(time_steps+1)
    flow2vt=zeros(time_steps+1)
    xBvt=zeros(time_steps+1)
    # For reactor 2

    T0_1invt[1]=T0
    T1vt[1]=T10
    xA1vt[1]=xA10
    xB1vt[1]=xB10
    heat1vt[1]=Q1_nom
    flow1vt[1]=F0

    T0_2invt[1]=T0
    T2vt[1]=T20
    xA2vt[1]=xA20
    xB2vt[1]=xB20
    heat2vt[1]=Q2_nom
    flow2vt[1]=F0
    xBvt[1]=(xB1vt[1]*flow1vt[1]+xB2vt[1]*flow2vt[1])/(flow1vt[1]+flow2vt[1])
    println("xBvt[1]=", xBvt[1])

    times[1]=0
    tt=1
    for tt=1:time_steps
        println("These are the inputs for MPC_solve")
        println("T=",T1vt[tt]," and ", T2vt[tt])
        println("Tin=",T0_1invt[tt]," and ", T0_2invt[tt])
        println("xB=",xB1vt[tt]," and ",xB2vt[tt])
        println("heat=",heat1vt[tt]," and ",heat2vt[tt])
        println("flow=",flow1vt[tt]," and ",flow2vt[tt])

        heat1vt[tt+1],flow1vt[tt+1],heat2vt[tt+1],flow2vt[tt+1]=MPC_solve(T0_1invt[tt],T0_2invt[tt],T1vt[tt],xA1vt[tt],xB1vt[tt],T2vt[tt],xA2vt[tt],xB2vt[tt],q_T,q_xA,q_xB,r_heat,r_flow,dt,P;heat_init=heat1vt[tt],flow_init=flow1vt[tt])
        if tt>=dist_time
            T0_1invt[tt]=T0+Dist_T10
            T0_2invt[tt]=T0+Dist_T20
        end
        #Solve the optimal flowvt and heatvt by the last dt measurements
        newstate1=MPC_step(T0_1invt[tt],T1vt[tt],xA1vt[tt],xB1vt[tt],heat1vt[tt+1],flow1vt[tt+1],dt)
        newstate2=MPC_step(T0_2invt[tt],T2vt[tt],xA2vt[tt],xB2vt[tt],heat2vt[tt+1],flow2vt[tt+1],dt)
        #getting the new measurements
        if tt>=dist_time
            T0_1invt[tt+1]=T0+Dist_T10
            T0_2invt[tt+1]=T0+Dist_T20
        else
            T0_1invt[tt+1]=T0
            T0_2invt[tt+1]=T0
        end
        T1vt[tt+1]=newstate1[1]
        xA1vt[tt+1]=newstate1[2]
        xB1vt[tt+1]=newstate1[3]
        T2vt[tt+1]=newstate2[1]
        xA2vt[tt+1]=newstate2[2]
        xB2vt[tt+1]=newstate2[3]
        xBvt[tt+1]=(flow1vt[tt+1]*xB1vt[tt+1]+flow2vt[tt+1]*xB2vt[tt+1])/(flow1vt[tt+1]+flow2vt[tt+1])
        times[tt+1]=times[tt]+dt
        println(xBvt[tt+1])
    end

    l=@layout [
        grid(3,1){0.2w} grid(4,2)
    ]
    plot(times,[xB1vt,xB2vt,xBvt,T1vt,T2vt,T0_1invt,T0_2invt,heat1vt,heat2vt,flow1vt,flow2vt],layout=l,legend=false,xlabel="Time (s)", ylabel=["xB1" "xB2" "xB" "Temp1(K)" "Temp2(K)" "T1_Input(K)" "T2_Input(K)" "Q1(KW)" "Q2(KW)" "F1(m^3/s)" "F2(m^3/s)"],gridalpha=0.75,minorgrid=true,minorgridalpha=0.1)
    # plot(times,[T1vt,xA1vt,xB1vt,heat1vt,flow1vt,T2vt,xA2vt,xB2vt,heat2vt,flow2vt],layout=(2,5),legend=false,xlabel="Time (s)", ylabel=["Temp1(K)" "xA1" "xB1" "Q1(KW)" "F1(m^3/s)" "Temp2(K)" "xA2" "xB2" "Q2(KW)" "F2(m^3/s)"],
    #     gridalpha=0.75,minorgrid=true,minorgridalpha=0.3)
    top_excel_file = out_dir * "\\2Rcode_T1_" * string(388.7) *"_T2_" * string(388.7) * "_T3_" * string(388.7) * "_xB1_" *string(0.11) * "_xB2_" *string(0.11) * "_xB3_" *string(0.11) * "_Tin1_" *string(300)* "_Tin2_" *string(300)*  "_Tin3_" *string(300)* "SetChange_xB_" * string(0.0) * "SetChange_T1_" * string(388.7) * "SetChange_T2_" * string(388.7) * "SetChange_T3_" * string(388.7) 
    # R2=R3
    df_MPC = DataFrame(times=vec(times), T01=vec(T0_1invt), T02=vec(T0_2invt), T03=vec(T0_2invt), T1initial=vec(T1vt), T2initial=vec(T2vt), T3initial=vec(T2vt), xB1initial=vec(xB1vt), xB2initial=vec(xB2vt), xB3initial=vec(xB2vt), xBtinitial=vec(xBvt), flowvt1=vec(flow1vt), flowvt2=vec(flow2vt), flowvt3=vec(flow2vt), heatvt1=vec(heat1vt), heatvt2=vec(heat2vt), heatvt3=vec(heat2vt), Performance_index=vec(b), xBt_PI=vec(b1), Tvt_PI=vec(b2), Fvt_PI=vec(b3), Qvt_PI=vec(b4), tt_stable=vec(fill(s[6],length(times))), Configuration_record=vec(record_configuration))

    CSV.write(top_excel_file* ".csv", df_MPC)
end

function MPC_step(T0_in,T_0,xA_0,xB_0,heat,flow,dt)
    # println("These are the inputs for MPC_step_all")
    # println("T=",T_0)
    # println("Tin=",T0_in)
    # println("xB=",xB_0)
    # println("heat=",heat)
    # println("flow=",flow)
    #Getting the measurement C and T, (steps forward in time by dt and takes new "measurements" by simulating the system with the decided control input)
    #Usually f(t,u) or in-place f(t,u,du),the latter is used here.
    f(y,p,t)=[flow/V*(T0_in-y[1])+(-d_H1*m/c_p*k1*exp(-E1/R_gas/y[1])*y[2])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/y[1])*y[3])+heat/rho/c_p/V,
        flow/V*(xA0-y[2])+(-k1*exp(-E1/R_gas/y[1])*y[2]),
        flow/V*(xB0-y[3])+k1*exp(-E1/R_gas/y[1])*y[2]+(-k2*exp(-E2/R_gas/y[1])*y[3])]
    #f definition: y[1]is the Temperature, y[2] is the fraction of xA, y[3] is the fraction for B)

    prob=ODEProblem(f,[T_0,xA_0,xB_0],(0.0,dt))
    # prob=ODEProblem(f,[T_0,xA_0,xB_0],(0.0,20*dt))
    #ODEProblem is a function from DifferentialEquations package (f,boundary,tspan) https://diffeq.sciml.ai/stable/tutorials/ode_example/
    #ODEProblem is defining a problem, to solve it, use solve() function
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    #For stiff problems art high tolerances, Rosenbrock23() is more efficient for small systems, TRBDF2 is better for large systems
    # println("Next measurement is: ", last(soln.u))
    return last(soln.u)
    # a=soln.t
    # A=Array(soln)
    # return a,A[1,:],A[2,:],A[3,:]
    # last() means get the last element of an ordered collection
    # println("T_step=",soln[1],"xA_step=",soln[2])
end

function findSS(T0_in,T10,xB10)
    println("Here is the function of calculating steady state")
    println(T0_in," and ", T10," and ",xB10)
    #This function solves the nonlinear system of equations to get a steady state input for a desired T/xB set point.
    flow_start=(-k1*exp(-E1/R_gas/T10)*(1-xB10)+(k2*exp(-E2/R_gas/T10)*xB10))*V/(-xB10)
    heat_start=-rho*c_p*V*(flow_start/V*(T0_in-T10)+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T10)*(1-xB10))+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T10)*xB10))
    # println("flow_start=",flow_start)
    # println("heat_start=",heat_start)
    function f!(F,x)
        F[1]=x[1]/V*(T0_in-T10)+(-d_H1*m/c_p*k1*exp(-E1/R_gas/T10)*x[3])+(-d_H2*m/c_p*k2*exp(-E2/R_gas/T10)*xB10)+x[2]/rho/c_p/V
        F[2]=x[1]/V*(xA0-x[3])+(-k1*exp(-E1/R_gas/T10)*x[3])
        F[3]=(x[1]/V*(xB0-xB10)+k1*exp(-E1/R_gas/T10)*x[3]+(-k2*exp(-E2/R_gas/T10)*xB10))
    end
    # Here assume xA10+xB10=1
    soln=nlsolve(f!,[flow_start,heat_start,1-xB10])
    xA_ss=soln.zero[3]
    heat_ss=soln.zero[2]
    flow_ss=soln.zero[1]
    return heat_ss,flow_ss,xA_ss
end
global out_dir = "G:\\My Drive\\Research\\GNN projects\\Data\\Parallel"
MPC_tracking(0,0,1,1e7,1e7,1e-5,1e7,90,1000,0;tmax=1500)