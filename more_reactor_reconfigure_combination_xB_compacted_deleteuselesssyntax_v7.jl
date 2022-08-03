# For any given number of reactors and potential configurations
# https://www.sciencedirect.com/science/article/pii/S000925090800503X?casa_token=aY6Jl0CMNX5AAAAA:JUSu3a5swBkQP8395S3Tfvg0XHZKA5THcWVmWFVhob7QOhQIER3YlNL0F7cW2IbdYC5hzNqg#fig6

using Plots, JuMP, DifferentialEquations, NLsolve, BenchmarkTools, Ipopt, MathOptInterface, Printf, ProgressBars, DelimitedFiles, Profile

function loadProcessData(N::Int,n::Array{Int,2})
    # global F0=9/3600/N #m^3/s
    global V=0.5/3 #m^3
    # global T0=[300 300 300] #K
    global T0=[300 300] #K
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
    # global Ts=[370;370;388.7] # will change with different input n and other initial conditions
    # global xBs=[0.055; 0.055; 0.11] # will change with different input n and other initial conditions
    global Ts=[388.7;388.7] # will change with different input n and other initial conditions
    global xBs=[0.11;0.11] # will change with different input n and other initial conditions
    # global xAs=[1-xBs[1];1-xBs[2];1-xBs[3]] # will change with different input n and other initial conditions
    global xAs=[1-xBs[1];1-xBs[2]] # will change with different input n and other initial conditions
    global F0=(-k1*exp(-E1/R_gas/Ts[1])*(1-xBs[1])+(k2*exp(-E2/R_gas/Ts[1])*xBs[1]))*V/(xB0-xBs[1])
    global Flow0=zeros(N,N+2)
    for i=1:N
        for j=1:N+2
            if n[i,j]==1
                Flow0[i,j]=F0
            end
        end
    end
    global Q_nom=zeros(N)
    global F_nom=zeros(N)
    Q_nom,F_nom=findSS_all(T0,Ts,xBs,n,Flow0)
    println("Parameters Loaded!")
    # println("Q = ", Q_nom," F = ", F_nom)
end

function MPC_solve(n::Array{Int,2},Flow,T0_inreal,T_0real,xA_0real,xB_0real,q_T,q_xA,q_xB,r_heat,r_flow,dt,P,N;heat_init=0,flow_init=0)
    # println("n=",n)
    global count

    K=round(Int,P/dt)

    MPC=Model(Ipopt.Optimizer)
    MOI.set(MPC, MOI.RawOptimizerAttribute("print_level"), 1)

    T0_in=T0_inreal
    T_0=T_0real
    xA_0=xA_0real
    xB_0=xB_0real

    # println("T0_in=",T0_in)
    # println("T_0=",T_0)
    # println("xA_0=",xA_0)
    # println("xB_0=",xB_0)

    heat_ss=zeros(N)
    flow_ss=zeros(N)

    # Only steady states for input streams from outside instead of other reactors
    # heat_ss,flow_ss=findSS_all(T0_in,Ts,xBs,n,Flow)
    heat_ss,flow_ss=findSS_all(T0_in,Ts,xBs,n,Flow0)

    # println("Heat_ss=",heat_ss)
    # println("Flow_ss=",flow_ss)

    for i=1:N
        Flow[i,N+1]=flow_ss[i]
    end
    # println("Flow=",Flow)

    xA_guess=zeros(N,K+1)
    xB_guess=zeros(N,K+1)
    T_guess=zeros(N,K+1)
    for i=1:N
        xA_guess[i,1]=xA_0[i]
        xB_guess[i,1]=xB_0[i]
        T_guess[i,1]=T_0[i]
    end
    xB_tot_guess=zeros(K+1)
    xB_tot_guess[1]=sum(n[i,N+2]*Flow0[i,N+2]*xB_guess[i,1] for i=1:N)/sum(n[i,N+2]*Flow0[i,N+2] for i=1:N)

    for k=1:K
        for i=1:N
            T_guess[i,k+1] = (1/V*(sum(n[j,i]*Flow0[j,i]*T_guess[j,k] for j=1:N) + n[i,N+1]*Flow0[i,N+1]*T0_in[i]- sum(n[i,j]*Flow0[i,j]*T_guess[i,k] for j=1:N) - n[i,N+2]*Flow0[i,N+2]*T_guess[i,k]) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_guess[i,k])*xA_guess[i,k])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_guess[i,k])*xB_guess[i,k]) + heat_ss[i]/rho/c_p/V)*dt + T_guess[i,k]
            # println("When k=",k+1," T_guess=",T_guess[:,k+1])
            xA_guess[i,k+1] = (1/V*(sum(n[j,i]*Flow0[j,i]*xA_guess[j,k] for j=1:N) + n[i,N+1]*Flow0[i,N+1]*xA0 - sum(n[i,j]*Flow0[i,j]*xA_guess[i,k] for j=1:N) - n[i,N+2]*Flow0[i,N+2]*xA_guess[i,k]) + (-k1*exp(-E1/R_gas/T_guess[i,k])*xA_guess[i,k]))*dt + xA_guess[i,k]
            xB_guess[i,k+1] = (1/V*(sum(n[j,i]*Flow0[j,i]*xB_guess[j,k] for j=1:N) - sum(n[i,j]*Flow0[i,j]*xB_guess[i,k] for j=1:N) - n[i,N+2]*Flow0[i,N+2]*xB_guess[i,k]) + k1*exp(-E1/R_gas/T_guess[i,k])*xA_guess[i,k] + (-k2*exp(-E2/R_gas/T_guess[i,k])*xB_guess[i,k]))*dt + xB_guess[i,k]
            # println("When k=",k+1," xB_guess=",xB_guess[:,k+1])
        end
        xB_tot_guess[k+1] = sum(n[i,N+2]*Flow0[i,N+2]*xB_guess[i,k] for i=1:N)/sum(n[i,N+2]*Flow0[i,N+2] for i=1:N)
    end
    # println("n=",n[:,N+2])
    # println("Flow[i,N+2]=",Flow[:,N+2])
    # println("xB_guess=",xB_guess)
    # println("xBt_guess=",xB_tot_guess)

    JuMP.@variables MPC begin
        # Q[i=1:N,k=0:K-1], (lower_bound=0.2*heat_ss[i], upper_bound=1.8*heat_ss[i],start=heat_ss[i])# Q of the reactors
        Q[i=1:N,k=0:K-1], (lower_bound=0, upper_bound=1.8*heat_ss[i],start=heat_ss[i])# Q of the reactors
        F[i=1:N,j=1:N+2,k=0:K-1], (lower_bound=n[i,j]*0.2*flow_ss[i], upper_bound=n[i,j]*1.8*flow_ss[i],start=Flow[i,j])# Flowrate between reactors
        # F[i=1:N,j=1:N+2,k=0:K-1;n[i,j]==1], (lower_bound=0.2*flow_ss[i], upper_bound=1.8*flow_ss[i],start=flow_ss[i])# Flowrate between reactors
        # F[i=1:N,j=1:N+2,k=0:K-1];n[i,j]==0, (lower_bound=0, upper_bound=1e-5,start=0)# Flowrate between reactors

        T[i=1:N,k=0:K], (lower_bound=T0[i],upper_bound=2000,start=T_guess[i,k+1])
        xA[i=1:N,k=0:K], (lower_bound=0, upper_bound=1,start=xA_guess[i,k+1])
        xB[i=1:N,k=0:K], (lower_bound=0, upper_bound=1,start=xB_guess[i,k+1])
        xBt[k=0:K], (lower_bound=0, upper_bound=1,start=xB_tot_guess[k+1])

        # m[k=0:K], (lower_bound=0, upper_bound=1,start=m_init)
    end

    @constraints MPC begin
        T_init[i=1:N], T[i,0]==T_0[i]
        xA_init[i=1:N], xA[i,0]==xA_0[i]
        xB_init[i=1:N], xB[i,0]==xB_0[i]
        xBt_init, xBt[0]==sum(n[i,N+2]*Flow[i,N+2]*xB_0[i] for i=1:N)/sum(n[i,N+2]*Flow[i,N+2] for i=1:N)
        MassB[i=1:N,k=0:K-1], sum(n[j,i]*F[j,i,k] for j=1:N) + n[i,N+1]*F[i,N+1,k]==sum(n[i,j]*F[i,j,k] for j=1:N) + n[i,N+2]*F[i,N+2,k]
    end

    @NLconstraints MPC begin
        Temp[i=1:N,k=0:K-1], T[i,k+1] == (1/V*(sum(n[j,i]*F[j,i,k]*T[j,k] for j=1:N) + n[i,N+1]*F[i,N+1,k]*T0_in[i]- sum(n[i,j]*F[i,j,k]*T[i,k] for j=1:N) - n[i,N+2]*F[i,N+2,k]*T[i,k]) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T[i,k])*xA[i,k])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T[i,k])*xB[i,k]) + Q[i,k]/rho/c_p/V)*dt + T[i,k]
        MoleFractionxA[i=1:N,k=0:K-1], xA[i,k+1] == (1/V*(sum(n[j,i]*F[j,i,k]*xA[j,k] for j=1:N) + n[i,N+1]*F[i,N+1,k]*xA0 - sum(n[i,j]*F[i,j,k]*xA[i,k] for j=1:N) - n[i,N+2]*F[i,N+2,k]*xA[i,k]) + (-k1*exp(-E1/R_gas/T[i,k])*xA[i,k]))*dt + xA[i,k]
        MoleFractionxB[i=1:N,k=0:K-1], xB[i,k+1] == (1/V*(sum(n[j,i]*F[j,i,k]*xB[j,k] for j=1:N) - sum(n[i,j]*F[i,j,k]*xB[i,k] for j=1:N) - n[i,N+2]*F[i,N+2,k]*xB[i,k]) + k1*exp(-E1/R_gas/T[i,k])*xA[i,k] + (-k2*exp(-E2/R_gas/T[i,k])*xB[i,k]))*dt + xB[i,k]
        # OutputMoleFraction[k=0:K-1], xBt[k+1] == sum(n[i,N+2]*F[i,N+2,k]*xB[i,k] for i=1:N)/sum(n[i,N+2]*F[i,N+2,k] for i=1:N)
        OutputMoleFraction[k=0:K-1], xBt[k+1] == sum(n[i,N+2]*F[i,N+2,k]*xB[i,k+1] for i=1:N)/sum(n[i,N+2]*F[i,N+2,k] for i=1:N)
        # OutputMoleFraction[k=1:K], xBt[k] == sum(n[i,N+2]*F[i,N+2,k]*xB[i,k] for i=1:N)/sum(n[i,N+2]*F[i,N+2,k] for i=1:N)
    end

    @objective(MPC,Min,sum(q_T*(T[i,k]-Ts[i])^2 for i=1:N for k=0:K)+sum(q_xB*(xBt[k]-xBs[1])^2 for k=0:K)+sum(r_heat*(Q[i,k]-Q[i,k-1])^2+r_flow*(n[i,j]*F[i,j,k]-n[i,j]*F[i,j,k-1])^2 for i=1:N for j=1:N+2 for k=1:K-1)+sum(r_heat*(Q[i,0]-Q[i,K-1])^2+r_flow*(n[i,j]*F[i,j,0]-n[i,j]*F[i,j,K-1])^2 for i=1:N for j=1:N+2))
    JuMP.optimize!(MPC)


    st=MathOptInterface.RawStatusString()
    if st=="INFEASIBLE_POINT"
        println(xA1_guess,xA2_guess)
        error("Solver infeasible, problem stopping")
    end
    # obj=getobjectivevalue(MPC) # works for Julia 1.15.3
    obj=objective_value(MPC) # works for Julia 1.17.2
    println("Obj in MPC=",obj)

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

    obj_T=sum(q_T*(results_T[i,k]-Ts[i])^2 for i=1:N for k=0:K)
    obj_xBt=sum(q_xB*(results_xBt[k]-xBs[1])^2 for k=0:K)
    obj_Q=sum(r_heat*(results_heat[i,k]-results_heat[i,k-1])^2 for i=1:N for k=1:K-1)+sum(r_heat*(results_heat[i,0]-results_heat[i,K-1])^2 for i=1:N)
    obj_F=sum(r_flow*(n[i,j]*results_flow[i,j,k]-n[i,j]*results_flow[i,j,k-1])^2 for i=1:N for j=1:N+2 for k=1:K-1)+sum(r_flow*(n[i,j]*results_flow[i,j,0]-n[i,j]*results_flow[i,j,K-1])^2 for i=1:N for j=1:N+2)

    # println("results_heat=",results_heat)
    # println("results_flow=",results_flow)
    # println("Obj_T= ",obj_T)
    # println("Obj_xBt= ",obj_xBt)
    # println("Obj_Q= ",obj_Q)
    # println("Obj_F= ",obj_F)

    println("soln_heat=",results_heat0)
    println("soln_flow=",results_flow0)
    return results_heat0, results_flow0

end

function MPC_tracking(n::Array{Int,2},Dist_T0,q_T,q_xA,q_xB,r_heat,r_flow,dt,P,
    dist_time;tmax=200,print=true,save_plots=false,plot_name="all_plots.png") # This is for continous disturbance on the (unstable) input temperature
    # (runs the moving horizon loop for set point tracking)
    # N=length(Dist_T0)
    # When testing continous disturbance system, the Dist_T0 contains the beginning point
    global N=size(n)[1]
    l=length(dist_time)
    # Check the length of disturbance vectors and dist_time vector are the same
    if size(Dist_T0)[2] == l
        println("The length of disturbance variables == the one of dist_time vector=",l)
    else println("The length of disturbance variables are not equal to the one of dist_time vector")
        return
    end
    loadProcessData(N,n)

    time_steps=round(Int,tmax/dt)

    global times=zeros(time_steps+1)
    pt=round(Int,P/dt)

    global R3IndiObj=zeros(N,time_steps)
    global MPCSolus=zeros(N,pt)
    global heatss=zeros(N,time_steps+1)
    global flowss=zeros(N,time_steps+1)
    global count=1

    global ObjValue=zeros(time_steps) # To storage the optimal objective value from MPC for each iteration

    global Tvt=zeros(N,time_steps+1)
    global xAvt=zeros(N,time_steps+1)
    global xBvt=zeros(N,time_steps+1)
    global heatvt=zeros(N,time_steps+1)
    global flowvt=zeros(N,N+2,time_steps+1)
    global xBtvt=zeros(1,time_steps+1)
    newstate=zeros(N,N)
    # Y=zeros(time_steps+1)
    global T0_invt=zeros(N,time_steps+1)

    # global recordFindSS=zeros()
    # global recordStepAll=zeros()

    T0_invt[:,1]=T0
    Tvt[:,1]=Ts[:]
    xAvt[:,1]=xAs[:]
    xBvt[:,1]=xBs[:]
    heatvt[:,1]=Q_nom
    flowvt[:,:,1]=Flow0
    xBtvt[1]=sum(n[i,N+2]*flowvt[i,N+2,1]*xBvt[i,1] for i=1:N)/sum(n[i,N+2]*flowvt[i,N+2,1] for i=1:N)
    times[1]=0
    tt=1

    for tt=1:time_steps
        # println("These are the inputs for MPC_solve")
        # println("T=",Tvt[:,tt])
        # println("Tin=",T0_invt[:,tt])
        # println("xB=",xBvt[:,tt])
        # println("heat=",heatvt[:,tt])
        # println("flow=",flowvt[:,:,tt])
        resultsheatvt,resultsflowvt=MPC_solve(n,flowvt[:,:,tt],T0_invt[:,tt],Tvt[:,tt],xAvt[:,tt],xBvt[:,tt],q_T,q_xA,q_xB,r_heat,r_flow,dt,P,N;heat_init=heatvt[1,tt],flow_init=flowvt[1,1,tt])
        for i=1:N
            heatvt[i,tt+1]=resultsheatvt[i]
            flowvt[i,:,tt+1]=resultsflowvt[i,:]
        end
        println("count(tt)=",count)

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

        newstate=MPC_step_all(T0_invt[:,tt],Tvt[:,tt],xAvt[:,tt],xBvt[:,tt],heatvt[:,tt+1],flowvt[:,:,tt+1],n,dt)
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

        xBtvt[tt+1]=sum(n[i,N+2]*flowvt[i,N+2,tt+1]*xBvt[i,tt+1] for i=1:N)/sum(n[i,N+2]*flowvt[i,N+2,tt+1] for i=1:N)
        times[tt+1]=times[tt]+dt
        count=count+1
    end

    # p[1]=plot(times,[T0_invt[1,:],xBvt[1,:],xBtvt,heatvt[1,:],flowvt[1,:],mvt],layout=(2,3),xlabel="Time (s)",label=["R1" "R1" false "R1" "R1" false], ylabel=["T_Input(K)" "xB" "Final output xB" "Q(KW)" "F1(m^3/s)" "Flow assignment "])
    # for i=2:N
    #     R=@sprintf("R%s",i)
    #     p[i]=plot!(times,[T0_invt[i,:],xBvt[i,:],heatvt[i,:],flowvt[i,:],xBtvt],layout=(2,3),label=[R R false R R false])
    # end
    # for i=1:N
    #     display(p[i])
    # end
    if print || save_plots
        p1=plot(times,transpose(T0_invt),xlabel="Time (s)",label=["R1" "R2" "R3"],ylabel="Input Temperature")
        p2=plot(times,transpose(xBvt),xlabel="Time (s)", label=["R1" "R2" "R3"],ylabel="Individual xB")
        p3=plot(times,transpose(xBtvt),xlabel="Time (s)", label=false,ylabel="Final Output xB(xB3)")
        p4=plot(times,transpose(heatvt),xlabel="Time (s)", label=["R1" "R2" "R3"],ylabel="Q (kW)")
        flow_plot=zeros(N,time_steps+1)
        for i=1:N
            flow_plot[i,:]=flowvt[i,N+1,:]
        end
        p5=plot(times,transpose(flow_plot),xlabel="Time (s)", label=["R1" "R2" "R3"],ylabel="F (m^3/s)")
        p6=plot(times,transpose(Tvt),xlabel="Time (s)",label=["R1" "R2" "R3"],ylabel="Reactor Temperature")
        p_all=plot(p1,p2,p3,p4,p5,p6,layout=(2,3),xtickfontsize=6,ytickfontsize=6,xguidefontsize=8,yguidefontsize=8)
        if print
            display(p_all)
        end
        if save_plots
            println("saving fig to $plot_name")
            savefig(plot_name)
        end
    end

    s = zeros(4)
    for t = 2:count
        s += [sum((xBtvt[t] - xBs[1])^2), sum((Tvt[i,t]-Ts[i])^2 for i=1:N), sum((flowvt[i,j,t] - flowvt[i,j,t-1])^2 for i=1:N for j=1:N+2),
                sum((heatvt[i,t] - heatvt[i,t-1])^2 for i=1:N)]
    end
    return s
    # savefig("3R_1P2S_xBs_0.055_0.055_0.11_DistT3_10_time_8_allprofiles.png")
    # savefig("3R_1P2S_xBs_0.08_0.08_0.11_NoDist_allprofiles.png")


    # l=@layout [
    #     grid(1,1){0.5w} grid(N,1)]
    #
    # p7=plot(times,transpose(xBtvt),xlabel="Time (s)",label=false,ylabel="Final Output xB(xB3)")
    # p8=plot(times,xBvt[1,:],xlabel=false,label=false,ylabel="xB1")
    # p9=plot(times,xBvt[2,:],xlabel=false,label=false,ylabel="xB2")
    # p10=plot(times,xBvt[3,:],xlabel="Time (s)",label=false,ylabel="xB3")
    # p_con=plot(p7,p8,p9,p10,layout=l,legend=false)
    # display(p_con)

    # ylabel_array=["Final output xB"]
    # label_array=["Mixer"]
    # Y=xBtvt
    # for i=1:N
    #     Y=cat(Y,xBvt[i,:],dims=(2,2))
    #     aa=[@sprintf("xB%s",i)]
    #     bb=[@sprintf("R%s",i)]
    #     ylabel_array=append!(ylabel_array,aa)
    #     label_array=append!(label_array,bb)
    # end
    # # println("ylabel=",ylabel_array)
    # p_con=plot(times,Y,layout=l,xlabel="Time (s)", ylabel=reshape(ylabel_array,1,N+1),label=reshape(label_array,1,N+1))

    # savefig("3R_1P2S_xBs_0.1_0.1_0.1_DistT3_10_time_4_Only_xB.png")
    # savefig("3R_1P2S_xBs_0.1_0.1_0.1_NoDist_Only_xB.png")

    # println("Objective value from MPC=")
    # show(stdout, "text/plain", ObjValue)
    # println("Real system objective value B")
    # show(stdout, "text/plain", b)
end

function MPC_step_all(T0_in,T_0,xA_0,xB_0,heat,Flow,n,dt) # Use one ODE solver to solve the whole system
    # println("These are the inputs for MPC_step_all")
    # println("T=",T_0)
    # println("Tin=",T0_in)
    # println("xB=",xB_0)
    # println("heat=",heat)
    # println("flow=",Flow)
    function odeodes!(du,u,p,t)
        for i=0:N-1 # N reactors in total
            du[3*i+1] = 1/V*(sum(n[j,i+1]*Flow[j,i+1]*u[3*(j-1)+1] for j=1:N) + n[i+1,N+1]*Flow[i+1,N+1]*T0_in[i+1] - sum(n[i+1,j]*Flow[i+1,j]*u[3*i+1] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*u[3*i+1]) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/u[3*i+1])*u[3*i+2])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/u[3*i+1])*u[3*i+3]) + heat[i+1]/rho/c_p/V  # Temperature of the i th reactor
            du[3*i+2] = 1/V*(sum(n[j,i+1]*Flow[j,i+1]*u[3*(j-1)+2] for j=1:N) + n[i+1,N+1]*Flow[i+1,N+1]*xA0 - sum(n[i+1,j]*Flow[i+1,j]*u[3*i+2] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*u[3*i+2]) + (-k1*exp(-E1/R_gas/u[3*i+1])*u[3*i+2])
            du[3*i+3] = 1/V*(sum(n[j,i+1]*Flow[j,i+1]*u[3*(j-1)+3] for j=1:N) - sum(n[i+1,j]*Flow[i+1,j]*u[3*i+3] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*u[3*i+3]) + k1*exp(-E1/R_gas/u[3*i+1])*u[3*i+2] + (-k2*exp(-E2/R_gas/u[3*i+1])*u[3*i+3])
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

    prob=ODEProblem(odeodes!,initial_vec,(0.0,dt))
    # prob=ODEProblem(odeodes!,initial_vec,(0.0,20*dt))
    soln=DifferentialEquations.solve(prob,Rosenbrock23())
    println("Next measurement is: ", last(soln.u))
    return last(soln.u)
    # a=soln.t
    # A=Array(soln)
    # return a,A[1,:],A[2,:],A[3,:]
    # return a,A[4,:],A[5,:],A[6,:]
    # return the time array and simulated T, xA, xB arrays

end

# Flows contains the info of all streams that connects
function findSS_all(T0_in,T_0,xB_0,n,Flow)
    # assume the connecting flows are constant
    # println("Here is the function of calculating steady state")
    # println(T0_in," and ",T_0," and ",xB_0," and ",n," and ",Flow)
    flow_start=zeros(N)
    heat_start=zeros(N)
    for i=0:N-1
        # flow_start[i+1] = -1/xA0*(sum(n[j,i+1]*Flow[j,i+1]*(1-xB_0[j]) for j=1:N) - sum(n[i+1,j]*Flow[i+1,j]*(1-xB_0[i+1])  for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*(1-xB_0[i+1]) + V*(-k1*exp(-E1/R_gas/T_0[i+1])*(1-xB_0[i+1]))) # assume all reactors must have an input stream
        # flow_start[i+1] = sum(n[i+1,j]*Flow[i+1,j] for j=1:N) + n[i+1,N+2]*Flow[i+1,N+2] - sum(n[j,i+1]*Flow[j,i+1]*(1-xB_0[j]) for j=1:N)
        flow_start[i+1] = sum(n[i+1,j]*Flow[i+1,j] for j=1:N) - sum(n[j,i+1]*Flow[j,i+1] for j=1:N) - 1/xB_0[i+1]*(-sum(n[j,i+1]*Flow[j,i+1]*xB_0[j] for j=1:N) + sum(n[i+1,j]*Flow[i+1,j]*xB_0[i+1]  for j=1:N) + V*(-k1*exp(-E1/R_gas/T_0[i+1])*(1-xB_0[i+1])+k2*exp(-E2/R_gas/T_0[i+1])*xB_0[i+1]))
        heat_start[i+1] = -rho*c_p*V*(1/V*(sum(n[j,i+1]*Flow[j,i+1]*T_0[j] for j=1:N) + n[i+1,N+1]*flow_start[i+1]*T0_in[i+1] - sum(n[i+1,j]*Flow[i+1,j]*T_0[i+1] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*T_0[i+1]) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_0[i+1])*(1-xB_0[i+1]))+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_0[i+1])*xB_0[i+1]))
        # P1=-rho*c_p*(sum(n[j,i+1]*Flow[j,i+1]*T_0[j] for j=1:N))
        P1=-rho*c_p*(n[i+1,N+1]*flow_start[i+1]*T0_in[i+1])
        # P2=-rho*c_p*( - sum(n[i+1,j]*Flow[i+1,j]*T_0[i+1] for j=1:N) )
        P2=-rho*c_p*(-n[i+1,N+2]*(sum(n[j,i+1]*Flow[j,i+1] for j=1:N)+n[i+1,N+1]*flow_start[i+1]-sum(n[i+1,j]*Flow[i+1,j] for j=1:N))*T_0[i+1])
        P3=-rho*c_p*(-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_0[i+1])*(1-xB_0[i+1]))+-rho*c_p*(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_0[i+1])*xB_0[i+1])
        # println("In ",i+1, " Reactor, the two part value are ", P1, " and ", P2, " and ", P3)
    end
    # println("flow_start=",flow_start)
    # println("heat_start=",heat_start)

    initial_vec2=zeros(3*N)

    for i=0:N-1
        initial_vec2[3*i+1]=flow_start[i+1]
        initial_vec2[3*i+2]=heat_start[i+1]
        initial_vec2[3*i+3]=1-xB_0[i+1]
    end

    function f!(du,u)
        for i=0:N-1 # N reactors in total
            du[3*i+1] = 1/V*(sum(n[j,i+1]*Flow[j,i+1]*T_0[j] for j=1:N) + n[i+1,N+1]*u[3*i+1]*T0_in[i+1] - sum(n[i+1,j]*Flow[i+1,j]*T_0[i+1] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*T_0[i+1]) + (-d_H1*mass/c_p*k1*exp(-E1/R_gas/T_0[i+1])*u[3*i+3])+(-d_H2*mass/c_p*k2*exp(-E2/R_gas/T_0[i+1])*xB_0[i+1]) + u[3*i+2]/rho/c_p/V
            du[3*i+2] = 1/V*(sum(n[j,i+1]*Flow[j,i+1]*u[3*(j-1)+3] for j=1:N) + n[i+1,N+1]*u[3*i+1]*xA0 - sum(n[i+1,j]*Flow[i+1,j]*u[3*i+3] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*u[3*i+3] + (-k1*exp(-E1/R_gas/u[3*i+1])*u[3*i+3]))
            du[3*i+3] = 1/V*(sum(n[j,i+1]*Flow[j,i+1]*xB_0[j] for j=1:N) - sum(n[i+1,j]*Flow[i+1,j]*xB_0[i+1] for j=1:N) - n[i+1,N+2]*Flow[i+1,N+2]*xB_0[i+1]) + k1*exp(-E1/R_gas/T_0[i+1])*(1-xB_0[i+1]) + (-k2*exp(-E2/R_gas/T_0[i+1])*xB_0[i+1])
        end
    end
    # u[3*i+1] is flow rate
    # u[3*i+2] is heating rate
    # u[3*i+3] is xA

    soln=nlsolve(f!,initial_vec2)
    xA_ss=zeros(N)
    heat_ss=zeros(N)
    flow_ss=zeros(N)
    for i=0:N-1
        xA_ss[i+1]=soln.zero[3*i+3]
        heat_ss[i+1]=soln.zero[3*i+2]
        flow_ss[i+1]=soln.zero[3*i+1]
    end
    return heat_ss,flow_ss
end

function permutate_weights(out_dir, disturbances)
    original_weights = [1,1e7,1e7,1e-5,1e7]
    powers_each_side = 2
    permutation_weights = zeros(Float64, (2*powers_each_side + 1,length(original_weights)))
    for i in 1:(2*powers_each_side + 1)
        for j in 1:length(original_weights)
            permutation_weights[i,j] = original_weights[j]*exp10(i - (powers_each_side + 1))
        end
    end
    display(permutation_weights)
    top_ten = fill(typemax(Float64), (10,10))

    num_permutations = 5^5 - 1 # iterate in base 5 through all possible permutations

    for i in ProgressBar(1735:1765)
    # for i in ProgressBar(0:num_permutations)
    # for i in 0:100
        base_five = string(i, base=5, pad=5)
        # println(base_five)
        current_weights = original_weights
        for j in 1:length(base_five)
            char = base_five[j]
            current_weights[j] = permutation_weights[parse(Int64, char) + 1, j]
        end
        q_T = current_weights[1]
        q_xA = current_weights[2]
        q_xB = current_weights[3]
        r_heat = current_weights[4]
        r_flow = current_weights[5]
        # println(current_weights)
        # discrepancies is an array of length 4 [qXb*dxB^2, qT*dT^2, r_flow*dFlow^2, r_heat*dHeat^2]
        discrepancies = MPC_tracking([0 0 1 1;0 0 1 1], disturbances,q_T,q_xA,q_xB,r_heat,r_flow,90,1000,[8 15];tmax=5000, print=false)
        # println("Discrepancies: $discrepancies")
        sum_discrepancies = sum(discrepancies)
        m = maximum(top_ten[:,10])
        if sum_discrepancies < m
            insert_index = findall(x -> x==m, top_ten[:,10])[1]
            top_ten[insert_index,1:5] = current_weights
            top_ten[insert_index,6:9] = discrepancies
            top_ten[insert_index,10] = sum_discrepancies
        end
        # for j in 1:10
        #     # println("$(top_ten[j,10]) $sum_discrepancies")
        #     is_full = maximum(top_ten[:,10]) != Inf
        #     # println("$(top_ten[:,10]) $is_full")
        #     if top_ten[j,10] > sum_discrepancies
        #         !is_full && top_ten[j,10] != Inf && continue
        #         top_ten[insert_index,1:5] = current_weights
        #         top_ten[insert_index,6:9] = discrepancies
        #         top_ten[insert_index,10] = sum_discrepancies
        #         break
        #     end
        # end
    end
    display(top_ten)
    println("writing top ten configurations to top_ten.txt")
    top_ten_file = out_dir * "\\top_ten.txt"
    touch(top_ten_file)
    file = open(top_ten_file, "w")
    writedlm(file, top_ten)
    close(file)

    originalDiscrepancy = MPC_tracking([0 0 1 1;0 0 1 1], disturbances,1,1e7,1e7,1e-3,1e9,90,1000,[8 15];tmax=5000) # no disturbance
    print("\nOriginal discrepancy: ")
    println(originalDiscrepancy[1] + originalDiscrepancy[4])
    return top_ten
end

function save_profile_images(inputMatrix, disturbances, out_dir)
    count = 1
    for row in eachrow(inputMatrix)
        println(row)
        q_T = row[1]
        q_xA = row[2]
        q_xB = row[3]
        r_heat = row[4]
        r_flow = row[5]
        image_name = join(row[1:5], "_")
        image_name = (out_dir * "\\Perm" * string(count) * "_" * image_name * ".png")
        MPC_tracking([0 0 1 1;0 0 1 1], disturbances,q_T,q_xA,q_xB,r_heat,r_flow,90,1000,[8 15];tmax=5000, print=false, save_plots=true, plot_name=image_name)
        count += 1
    end
end

out_dir = "C:\\Users\\sfay\\Documents\\Outputs"
disturbances = [6 6;0 0]
top_ten = permutate_weights(out_dir, disturbances)
# top_ten_hardcoded = [0.01	100000.0	1.0e11	0.0001	1.0e6	4.5902784576657645e-6	44.67795862520155	2.13733858592185e-6	7402.168692236633	4.5902784576657645e-6
#                     0.01	100000.0	1.0e11	0.001	1.0e9	0.005476772796985585	214532.18323354874	4.918302353195665e-6	349513.81350522436	0.005476772796985585
#                     0.01	100000.0	1.0e10	1.0000000000000002e-6	1.0e9	0.007322597410352353	109371.49364775658	3.485568289844138e-5	2.62339968185887e6	0.007322597410352353
#                     0.01	100000.0	1.0e11	1.0000000000000002e-6	1.0e8	0.007588939087214402	158621.9501849207	1.551582794482276e-5	1.0411854487997093e6	0.007588939087214402
#                     0.01	100000.0	1.0e11	1.0000000000000001e-7	1.0e9	0.010370141909242752	125577.67343643718	4.315955026410419e-5	3.006275735352627e6	0.010370141909242752
#                     0.01	100000.0	1.0e11	1.0e-5	1.0e7	0.010539706883296545	210056.02263620676	1.5477335065790005e-5	1.7222644535647202e6	0.010539706883296545
#                     0.01	100000.0	1.0e10	0.001	1.0e7	0.015604859673570717	145742.95127396716	1.950668784567637e-5	2.104809944045778e6	0.015604859673570717
#                     0.01	100000.0	1.0e10	0.001	1.0e8	0.016582382266799704	133381.66374726084	2.601153486813718e-5	2.0605073022839876e6	0.016582382266799704
#                     0.01	100000.0	1.0e11	0.0001	1.0e7	0.018010298125712493	261604.21391037246	1.1886056482999412e-5	943948.8612134649	0.018010298125712493
#                     0.01	100000.0	1.0e11	0.001	1.0e6	0.023132865203311474	283914.8446502505	1.2595620665598227e-5	967729.0801870388	0.023132865203311474]
out_dir = "C:\\Users\\sfay\\Documents\\Outputs\\Images\\"
save_profile_images(top_ten, disturbances, out_dir)

# MPC_tracking([0 0 1 1;0 0 1 1], [0 0;0 0],1,1e7,1e7,1e-3,1e9,90,1000,[8 15];tmax=5000) # no disturbance
# MPC_tracking([0 0 1 1;0 0  1 1], [10 10; 0 0],1,1e7,1e7,1e-3,1e9,90,1000,[8 15];tmax=5000) # disturbance on the first R
# MPC_tracking([0 0 1 1;0 0 1 1], [0 0;10 10],1,1e7,1e7,1e-3,1e9,90,1000,[8 15];tmax=5000) # disturbance on the second R
