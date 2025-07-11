# https://lanl-ansi.github.io/PowerModelsONM.jl/stable/tutorials/Use%20Case%20Examples.html
# Use case example for ieee13 

begin 
    import DataFrames as DF
    import CSV 
    import PowerModelsONM as ONM
end
onm_path = "G:\\My Drive\\Research\\ModularUnit"
ieee13_data = ONM.prepare_data!(
    Dict{String, Any}(
        "network" => joinpath(onm_path, "network.ieee13mod.dss"),
		"settings" => joinpath(onm_path, "ieee13_settings.json"),
		"events" => joinpath(onm_path, "ieee13_events.json"),
    )
)

ieee13_mn = deepcopy(ieee13_data["network"])

ieee13_mip_solver = ONM.build_solver_instances(;
solver_options=ieee13_data["settings"]["solvers"])["mip_solver"]

#=======
Use-Case I: Radiality constraint comparison
In the following use-case, we explore the effect of the radiality constraint on the number of loads shed throughout the time-series.

The radiality constraint implemented by default in PowerModelsONM is based on a multi-commodity flow formulation as defined in S. Lei, C. Chen, Y. Song, and Y. Hou, “Radiality constraints for resilient reconfiguration of distribution systems: Formulation and application to microgrid formation,” IEEE Transactions on Smart Grid, vol. 11, no. 5, pp. 3944–3956, 2020.

First we obtain the results for the case where radiality is enforced, which is the default state 
=======#
result_rad_enabled = ONM.optimize_switches(
	ieee13_mn,
	ieee13_mip_solver;
	formulation=ONM.LPUBFDiagPowerModel,
	algorithm="rolling-horizon",
	problem="block"
)

#====
Next, to obtain the results for the case where radiality is not enforced, we need to set the option 'options/constraints/disable-radiality-constraint' to false, which we can do with the set_setting! helper function which will return the multinetwork data structure.
====#
ieee13_mn_rad_disabled = ONM.set_setting!(
	deepcopy(ieee13_data),
	("options","constraints","disable-radiality-constraint"),
	true
)

result_rad_disabled = ONM.optimize_switches(
	ieee13_mn_rad_disabled,
	ieee13_mip_solver;
	formulation=ONM.LPUBFDiagPowerModel,
	algorithm="rolling-horizon",
	problem="block"
)

