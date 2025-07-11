using LinearAlgebra, Plots, NLsolve

# Parameters:
R_gas = 8.314 # Gas constant [J / (mol K)]
Vr = 2.675E-6 # Reactor volume [m3]
N = 200
dt = 1 # seconds
W_cat = 5E-3 # Catalyst Weight [kg]

# Vessel parameters
P_reactor = 1E5 # Reactor Pressure [Pa]
T = 1000 # Reactor (Catalyst) Temperature [K]
C_reactor = P_reactor / (R_gas * T) # Charged concentration, Ideal gas [mol / m3]

# Debugging - all concentrations start at 3%, except feed component, which are 1/3 the inlet flow composition [mol/m3]
C_CO_init = C_reactor * 0.033
C_H2_init = C_reactor * 0.033
C_CH4_init = C_reactor * 0.225
C_H2O_init = C_reactor * 0.675
C_CO2_init = C_reactor * 0.033

# Feed flow conditions
F_CH4 = 25 * (7.45 * 10^-7) * (1 / 300) # [mol / s]
F_H2O = 75 * (7.45 * 10^-7) * (1 / 300) # [mol / s]
F_H2 = 0 * (7.45 * 10^-7) * (1 / 55) # [mol / s]
F_T0 = sum([F_CH4, F_H2O, F_H2]) # [mol / s]

C_CO_vec = zeros(N + 1)
C_H2_vec = zeros(N + 1)
C_CH4_vec = zeros(N + 1)
C_H2O_vec = zeros(N + 1)
C_CO2_vec = zeros(N + 1)

C_CO_vec[1] = C_CO_init
C_H2_vec[1] = C_H2_init
C_CH4_vec[1] = C_CH4_init
C_H2O_vec[1] = C_H2O_init
C_CO2_vec[1] = C_CO2_init


# Equilibrium Constant for SMR [bar2] (Abbas 2017)
K_I(T) = exp((-26830 / T) + 30.114)

# Equilibrium Constant for WGS [unitless] (Abbas 2017)
K_II(T) = exp((4400 / T) - 4.036)

# Activation energy and preexponential factor for SMR (Abbas 2017)
E1 = 257.01 * 1E3 # [J/mol]
ko1 = 5.19E9 * 1E3 # mol bar^(0.5) / (kg_cat s)

# Activation energy and preexponential factor for WGS (Abbas 2017)
E2 = 236.70 * 1E3 # [J/mol]
ko2 = 1.32E10 * 1E3 # mol bar^(0.5) / (kg_cat s)

# Kinetic rate constant for SMR
k_1(T) = ko1 * exp(-E1 / (R_gas * T))

# Kinetic rate constant for WGS
k_2(T) = ko2 * exp(-E2 / (R_gas * T))

# Heat of adsorption for each species [J / mol] (Xu 1989)
ΔH_CO = -70.65 * 1E3
ΔH_H2 = -82.90 * 1E3
ΔH_CH4 = -38.28 * 1E3
ΔH_H2O = 88.68 * 1E3

# Reference adsorption constant for each species (Xu 1989)
K_o_CO = 8.32E-5 # [bar]
K_o_H2 = 6.12E-9 # [bar]
K_o_CH4 = 6.65E-4 # [bar]
K_o_H2O = 1.77E5 # [unitless]

# Adsorption constant for each species (Xu 1989)
K_CO(T) = K_o_CO * exp(ΔH_CO / (R_gas * T)) # [bar]
K_H2(T) = K_o_H2 * exp(ΔH_H2 / (R_gas * T)) # [bar]
K_CH4(T) = K_o_CH4 * exp(ΔH_CH4 / (R_gas * T)) # [bar]
K_H2O(T) = K_o_H2O * exp(ΔH_H2O / (R_gas * T)) # [unitless]

# Definition of partial pressure [bar]
p(C, T) = (C * R_gas * T) / (1E5)

# Unit less term used in reaction kinetics (Abbas 2017)
Ω(T, C_CO, C_H2, C_CH4, C_H2O) = 1 + K_CO(T) * p(C_CO, T) + K_H2(T) * p(C_H2, T) + K_CH4(T) * p(C_CH4, T) + K_H2O(T) * (p(C_H2O, T) / p(C_H2, T))

# Rates of reaction SMR, WGS [mol / (kg_cat s)]
R_SMR(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) = ((k_1(T) / ((p(C_H2, T))^(2.5))) * (p(C_CH4, T) * p(C_H2O, T) - (((p(C_H2, T)^3) * p(C_CO, T)) / K_I(T))) * (1 / (Ω(T, C_CO, C_H2, C_CH4, C_H2O)^2)))
R_WGS(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) = ((k_2(T) / ((p(C_H2, T)))) * (p(C_CO, T) * p(C_H2O, T) - (((p(C_H2, T)) * p(C_CO2, T)) / K_II(T))) * (1 / (Ω(T, C_CO, C_H2, C_CH4, C_H2O)^2)))

# Outlet volumetric flow rate [m3 / s]
q(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) = (F_T0 + 2 * (R_SMR(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2)) * W_cat) / (P_reactor / (R_gas * T))

R_SMR_vec = zeros(N)
R_WGS_vec = zeros(N)

for i = 1:N
    C_CO = C_CO_vec[i]
    C_H2 = C_H2_vec[i]
    C_CH4 = C_CH4_vec[i]
    C_H2O = C_H2O_vec[i]
    C_CO2 = C_CO2_vec[i]

    # Calculate rates
    R_SMR_rate = R_SMR(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2)
    R_WGS_rate = R_WGS(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2)

    # Concentration updates
    C_CO_vec[i+1] = C_CO + ((1 / Vr) * ((R_SMR_rate - R_WGS_rate) * W_cat - q(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) * C_CO)) * dt
    C_H2_vec[i+1] = C_H2 + ((1 / Vr) * (F_H2 + (3 * R_SMR_rate + R_WGS_rate) * W_cat - q(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) * C_H2)) * dt
    C_CH4_vec[i+1] = C_CH4 + ((1 / Vr) * (F_CH4 - R_SMR_rate * W_cat - q(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) * C_CH4)) * dt
    C_H2O_vec[i+1] = C_H2O + ((1 / Vr) * (F_H2O - (R_SMR_rate + R_WGS_rate) * W_cat - q(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) * C_H2O)) * dt
    C_CO2_vec[i+1] = C_CO2 + ((1 / Vr) * (R_WGS_rate * W_cat - q(T, C_CO, C_H2, C_CH4, C_H2O, C_CO2) * C_CO2)) * dt

    R_SMR_vec[i] = R_SMR_rate
    R_WGS_vec[i] = R_WGS_rate

end


plot(C_CO_vec, label="CO", lw=4)
plot!(C_CH4_vec, label="CH4", lw=4)
plot!(C_H2O_vec, label="H2O", lw=4)
plot!(C_CO2_vec, label="CO2", lw=4)
plot!(C_H2_vec, label="H2", lw=4)


