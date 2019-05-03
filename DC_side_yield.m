% Creating a script that calculates the overall performance of a PV module cascaded with the temperature calculator

% General data
n = 4400; % number of panels
e_top = 0.84; % top glass emissivity
e_bottom = 0.893; % bottom surface emissivity
R = 0.1; % Reflectivity of the module
T_M = 20; % Initial module temperature - for loop to start
mod_pos = 1; % module position [m]

% Environmental conditions
G_M = 700; % irradiance [W/m^2]
w = 5; % wind speed [m/s]
Ta = 15; % Ambient temperature
cloud_cover = 4; % [okta]
T_gr = 13 + 273.15;
Pr = 0.71;
p_atm = 1; % [bar]

% Module data
eta_STC = 0.198; % module efficiency
A = 0.798 * 1.58; % m^2
Voc = 53.2;
Isc = 6.03;
Voc_mpp = 44.3;
Isc_mpp = 5.65;
P_mpp = 250; % Wp
k_v = -0.235;
k_i = 0.055;
k_P = -0.258;
k_eta = 0; % what's the correct value here?
T_NOCT = 44 + 273.15; % Kelwin
kappa = -0.0035; % /deg C --> for c-Si modules (p.331)

%% Calculations
% assignment 3 code should calculate module temperature at stc and under given irradiance
Voc_temp = Voc_mpp + k_v * (Tm - T_STC);
Isc_temp = Isc_mpp + k_i * (Tm - T_STC);
Pmpp_temp = P_mpp + k_P * (Tm - T_STC);
eta_temp = eta_STC + k_eta * (Tm - T_STC);