% Creating a script that calculates the overall performance of a PV module cascaded with the temperature calculator
clear all
clc

tic
% run('Module_Temperature_Calculator.m');

% Physical constants
kB = 1.38*10^(-23);
q = 1.602*10^(-19);
Pr = 0.71;

% General data
nm = 4400; % number of modules
emis_top = 0.84; % top glass emissivity
emis_bot = 0.893; % bottom surface emissivity
R = 0.1; % Reflectivity of the module
mod_pos = 1; % module position [m]
% assignment 3 code should calculate module temperature at stc and under given irradiance
Tstc = 273.15 + 25; % temperature under Standard Test Conditions [K]
Tm = 312.91; % Temperature of the module under the given irradiance [K]
% HERE THE TEMP. CALCULATOR SHOULD BE INCLUDED (when it is correct)

% Environmental conditions
Gm = 1000; % irradiance [W/m^2]
w = 5; % wind speed [m/s]
Ta = 273.15 + 15; % Ambient temperature
cloud_cover = 4; % [okta]
T_gr = 273.15 + 13; % temperature of the ground
% p_atm = 1; % [bar] What units should I use?

% Module data
eta_stc = 0.194; % module efficiency
Am = 0.798 * 1.58; % m^2
Voc = 53;
Isc = 5.86;
Voc_mpp = 44.3;
Isc_mpp = 5.54;
Pnoct = 187.3; % W
P_mpp = 245; % Wp
k_v = -0.235; % %/C
k_i = 0.055; % %/C
k_P = -0.258; % %/C
% k_eta = 0; % what's the correct value here?
T_NOCT = 273.15 + 44; % Kelwin
kappa = -0.0035; % /deg C --> for c-Si modules (p.331)
n = 1.5; % ideality factor [-]

%% Calculations
% Do I need the below variables? Can I just use eq. 20.20 to avoid using equations 20.21 - 20.24?
% These are slide 14 calculations - G = const., Tm = variable
% Voc_temp = Voc_mpp + k_v * (Tm - Tstc);
% Isc_temp = Isc_mpp + k_i * (Tm - Tstc);
% Pmpp_temp = P_mpp + k_P * (Tm - Tstc);
% eta_temp = Pmpp_temp / (Am * 1000);

% Write equation 20.8 for reference!

% Now I want both irradiance and temperature to vary

FF = (Voc_mpp * Isc_mpp) / (Voc * Isc);
Voc_25C = Voc + ((n*kB*293.15)/(q)) * log(Gm/1000);
Isc_25C = Isc * Gm / 1000;
Pmpp_25C = FF * Voc_25C * Isc_25C;
eta_25C = Pmpp_25C / (Am * Gm) ; % eq. 20.24
eta = eta_25C * (1 + kappa*(Tm - 298.15)); % eta(Tm,Gm) eq. 20.25

Pstc = eta_stc * 1000 * Am;
Pdc = eta * Gm * Am; % this eta is temperature dependent

P_inst = nm * P_mpp;
P_real_inst = nm * Pdc;

fprintf('The rated power of the Blackfriars Bridge Power Plant is %d kW, but its real power is %5.2f kW \n', P_inst/1000, P_real_inst/1000);

toc

%% ISSUES
% Q: How can I determine the ideality factor of my PV module?
% A: We assume it to be 1.5 - vide p.331

% Q: Why does the equation 20.16 give me different results than the module
% specification? Equation gives eta=0.28, spec. eta = 0.198
% A: That's because in eq. 20.16 temperature is not included.

% Is the saturation current calculated correctly? Did I use the correct temperature?
% Tm - should be imported from the previous script, for now it is set for STC
% Tstc - is it module temperature or ambient temperature?
% slide 16 (L4) - do I need those equations?

% Simulation with wind speed vs energy yield is required.