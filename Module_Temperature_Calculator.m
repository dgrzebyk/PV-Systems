% NOCT conditions: 800 W/m^2, w = 1 m/s, Ta= 20 deg C, module mounted on a rack 
clear all
clc
%% Inputs

G_M = 750; % W/m^2 --> THIS CHANGES ACCORDING TO TIME, is an output of the solar calculator
Ta = 273.15 + 20; % Ambient temperature [K]
T_gr = 273.15 + 30; % temperature of the ground [K]
cloud_cover = 0; % [okta] (0 is clear sky and 8 is completely cloudy sky)
w = 4; % wind speed [m/s]

eta = 0.18; % module efficiency at STC
T_INOCT = 273.15 + 54; % [K] - is this module temperature under INOCT?
theta_M = 30; % tilt angle of a module! [degrees]
L = 1.5; % length of the rectangle cross section
W = 1; % width of the rectangle cross section
reflectivity = 0.1; % reflectivity of the module
emis_top = 0.84;
emis_back = 0.89;
k = 0.026; % heat conductivity of air [W*m^-1*K^-1]
% inst_conf = 0; % ? rack mount / direct mount / standoff
sigma = 5.6704 * 10 ^ (-8); % ?
v = 1.837*10^(-5); % viscosity of air [kg * m^-1 * s^-1]
T_M = 293.15; % [K] Just for first iteration! Temperature of the module

%% Calculations

if cloud_cover > 6
    T_sky = T_ambient;
else
    T_sky = 0.0552 * Ta ^ (3/2);
end

%ym = 1; % module height - given
%yr = 1; % anemometer height - given
%wr = 10; % wind speed registered by anemometer - given
%w = wr * (ym / yr)^0.2; % wind speed adjusted to solar module height
% Create a loop to show the dependence of module yield on its height or
% wind speed

alpha = (1 - reflectivity)*(1 - eta); % absorptivity of a solar module

for i = 1:3
    
    hr_sky = emis_top * sigma * (T_M^2 + T_sky^2) * (T_M + T_sky);
    hr_gr = emis_back * sigma * (T_M^2 + T_gr^2) * (T_M + T_gr);
    Dh = (2*L*W) / (L + W);
    Gr = ((9.8 * (1/Ta) * (T_M - Ta) * Dh^3) / v^2) * sind(theta_M); % Grashof number
    Pr = 0.708; % Prandtl number for air [-]
    Nu = 0.21 * (Gr * Pr)^0.32; % Nusselt number - expresses the ratio between conductive and convective heat transfer
    Re = (w*Dh)/v; % Reynolds number
    %h_lam_forced = w^0.5;
    %h_turb_forced = w^0.8;

    h_lam_forced = ((0.86 * Re^(-0.5)) / (Pr^0.67)) * 1.204 * 1005 * w;
    h_turb_forced = ((0.028 * Re^(-0.2)) / (Pr^0.4)) * 1.204 * 1005 * w;
    
    if w <= 3
        h_forced = h_lam_forced;
    else
        h_forced = h_turb_forced;
    end
    
    h_free = Nu * k / Dh;
    hc_top = (h_forced ^ 3 + h_free ^ 3)^(1/3); % it is called also h_mixed

    R = (alpha * G_M - hc_top * (T_INOCT - Ta) - emis_top * sigma * (T_INOCT^4 - T_sky^4) ) / (hc_top * (T_INOCT - Ta) + emis_top * sigma * (T_INOCT^4 - T_sky^4) );
    
    hc_bot = R * hc_top; % eq. G.17

    hc = hc_top + hc_bot;

    T_M = (alpha * G_M + hc * Ta + hr_sky * T_sky + hr_gr * T_gr) / (hc + hr_sky + hr_gr);
    
    R2 = (hc_bot * (T_INOCT - Ta) + emis_back * sigma * (T_INOCT^4  - T_gr^4)) / (hc_top * (T_INOCT - Ta) + emis_back * sigma * (T_INOCT^4 - Ta^4));
    
    % Outputs for the user
    fprintf('Iteration: %d \n',i);
    fprintf('h_forced = %5.4f \n',h_forced);
    fprintf('h_free = %5.4f \n',h_free);
    fprintf('hc_T = %5.4f \n',hc_top);
    fprintf('hc_B = %5.4f \n',hc_bot);
    fprintf('hr_sky = %5.4f \n',hr_sky);
    fprintf('hr_gr = %5.4f \n',hr_gr);
    fprintf('T_M = %5.4f \n',T_M);
    %fprintf('Module temperature is: %5.1f Kelwins \n',T_M);
    fprintf(' \n');
    
end

%% ISSUES
% h_forced is wrong
% R shouldn't be negative!!!