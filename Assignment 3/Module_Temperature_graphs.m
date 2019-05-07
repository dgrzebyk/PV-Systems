% NOCT conditions: 800 W/m^2, w = 1 m/s, Ta= 20 deg C, module mounted on a rack
clear all
clc
%% Inputs

G_M = 750; % W/m^2 --> THIS CHANGES ACCORDING TO TIME, is an output of the solar calculator
Ta = 273.15 + 20; % Ambient temperature [K]
T_gr = 273.15 + 30; % temperature of the ground [K]
cloud_cover = 0; % [okta] (0 is clear sky and 8 is completely cloudy sky)

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
sigma = 5.6704 * 10 ^ (-8); % ?ARGHH!!! 5.6704*10^(-8)
d_v = 1.837*10^(-5); % %kg/(m*s) dynamic viscocity of air
v = d_v/1.204; %kinematic viscocity
R_noct = 0.2929;

% First iteration variables
Tm = 293.15; % [K] Just for first iteration! Temperature of the module
R = R_noct; % just for the first iteration

%% Calculations

if cloud_cover > 6
    T_sky = T_ambient;
else
    T_sky = 0.0552 * Ta ^ (3/2);
end

alpha = (1 - reflectivity)*(1 - eta); % absorptivity of a solar module

for j = 1:5
    for w = 1:10
        for i = 1:3
            
            G_M = 200 * j;
            hr_sky = emis_top*sigma*(Tm^2+T_sky^2)*(Tm+T_sky);
            hr_gr = emis_back * sigma * (Tm^2 + T_gr^2) * (Tm + T_gr);
            Dh = (2*L*W) / (L + W);
            Gr = ((9.8 * (1/Ta) * (Tm - Ta) * Dh^3) / v^2) * sind(theta_M); % Grashof number
            Pr = 0.708; % Prandtl number for air [-]
            Nu = 0.21 * (Gr * Pr)^0.32; % Nusselt number - expresses the ratio between conductive and convective heat transfer
            Re = (w*Dh)/v; % Reynolds number

            h_lam_forced = ((0.86 * Re^(-0.5)) / (Pr^0.67)) * 1.204 * 1005 * w;
            h_turb_forced = ((0.028 * Re^(-0.2)) / (Pr^0.4)) * 1.204 * 1005 * w;

            if w <= 3
                h_forced = h_lam_forced;
            else
                h_forced = h_turb_forced;
            end

            h_free = Nu * k / Dh;
            hc_top = (h_forced ^ 3 + h_free ^ 3)^(1/3); % it is called also h_mixed

            % WHY DO WE USE R_NOCT AT ALL ITERATIONS?
            % R = (alpha*G_M-hc_top*(T_INOCT-Ta)-emis_top*sigma*(T_INOCT^4-T_sky^4)) / (hc_top*(T_INOCT-Ta)+emis_top*sigma*(T_INOCT^4-T_sky^4) ); %G.16

            hc_bot = R_noct * hc_top; % eq. G.17

            hc = hc_top + hc_bot;

            Tm = (alpha * G_M + hc * Ta + hr_sky * T_sky + hr_gr * T_gr) / (hc + hr_sky + hr_gr);
            Tm_r(w,i,j) = Tm; % add G_M
        end
    end
end

    
scatter(1:10,Tm_r(:,i,1)-273.15,25,[0 0.4470 0.7410],'filled')
hold on
scatter(1:10,Tm_r(:,i,2)-273.15,25,[0.8500 0.3250 0.0980],'filled')
scatter(1:10,Tm_r(:,i,3)-273.15,25,[0.9290 0.6940 0.1250],'filled')
scatter(1:10,Tm_r(:,i,4)-273.15,25,[0.4940 0.1840 0.5560],'filled')
scatter(1:10,Tm_r(:,i,5)-273.15,25,[0.4660 0.6740 0.1880],'filled')
legend('Gm = 200','Gm = 400','Gm = 600','Gm = 800','Gm = 1000')
xlabel('wind speed [m/s]')
ylabel('Module temperature [°C]')
title('Module temperature vs wind speed')

%% ISSUES
% kinematic viscosity of air!!!
% line 26 - R_noct issue. Why do we use it?