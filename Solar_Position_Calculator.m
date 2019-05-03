clear all
clc

phi = input('Insert latitude [deg]. Positive for N, negative for S \n');
lambda = input('Insert longitude [deg]. Positive for E, negative for W \n');
t = input('Insert date as a vector [Y,M,D,H,M,S] \n');
time_zone = input('Please insert your UTC time zone: \n');
date_1 = datetime(t) + hours(time_zone)*(-1);
date_2 = datetime(2000,1,1,12,0,0);

% Time should be automatically adjusted to UTC - to be added

time_difference = date_1 - date_2; % hours,minutes,seconds
D = hours(time_difference)/24;
% Mean longitude
q = 280.459 + 0.98564736 * D;
% Normalization
while q > 360
    q = q - 360;
end
% Mean anomaly of the Sun
g = 357.529 + 0.98560028 * D;
% Normalization
while g > 360
    g = g - 360;
end
% Ecliptic longitude of the Sun
lambda_s = q + 1.915 * sind(g) + 0.020*sind(2*g);
% First rotation
epsilon = 23.429 - 0.00000036 * D;
T = D/36525;
GMST = 18.697374558 + 24.06570982441908 * D + 0.000026 * T^2;

% GMST Normalization
counter = 0;
while GMST > 24
    GMST = GMST - 24;
    counter = counter + 1;
end
% Longitude of the observer
lambda_0 = 1;
% Local mean sidereal time
theta_L = GMST * 15 + lambda;

v_s = -sind(theta_L)*cosd(lambda_s)+cosd(theta_L)*cosd(epsilon)*sind(lambda_s);
funny_sign = -sind(phi)*cosd(theta_L)*cosd(lambda_s)-(sind(phi)*sind(theta_L)*cosd(epsilon)-cosd(phi)*sind(epsilon))*sind(lambda_s);

if funny_sign > 0 && v_s > 0
    As = atand(v_s / funny_sign);
elseif funny_sign < 0
    As = atand(v_s / funny_sign) + 180;
else
    As = atand(v_s / funny_sign) + 360;
end

sin_as = cosd(phi)*cosd(theta_L)*cosd(lambda_s)+(cosd(phi)*sind(theta_L)*cosd(epsilon)+sind(phi)*sind(epsilon))*sind(lambda_s);
as = asind(sin_as);

fprintf('Altitude is: %5.2f degrees \n', as)
fprintf('Azimuth is: %5.2f degrees \n',As)