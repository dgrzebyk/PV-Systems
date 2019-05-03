% We assume 1m^2 PV module, 60 deg tilt, orientation: South

tic

phi = 52.0116;%input('Insert latitude [deg]. Positive for N, negative for S \n');
lambda = 4.3571;%input('Insert longitude [deg]. Positive for E, negative for W \n');
t = [2017,1,1,0,0,0];%input('Insert date as a vector [Y,M,D,H,M,S] \n');
time_zone = 2;%input('Please insert your UTC time zone: \n');
date_1 = datetime(t) + hours(time_zone)*(-1);
date_2 = datetime(2000,1,1,12,0,0);

% Creation of empty variables
As = zeros(8760,1);
as = zeros(8760,1);
Ie_dir = zeros(8760,1);
Ie_global = zeros(8760,1);
G_direct = zeros(8760,1);

% Determining input constants
Ie0 = 1361; % W/m^2
c = 0.14;
h = 0; % input('Please insert the altitude of the observer: \n');
% Am = 180; % azimuth of the module where 0 is North

for Am = 0:15:360
    fprintf('Am = %d \n',Am);
    for theta_M = 0:10:90 % module tilt angle

        am = 90 - theta_M; % source: slide 40

        for i = 1:8760

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

            % Local mean sidereal time
            theta_L = GMST * 15 + lambda;

            v_s = -sind(theta_L)*cosd(lambda_s)+cosd(theta_L)*cosd(epsilon)*sind(lambda_s);
            funny_sign = -sind(phi)*cosd(theta_L)*cosd(lambda_s)-(sind(phi)*sind(theta_L)*cosd(epsilon)-cosd(phi)*sind(epsilon))*sind(lambda_s);

            if funny_sign > 0 && v_s > 0
                As(i) = atand(v_s / funny_sign);
            elseif funny_sign < 0
                As(i) = atand(v_s / funny_sign) + 180;
            else
                As(i) = atand(v_s / funny_sign) + 360;
            end

            sin_as = cosd(phi)*cosd(theta_L)*cosd(lambda_s)+(cosd(phi)*sind(theta_L)*cosd(epsilon)+sind(phi)*sind(epsilon))*sind(lambda_s);
            as(i) = asind(sin_as);

            % When the Sun is below horizon we set solar altitude to be zero
            if as(i) < 0
                as(i) = 0;
            end

            %theta_M(i) = 90 - as(i);
            AM =( 1 / ( sind(as(i)) + 0.50572 * (6.07995 + as(i))^(-1.6364) ));
            Ie_dir = Ie0 * ( (1-c*h) * 0.7 ^ AM ^ 0.678 + c * h ); % This is DNI
%            Ie_global(i) = 1.1 * Ie_dir;

            AOI = acosd( cosd(am)*cosd(as(i))*cosd(Am-As(i)) + sind(am)*sind(as(i)) );

            % Direct Irradiance calculation
            G_direct(i) = Ie_dir * cosd(AOI);

            if G_direct(i) < 0
                G_direct(i) = 0;
            end

            date_1 = date_1 + hours(1);
        end

        tot_irr = sum(G_direct); % annual energy yield over a year kWh
        results(Am/15+1,theta_M/10+1) = tot_irr*0.001;
    end
    date_1 = date_1 - hours(8760);
end
%[x,y] = max(results);

%fprintf('Maximal irradiance is incident for tilt angle of: %d degrees \n',y);

%var = results(1:25,:)';
contourf(0:15:360,0:10:90,results)

toc