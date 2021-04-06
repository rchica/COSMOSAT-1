%% COSMOSAT-1 ADCS Team %% 
% 23/02/21

%% Orbit design script
% This script provides an interface to trade-off orbital elements 
% for the COSMOSAT-1 mission. 

% Orbits in consideration are sun-synchronous from 300 km to 800 km (budget
% constraints).

% All units are in S.I.

% Github: https://github.com/cosmos-urjc/COSMOSAT-1.git

%% Initial data and orbit requirements 
%Earth data
mu = 3.986e14;                                  %Earth gravitational parameter
J2 = 1.08263e-3;                                %Second zonal harmonic of the Earth
a_e = 6378.14e3;                                %Mean Earth radius
SolarYear = 365.242199;                         %Solar year
tau = (3600*24)*(SolarYear/(1+SolarYear));      %Sidereal day
Earth_Omega = (2*pi/SolarYear)/(3600*24);       %Earth Mean motion
Sun_Omega = (360)/SolarYear;                    %Earth rate relative to vernal equinox
 
%Time relationships 
vernalTime = datetime(2022,3,20,0,37,0,0);      %Day of the vernal equinox
launchTime = datetime(2022,7,7,12,0,0,0);       %Launch day
vernalTime = datenum(vernalTime);               %Day of the vernal equinox
time = datenum(launchTime);                     %Launch day

%Sun requirements (this must change as a function of the launch date. See Sun's analemma)
delta = deg2rad(2);                             %Solar declination
eps = 0;                                        %Equation of time

%Orbit requirements 
rep_days = 1;                                   %Days between identical groundtracks
rev_day = 15;                                   %Needed revolutions per day
LTANd = 12;                                     %Mean Local Time of the Ascending Node

%% Altitude selection
P = 86400/rev_day;                              %Nodal period
a_d = (sqrt(mu)*P/(2*pi))^(2/3);                %Desired orbital altitude 

%% Inclination selection
%General inclination evaluation
dOmega = Earth_Omega;                           %Orbital precession rate
dh = 100;                                       %Altitude step
a_max = 1000e3;                                 %Orbital altitude (circular orbit) over the geoid
a = a_e:dh:a_e+a_max;                           %Orbital altitude
e = 0;                                          %Orbital eccentricity
omega = pi/2;                                   %Argument of perigee for frozen orbits

%Inclination relationship with the orbital altitude
n = sqrt(mu./a.^3);                             %Orbital mean motion
p = a*(1-e^2);                                  %Semilatus rectum of the orbit
i = acos((-2*dOmega*p.^2)./(3*J2*a_e^2.*n));    %Inclination 

%Desired inclination
n_d = sqrt(mu/a_d^3);                           %Orbital mean motion
p_d = a_d*(1-e^2);                              %Semilatus rectum of the orbit
i_d = acos((-2*dOmega*p_d^2)/(3*J2*a_e^2*n_d)); %Desired inclination

%Perturbed desired inclination
tol = 1e-10;
error = 1;
while (error >= tol)
    %J2 perturbed mean motion
    n_p  = n_d*(1+(3/2)*J2*(a_e/a_d)^2*sqrt(1-e^2)*(1-(3/2)*sin(i_d)^2));
    i_p = acos((-2*dOmega*p_d^2)/(3*J2*a_e^2*n_p));
    error = abs(i_p-i_d);
    i_d = i_p;
end

%% RAAN selection
RAAN_d = (360/24)*LTANd-180+(Sun_Omega)*(time-vernalTime);
RAAN_d = mod(RAAN_d, 360);

%% Fundamental interval 
dL = 360*(rep_days/rev_day);                    %Earth angle between adcent groundtracks

%% Solar angle 
s = [cos(delta)*cos(eps); cos(delta)*sin(eps); sin(delta)];     %Sun vector 
h = [sin(i_d)*sin(RAAN_d); -sin(i_d)*cos(RAAN_d); cos(i_d)];    %Angular momentum vector 
beta = asin(dot(h,s));                                          %Solar angle 

%% Shadow time 
eta = asin(a_e/a_d);                        
nu = 2*acos(cos(eta)/cos(beta));            
dTimeShadow = rad2deg(nu)/360*P;                                         %Spent time in shadow

%% Results
fprintf("LTAN: %.4f h \n", LTANd);
fprintf("Orbital altitude: %.4f km \n", (a_d-a_e)/10^3);
fprintf("Orbital inclination: %.8f deg \n", rad2deg(i_d));
fprintf("RAAN: %.8f deg \n", RAAN_d);
fprintf("Time in shadow: %.4f min \n", dTimeShadow/60);

figure(1) 
hold on
plot((a-a_e)/1000, rad2deg(i), 'b');
plot((a_d-a_e)/1000, rad2deg(i_d), 'or'); 
hold off 
grid on
xlabel('Orbital altitude over the geoid (km)'); 
ylabel('Orbital inclination (deg)'); 
legend('Sun-synchronous orbit', 'Design point for COSMOSAT-1');
title('Orbit design point for COSMOSAT-1');
