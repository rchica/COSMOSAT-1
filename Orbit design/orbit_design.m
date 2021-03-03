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
tau = (3600*24)*(365.242199/(1+365.242199));    %Sidereal day
Earth_Omega = (pi/180)*(360)/tau;               %Mean sidereal motion of the Earth

%Sun requirements 
delta = deg2rad(2);                             %Solar declination
eps = 0;                                        %Equation of time

%Orbit requirements 
rep_days = 1;                                   %Days between identical groundtracks
rev_day = 14;                                   %Needed revolutions per day
LTANd = 12;                                     %Mean Local Time of the Ascending Node

%% Altitude selection
P = 86400/rev_day;                              %Nodal period
a_d = (mu*(P/2*pi)^2)^(1/3);                    %Desired orbital altitude 

%% Inclination selection
%General inclination evaluation
dOmega = Earth_Omega;                           %Orbital precession rate
dh = 100;                                       %Altitude step
a_max = 1000e3;                                 %Orbital altitude (circular orbit) over the geoid
a = a_e:dh:a_e+a_max;                           %Orbital altitude
n = sqrt(mu./a.^3);                             %Mean orbital motion
p = a;                                          %Semilatus rectum of a circular orbit
i = acos(dOmega./(-(2/3)*J2*n.*(a_e./p).^2));   %Inclination relationship with the orbital altitude

%Specific orbit determination
n_d = sqrt(mu/a_d^3);                           %Desired mean orbital motio
p_d = a_d;                                      %Desired semilatus rectum of a circular orbit
i_d = acos(dOmega/(-(2/3)*J2*n_d*(a_e/p_d)^2)); %Desired inclination

%% Fundamental interval 
dL = 360*(rep_days/rev_day);                    %Earth angle between adcent groundtracks

%% Solar angle 
s = [cos(delta)*cos(eps); cos(delta)*cos(eps); sin(delta)];     %Sun vector 
h = [sin(i)*sin(MOmega); -sin(i)*cos(MOmega); cos(i)];          %Angular momentum vector 
beta = asin(dot(h, s));                                         %Solar angle 

%% Shadow time 
eta = asin(a_e/a_d);                        
nu = 2*acos(cos(eta)/cos(beta));            
dTimeShadow = nu/360*P;                                         %Spent time in shadow

%% Results
fprintf("LTAN: %.f \n", LTANd);
fprintf("Orbital altitude: %.f km \n", a_d/10^3);
fprintf("Orbital inclination: %.f \n", i_d);
fprintf("RAAN: %.f \n", RAAN);
fprintf("Time in shadow: %.f h\n", dTimeShadow/3600);

%%
figure(1) 
hold on
plot(a/1000, rad2deg(i), 'b');
plot(a_d/1000, rad2deg(i_d), 'or'); 
hold off 
grid on
xlabel('Orbital altitude over the geoid (km)'); 
ylabel('Orbital inclination (deg)'); 
legend('Sun-synchronous orbit', 'Design point for COSMOSAT-1');
title('Orbit design point for COSMOSAT-1');
