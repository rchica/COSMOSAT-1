%% CosmosSat-1 ADCS Team %% 
% 23/02/21

%% Preliminary disturbance study
% This script provides an interface to study attitude and orbital
% disturbances for the CosmoSat-1 mission. 

% Orbits in consideration are sun-synchronous from 300 km to 800 km.
% Circular orbits are first assume as only altitude is interesting for
% disturbances computation. Elliptic orbits may be characterized by an
% worst-case and best-case, ranging in the orbit altitude range proposed.

% All units are in S.I.

%% General setup 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration tolerances

%% Initial data 
%Earth characteristics 
mu = 3.986e14;          %Earth gravitational parameter 
Re = 6371e3;            %Earth mean radius

% CosmosSat-1 characteristics 
mass = ;                %Spacecraft mass
I = [];                 %Inertia dyadic
max_distance = [];      %Maximum distance to the center of gravity
max_area = [];          %Maximum exposed area
mag_moment = ;          %Residual magnetic moment of the spacecraft

% Orbit characteristics 
hmax = 800e3;           %Maximum orbit altitude
dh = 1;                 %Altitude step
h = Re:dt:(Re+hmax);    %Orbital altitude range    

% Earth magnetic field characteristics

% Atmosphere model

%% Orbit integration 
% Time span 
tmax = 1e3;             %Maximum integration time (s)
dt = 1;                 %Time step (s)
tspan = 0:dt:tmax;      %Integration span (s)

%% Orbital and attitude disturbances computation

%% Results