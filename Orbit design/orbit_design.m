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
tau = (3600*24)*(365.242199/(1+365.242199)); 
Earth_Omega = (2*pi)/tau;

%Orbit requirements 
repeating_days = 1; 