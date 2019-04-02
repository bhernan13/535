function dp = pressuredrop(mdot,porosity,dpipe,dpeb,dt,velocity,density_eth,viscosity_eth,npeb)

% Created: February 17, 2019
% Code Author: Madie Melcer
% Code to calculate the pressure drop across a slice of our pebble beds
% given inputs of mass flow rate, porosity, pipe diameter, pebble diameter,
% time step, ethane velocity, ethane density, ethane viscosity, and number
% of pebbles per layer. Uses pressure drop correlation from paper 'Analysis of Pressure
% Drop and Heat Transfer of a Pebble-Bed-Storage Heater for a Hypersonic
% Wind Tunnel' by D. E. Randall and A. Bedford in 1961 from Sandia. 
% Inputs are in metric units.
%% Load conversions
global conv

%% setup
D = dpipe/conv.ft2m; % ft
d = dpeb/conv.ft2m; % ft
rho_eth = density_eth*conv.kgm32lbft3; % lb/ft^3
mu_eth = viscosity_eth*conv.Pas2lbfts; % lb/ft-s
g = 32.2; % ft/s^2

%% calculations
dx = dpeb; % ft
A = pi/4*D^2; % ft^2
w = mdot/conv.lbm2kg; % lb/s, mass flow rate
U = velocity/conv.ft2m; % ft/s
vol_bed = pi*(D/2)^2*dx; % for this slice, ft^3
npeb_in_slice = npeb*(dx/d);
S = npeb_in_slice*4*pi*(d/2)^2/vol_bed; % 1/ft, area of pebbles per unit 
    % volume of bed
Re2 = rho_eth*U/(mu_eth*S); % Reynolds number for pebble bed based on 
    % pebble surface area per unit volume of bed
dp = 2.4*(Re2^-0.1)*((1-porosity)/porosity)*(w^2/(rho_eth*g*d*(A*porosity)^2)*dx);

