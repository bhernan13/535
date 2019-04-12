function dp = pressuredrop(mdot,IDpipe,dpeb,press,temp,fluid)

% Created: February 17, 2019
% Code Author: Madie Melcer
% Modified: Brian Hernan April 12, 2019
% Streamlined inputs for use in main pebble bed analysis loop.
% Removed 'global conv' conversions and switched to convert.m
% Call pebpack within loop for consistency and fewer inputs
%
% Code to calculate the pressure drop across a slice of our pebble beds
% given inputs of mass flow rate, porosity, pipe diameter, pebble diameter,
% time step, ethane velocity, ethane density, ethane viscosity, and number
% of pebbles per layer. Uses pressure drop correlation from paper 'Analysis of Pressure
% Drop and Heat Transfer of a Pebble-Bed-Storage Heater for a Hypersonic
% Wind Tunnel' by D. E. Randall and A. Bedford in 1961 from Sandia. 
% Inputs:
% mdot - massflow, kg/s
% porisity - portion of pipe not filled with pebbles, void volume/volume
% dpipe - Inner Pipe Diameter, in
% dpeb - Pebble Diameter, in
% Press - Local Pressure, Pa
% Temp - Local Temperature, K
% fluid - Refpropm library string of fluid and .mix or .fld suffix
%
% Outputs:
% dp - Pressure drop in psi over a given section of the bed
%
%% Setup
g = 32.2; % ft/s^2
IDpipe = convert.length(IDpipe,'m','ft');
dpeb = convert.length(dpeb,'m','ft'); % ft
rho_eth = convert.density(refpropm('D','T',temp,'P',convert.pressure(press,'pa','kPa'),fluid),'kg/m3','lbm/ft3');
mu_eth = convert.dynamicvisc(refpropm('V','T',temp,'P',convert.pressure(press,'pa','kPa'),fluid),'Pa-s','lb/ft-s');
w = convert.mass(mdot,'kg','lbm');

%% Calculations
dx = dpeb;
A = pi/4*IDpipe^2; % ft^2
U =  w / (rho_eth*A); % ft/s
[npeb_in_slice, porosity, ~ , vol_bed] = pebpack(dpeb,IDpipe);
% 1/ft, area of pebbles per unit volume of bed
S = npeb_in_slice*4*pi*(dpeb/2)^2/vol_bed; 
% Reynolds number for pebble bed based on pebble surface area per unit volume of bed
Re2 = rho_eth*U/(mu_eth*S); 
% Pressure drop over given section
dp = 2.4*(Re2^-0.1)*((1-porosity)/porosity)*(w^2/(rho_eth*g*dpeb*(A*porosity)^2)*dx);

