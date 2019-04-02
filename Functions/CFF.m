function [ CFF ] = CFF( T, P, mdot, Area, fluid )
% Function to calculate the real gas critical flow factor for sonic nozzle
% venturis. Equations taken from "Real Gas Effects and Critical Flow
% Through Nozzles and Tabulated Thermodynamic Properties" by Robert C.
% Johnson, NASA Lewis Research Center, Cleveland, OH
% Code Author: Brian Hernan
% Created: January 21, 2019
%% Unit Conversions
conv_in_m = 0.0254;

%%
Ru = 8314;

%% Obtain Preoperties from REFPROP
rho = refpropm('D', 'T', T, 'P', P/1e3, fluid);    
a = refpropm('A', 'T', T, 'P', P/1e3, fluid); 
gamma = refpropm('K', 'T', T, 'P', P/1e3, fluid);
MW = refpropm('M', 'T', T, 'P', P/1e3, fluid);
%% Flow velocity and mach
v = mdot / (rho*Area*(conv_in_m^2));
M = v/a;

%% Stagnation properties upstream of nozzle
Pt_venturi = P * ( 1 + ((gamma-1)/2) * M^2) ^ (gamma / (gamma-1));
Tt_venturi = T * ( 1 + ((gamma-1)/2) * M^2 );

%% Calculate Mass flowrate per unit area
G = rho*v;

%% Calculate Critical Flow Factor (CFF or Cstar)
CFF = (G * sqrt( (Ru/MW) * Tt_venturi) ) / Pt_venturi;
end

