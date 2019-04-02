function [ Dt_sonic_eng, v, M, Pt_venturi ] = Sonic_Venturi( T, P, mdot, Area, fluid )
% Author: Brian Hernan 
% Created: February 2018
% Function to size sonic venturis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Inputs: 
% Gas temperature - T [K]
% Gas pressure - P [Pa]
% Gas mass flowrate - mdot [kg/s]
% Criticial Flow Factor  - Cstar - refprop output in code
% Tube cross sectional area - Area [in^2]
% Fluid - used to get REFPROP properties. Enter as a string.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Outputs:
% Venturi throat diameter [in] - Dt_venturi_eng
% Gas velocity - v [m/s]
% Gas mach - M 

% Unit Conversions
conv_in_m = 0.0254;
conv_m_in = 1/0.0254;

Ru = 8314;

% Get gas Properties from REFPROP. Pressure entered in kPa. Used in
% sonic venturi sizing.
for i = 1:length(P)
    rho = refpropm('D', 'T', T, 'P', P(i)/1e3, fluid);    
    a = refpropm('A', 'T', T, 'P', P(i)/1e3, fluid); 
    gamma = refpropm('K', 'T', T, 'P', P(i)/1e3, fluid);
    MW = refpropm('M', 'T', T, 'P', P(i)/1e3, fluid);
    Cstar = refpropm('~','T',T,'P',P(i)/1e3, fluid);
end
    
% Flow velocity and mach
v = mdot / (rho*Area*(conv_in_m^2));
M = v/a;

% Stagnation properties upstream of venturi
Pt_venturi = P * ( 1 + ((gamma-1)/2) * M^2) ^ (gamma / (gamma-1));
Tt_venturi = T * ( 1 + ((gamma-1)/2) * M^2 );

% Calculate Sonic venturi throat size
At_sonic = mdot * sqrt(((Ru/MW)*Tt_venturi)) / (Pt_venturi*Cstar);
Dt_sonic = sqrt((4*At_sonic) / pi);
Dt_sonic_eng = Dt_sonic*conv_m_in;




