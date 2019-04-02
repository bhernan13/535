function [ output ] = regulator_Sizing( phase, units, fluid, solve, P1, P2, T, Var4   )
%% Regulator and Valve Cv and Flowrate Sizing
% Equations taken from "Tescom Flow Formulas for Computing Gas and Liquid
% Flow Through Valves and Regulators" Reference.
% Code Author: Brian Hernan
% Created: January 19,2019
%% Function Inputs
%
% phase - Input Options: liquid, gas
%      Descritpion: String describing the phase of the fluid (liquid or gas)
% units - Input Options: english, metric, volumetric, standard
%      Description: Which system will units be entered in. 
%      Units must be entered in base units:
%      Metric, kg/s, pa
%      English: lbm/s, psi
%      VolumetricE: ft3/s, psi - To be added later
%      VolumetricM: m3/s, pa - To be added later
%      Standard: SCFM, psi - To be added later
% fluid - Input Options: Any REFPROP fluid or mixture
%      Description: strong containing which fluid is being used
% solve - Input Options: Coefficient or Flowrate
%      Description: String detemrines if solving for flow coefficient (Cv)
%      or for flowrate (Ql, Qg, Q, M). Flowrate input solves flowrate,
%      Coefficient input solves Cv.
% P1 - upstream regulator pressure, units must match above
% P2 - downstream regulator pressure, units must match above
% T - Temperature fluid is expected to be at. Most likely ambient (294K)
% Var4 - Input flowrate (in above one of the above unit systems) if solving
% for flow coefficient. Input flow coefficient if solving for flowrate
% 
%% Function Outputs
% 
% Output -  Outputs either flowrate (kg/s) or flow coefficient, Cv
%% Variables Definition
% Names and definitions correspond to the Tescom Flow Formulas reference
%
% C_v - Flow coefficient for regulators and valves that expresses flow 
% capabilities of a unit at full open condition. For liquids, this 
% coefficient is defined as the flow of water at 60? F in gallons per
% minute at a pressure drop of one psig. For gases, this coefficient is 
% defined as the flow of air at standard conditions in standard cubic feet 
% per minute for each psig of inlet pressure.
%
% S_L - Specific gravity of liquids relative to water, both at a standarf
% temperature of 60F. (Specific gravity of water = 1.0 @ 60F)
%
% S_g - Specific gravity of a gas relative to air; equals the ratio of the
% molecular weight of the gas to that of air. 
%(Specific gravity of air = 1.0 @ 60F)
%
% psia - Absolute pressure which is gauge pressure (psig) plus 14.7
% (atmospheric)
%
% P - Line pressure (psia)
% P1 - Inlet Pressure (psia)
% P2 - Outler pressure (psia)
% dP - Differential pressure (P1-P2, psia)
% Q_L - Liquid flow in gallons per minute(GPM)
% Q_g - Gas flow in standard cubic feet per minute(SCFM). At standard
% conditions of 60F and 14.7 psia.
% Q - Volumetric flowrate in cubic feet per minute (CFM)
% M - Mass flow rate in pounds per minute(lbs/min)

%% Conversions
conv_m3_ft3 = 35.3147; % ft3 per m3
conv_ft3_GAL = 7.48052; % ft3 per GAL
conv_kg_lbm = 2.20462; % lbs per kg
conv_psi_pa = 6894.76; % Pa per psi
conv_kgm3_slugft3 = 0.00194032; % slug ft3 per kg/m3

switch units
    case 'metric'
    P1 = P1/conv_psi_pa; % psi
    P2 = P2/conv_psi_pa; % psi
    dP = P1-P2; % psi
    if dP < 0
        error('Your downstream pressure is greater than your upstream pressure P2>P1. Please check your input pressures. Check Line 75 and inputs #5 and #6. :)');
    end
    if strcmp(solve,'coefficient')
        Density = ( refpropm('D','T',T,'P',P1*conv_psi_pa/1e3,fluid) ) * conv_kg_lbm / conv_m3_ft3; % lbm/ft3
        M = Var4 * conv_kg_lbm * 60; % lbm/min
    else
        Density = refpropm('D','T',T,'P',P1*conv_psi_pa/1e3,fluid) * conv_kg_lbm / conv_m3_ft3; % lbm/ft3
        Cv = Var4;
    end
    case 'english'
    P1 = P1; % psi
    P2 = P2; % psi
    dP = P1-P2; % psi
    if dP < 0
        error('Your downstream pressure is greater than your upstream pressure P2>P1. Please check your input pressures. Check Line 90 and inputs #5 and #6. :)');
    end
    if strcmp(solve,'coefficient')
        Density = refpropm('D','T',T,'P',P1*conv_psi_pa/1e3,fluid) * conv_kg_lbm / conv_m3_ft3; % lbm/ft3
        M = Var4 * 60; % lbm/min
    else 
        Density = refpropm('D','T',T,'P',P1*conv_psi_pa/1e3,fluid) * conv_kg_lbm / conv_m3_ft3; % lbm/ft3
        Cv = Var4; 
    end        
end

% Reference Point
Tref = 294; %Kelvin
Pref = 101325; % Pa

switch phase
    % Cv or Flowrate calculation for liquids
    case 'liquid'
        Density_H2O = refpropm('D','T',Tref,'P',Pref/1e3,'water') * conv_kg_lbm / conv_m3_ft3;
        S_L = Density / Density_H2O;
        switch solve
            case 'flowrate'
                Q_L = ( Cv*sqrt(dP) ) / sqrt(S_L); % SCFM
                Q_L_h2o = Q_L / (sqrt(1/S_L)); % SCFM of water
                M = ( Q_L / (conv_ft3_GAL * 60) ) * Density; % lbm/s
                output = {Q_L, 'SCFM', fluid;
                          Q_L_h2o, 'SCFM', 'water';
                          M, 'lbm/s', fluid};
            case 'coefficient'
                Q_L = ( M / Density ) * conv_ft3_GAL; % SCFM
                Q_L_h2o = Q_L / (sqrt(1/S_L)); % SCFM of water
                Cv = (Q_L * sqrt(S_L) ) / sqrt(dP);   
                output = {'SCFM', Q_L, fluid;
                          'SCFM', Q_L_h2o, 'water';
                          'Flow Coefficient', Cv, fluid};
        end
    case 'gas'
        Density = refpropm('D','T',T,'P',P1,fluid);
        MM_air = refpropm('M','T',Tref,'P',Pref,'air.mix');
        MM = refpropm('M','T',Tref,'P',Pref,fluid);
        S_g = MM/MM_air;
        switch solve
            case 'flowrate'
                if P1 >= 2*P2
                    Q_g = ( Cv*P1 ) / ( 2*sqrt(S_g) ); % SCFM
                    Q = Q_g*14.7 / P1; % CFM
                    M = Q_g * (1/60) * Density; % lbm/s
                    Q_g_air = Q_g / (sqrt(1/S_g)); % SCFM of air
                elseif P1 < 2*P2
                    Q_g = ( Cv * sqrt((dP*P2)) ) / ( sqrt(S_g) ); % SCFM
                    Q = Q_g*14.7 / P1; % CFM
                    M = Q_g * (1/60) * Density; % lbm/s
                    Q_g_air = Q_g / (sqrt(1/S_g)); % SCFM of air
                end
                output = {Q_g, 'SCFM', fluid;
                          Q, 'CFM', fluid;
                          M, 'lbm/s', fluid;
                          Q_g_air, 'SCFM', 'air'}; 
            case 'coefficient'
                if P1 >= 2*P2
                    Q_g = (M / Density) * conv_ft3_GAL; % SCFM
                    Q = Q_g*14.7 / P1; % CFM
                    Q_g_air = Q_g / (sqrt(1/S_g)); % SCFM of air
                    Cv = ( Q_g * 2*sqrt(S_g) ) / P1; 
                elseif P1 < 2*P2
                    Q_g = (M / Density) * conv_ft3_GAL; % SCFM
                    Q = Q_g*14.7 / P1; % CFM
                    Q_g_air = Q_g / (sqrt(1/S_g)); % SCFM of air
                    Cv = ( Q_g*sqrt(S_g) ) / ( sqrt(dP*P2) ); 
                end
                output = {Q_g, 'SCFM', fluid;
                          Q, 'CFM', fluid;
                          Cv, 'Flow Coefficient', fluid;
                          Q_g_air, 'SCFM', 'air'};
        end

end

