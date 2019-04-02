function [ d_throat, P_set, v ] = cavVenturi( P1, T, mdot, fluid, area, setpressure, P_min )
% Function for calculating throat size of cavitating venturis. Currently
% takes metric units only. Conversions for english added later.

P_vap = refpropm('P', 'T', T, 'Q', 1, fluid)*1e3;
rho = refpropm('D', 'T', T, 'P', P1/1e3, fluid);
Cd = 0.97;
dP = P1 - P_vap;

if P_min == 0
    P_min = 543*6894.76; % vapor pressure of fluid - here it's ethane
end

A_throat = mdot ./ (sqrt(2.*rho.*dP));

d_throat = real(sqrt((4.*A_throat)./pi));

if strcmp(setpressure,'yes')
    P_set = ( (mdot./(A_throat.*Cd) ).^2 .* (1./(2.*rho)) ) + P_vap;
    if P_set < P_min
        indexPlow = find(P_set <= P_min);
        % Calculate dP as the difference between vapor pressure and
        % saturation pressure for a saturate liquid
        dP_adjust = P_min - P_vap;
        % Calculate new throat area using minimum set pressure
        A_throatadjust = mdot ./ sqrt(2.*rho(indexPlow).*dP_adjust);
        % Store new throat areas in 
        A_throat(indexPlow) = A_throatadjust;
    end
        
else
    P_set = NaN;
end

v = mdot ./ (area.*rho);

end

