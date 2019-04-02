Functions%***********************************
%Critical Flow Calculation Subroutine
%Created: April 16, 2004 Purdue University
%Author: Matthew Long
%Last Modified: April 16, 2004 Purdue University
%Last Modified By: Matthew Long
%Description:  This subroutine calculates the mass flow rate for a given
%critical flow (sonic) nozzle based on compressible flow equation found in
%"R.W. Miller, 'Flow Measurement Engineering Handbook', 2nd Ed, 1989" as
%This is the source used in flow calculations by Stuart Cole at Flowmaxx 
%Engineering
%
%Inputs:
%R=Specific Gas Constant: Ru/M, ft*lbf/lbm*R
%k=isentropic exponent* (*note for pressures above 100 psi, gamma is NOT a
%   good estimate of k).
%Pc=Critical Pressure of gas, psi
%Tc=Critical Temperature of gas, R
%P1=Inlet Static Pressure, psi
%T1=Inlet Static Temperature, R
%Cd=Critical Flow Nozzle discharge coefficient based on curve fit from 
%   nozzle geometry
%d1=inlet diameter, in
%d2=nozzle throat diameter, in
%
%Output:
%m=mass flow rate, lbm/sec
function [m]=critical_flow(R,k,Pc,Tc,P1,T1,Cd,d1,d2)

%****Calculate Compressibility factor Z
P2=P1.*(1+(k-1)./2).^(-k./(k-1)); %Throat Static Pressure
rc=P2./P1; %Critical Pressure Ratio
Pr=P1./Pc; %Reduced Pressure of gas
Tr=T1./Tc; %Reduced Temperature of gas
B0=0.083-0.422./Tr.^1.6; %First Virial coefficient
Z=1+B0.*Pr./Tr; % Compressiblity Factor

%Calculate compressibility factor for Gox
if (R>48.4 & R<49),
    P=linspace(1000,2000,11);
    Z1=[0.96067 0.95766 0.95486 0.95227 0.94988 0.94772 0.94579 0.94409 0.94262 0.94138 0.94039];
    [C,S]=polyfit(P,Z1,2);
    Z=polyval(C,P1);
end
%****Calculate Compressibility factor Z

beta=d2/d1; %Beta
gc=32.174; %Constant lbm*ft/lbf*s^2
A2=0.25*pi*d2^2; %in^2 Throat area

num=2.*gc.*rc.^(2./k).*(k./(k-1)).*(1-rc.^((k-1)./k)); 
den=Z.*R.*T1.*(1-beta.^4.*rc.^(2./k));

m=Cd.*A2.*P1.*(num./den).^0.5; %lbm/sec Mass flow rate