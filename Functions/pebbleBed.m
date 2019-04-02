function Outputs = pebbleBed(dx,dt,bedLength,runTime,properties,inletPress,...
              inletTemp,Setpoint,mdot,porosity,pebbleDiameter,...
              bedID,N,fluid,material,options,outputs)
% Code Author: Brian Hernan
% Date: February 26, 2019
% Code to calculate sizing and performance of a pebble bed heater.Fixed 
% delta t, and fixed delta x, with varying numbers of dx steps advanced per
% time step (will be accounted for by factoring in velocity and pressure 
% chagnes). Pressure drop is accounted for. Determines the number of 
% spheres (balls,pebbles) per layer in the pipe based on how many circles 
% of a given diameter can fit into a larger circle. Assumed that each layer 
% of circles are stacked discretly and no pebble from a previous layer 
% the layer above or below. The boundary of the pebble bed pipe will be
% affects treated as adiabatic and at a constant set temeprature. Pebble 
% material properties vary as a function of temperature based off of the 
% best available data. Fluid properties will be obtained at relevant 
% conditions from REFPROP or from imported NIST tables.
%% Stop debug 
dbstop if error
%% Load Conversions
conversions;
%% Store inputs
inputs{1,1} = 'nargin'; inputs{1,2} = nargin;
inputs{2,1} = 'nargout'; inputs{2,2} = nargout;
inputs{3,1} = 'Fluid'; inputs{3,2} = fluid;
inputs{4,1} = 'Pebble Material'; inputs{4,2} = material;
inputs{5,1} = 'Pebble Diameter [in]'; inputs{5,2} = pebbleDiameter/conv.in2m;
inputs{6,1} = 'Bed Length [in]'; inputs{6,2} = bedLength/conv.in2m;
inputs{7,1} = 'Bed ID [in]'; inputs{7,2} = bedID/conv.in2m;
inputs{8,1} = 'Inlet Temp [F]'; inputs{8,2} = (inletTemp/conv.r2k) - conv.f2r;
inputs{9,1} = 'Inlet Pressure [psi]'; inputs{9,2} = inletPress/conv.psi2pa;
inputs{10,1} = 'Bed Setpoint [F]'; inputs{10,2} = (Setpoint/conv.r2k)-conv.f2r;
inputs{11,1} = 'Mass Flow [lbm/s]'; inputs{11,2} = mdot / conv.lbm2kg;
inputs{12,1} = 'Test Duration'; inputs{12,2} = runTime;
%% X and t vectors
x = [0:dx:bedLength];
t = [0:dt:runTime];
%% Initialize Vectors
q_fluid = zeros(length(x),length(t));
q_peb = zeros(length(x),length(t));
rho = zeros(length(x),length(t));
Re1 = zeros(length(x),length(t));
Pr = zeros(length(x),length(t));
k_fluid = zeros(length(x),length(t));
hf = zeros(length(x),length(t));
k_peb = zeros(length(x),length(t));
H = zeros(length(x),length(t));
T_fluid = zeros(length(x),length(t));
T_peb = zeros(length(x),length(t));
V = zeros(length(x),length(t));
mu = zeros(length(x),length(t));
Cp_fluid = zeros(length(x),length(t));
P_fluid = zeros(length(x),length(t));
Cp_peb = zeros(length(x),length(t));
dP = zeros(length(x),length(t));
%% Select Properties for Chosen Pebble Material
materialSelect = strcmpi(properties(:,1),material);
if materialSelect(:) == 0
    error('Material data is not entered for chosen material. Choose a material with existing data or enter new material data.');
end
indexMaterial = find(materialSelect ~= 0);
Data_k = properties{indexMaterial,2};
Data_Cp = properties{indexMaterial,3};
Data_rho = properties{indexMaterial,4};
%% Set initial values
i = 1; j = 1;
P_fluid(i,:) = inletPress;
T_fluid(j,:) = inletTemp;
T_peb(:,i) = Setpoint;
m_slug = mdot*dt;
layerVolume = (pi/4)*(bedID^2)*pebbleDiameter;
xsectionArea = (pi/4)*(bedID^2);
m_peblayer = (1-porosity)*layerVolume*Data_rho;
%Surface area per unit depth, where dx is a unit depth. 
if options{4,2}
    areaDepth = dx*((6*xsectionArea*(1-porosity)) / pebbleDiameter);
else
    areaDepth = (bedLength/pebbleDiameter)*(4*0.25*pi*pebbleDiameter^2)*N; %0.9399;
end
%% Pebble Bed 1 Analysis
% fprintf('\nBeginning Bed 1 Analysis...\n');
tic
for i = 1:length(t)
%     disp('Hi, meet me, Bug.');
% if i == i
%     pause
% end
   for j = 1:length(x)
%        if T_fluid(j,i) < 300
%            pause
%        end
       % Density of ethane 
       rho(j,i) = refpropm('D','T',T_fluid(j,i),'P',P_fluid(j,i)/1e3,fluid);
       % Viscocity of ethane
       mu(j,i) = refpropm('V','T',T_fluid(j,i),'P',P_fluid(j,i)/1e3,fluid);
       % Prandtl number of ethane
       Pr(j,i) = refpropm('^','T',T_fluid(j,i),'P',P_fluid(j,i)/1e3,fluid);
       % Thermal conductivity of ethane at current conditions. 
       k_fluid(j,i) = (refpropm('L','T',T_fluid(j,i),'P',P_fluid(j,i)/1e3,fluid));
       % Thermal conductivity of pebbles at current temp. Interpolated from
       % table data. Interpolation between data points. No extrapolation.
       % If temperature is above upper limit, value is fixed at
       % corresponding, same for if temeprature is lower
       if T_peb(j,i) > max(Data_k(1,:))
           k_peb(j,i) = Data_k(2,end);
       elseif T_peb(j,i) < min(Data_k(1,:))
           k_peb(j,i) = Data_k(2,1);
       else
           k_peb(j,i) = interp1(Data_k(1,:),Data_k(2,:),T_peb(j,i),'linear');
       end
       % Specific heat of ethane
       Cp_fluid(j,i) = refpropm('C','T',T_fluid(j,i),'P',P_fluid(j,i)/1e3,fluid);
       % Specific heat of Pebble Material. Curve fit for SS304.
       % Interpolation between data points for all other materials. If
       % temperature exceeds that max temperature for data the value is set
       % to the value corresponding to the max temperature of the data. No
       % extrapolation.
       if strcmpi(material,'SS304')
           Cp_peb(j,i) = Data_Cp(T_peb(j,i));
       else
           if T_peb(j,i) > max(Data_Cp(1,:))
               Cp_peb(j,i) = Data_Cp(2,end);
           elseif T_peb(j,i) < min(Data_Cp(1,:))
               Cp_peb(j,i) = Data_Cp(2,1);
           else
               Cp_peb(j,i) = interp1(Data_Cp(1,:),Data_Cp(2,:),T_peb(j,i),'linear');
           end
       end
       % Current ethane velocity
       V(j,i) = mdot / (rho(j,i)*xsectionArea);
       % Reynolds number for film heat transfer of the pebble
       Re1(j,i) = (rho(j,i)*V(j,i)*pebbleDiameter) / mu(j,i);
       % Film heat transfer coefficient of pebbles
       hf(j,i) = (0.56*((k_fluid(j,i)*conv.km2keng)/(pebbleDiameter/conv.in2m))*(Re1(j,i)^0.6)*(Pr(j,i)^0.33)) / conv.hm2heng;
       % Heat transfer coefficient from pebbles to ethane
       H(j,i) = 1 / ( (1/hf(j,i)) + ( pebbleDiameter / (2*k_peb(j,i)) ) );
       % Calculate residence time of the fluid slug 
       % Ethane temperature exiting the current layer at current timestep
       T_fluid(j+1,i) = T_fluid(j,i) + ( (H(j,i)*dt*areaDepth) / (m_slug * Cp_fluid(j,i)) ) * (T_peb(j,i) - T_fluid(j,i));
       % Pebble temperature of each layer at current timesstep
       T_peb(j,i+1) = T_peb(j,i) - ( (H(j,i)*dt*areaDepth) / (m_peblayer * Cp_peb(j,i)) ) * (T_peb(j,i) - T_fluid(j+1,i));
       % just a pause for debugging
%        disp('Hi, Bug here, not leaving.');
       % Pressure drop through given layer
       dP(j,i) = (pressuredrop(mdot,porosity,bedID,pebbleDiameter,dt,V(j,i),rho(j,i),mu(j,i),N))*conv.psi2pa;
       % Pressure 
       P_fluid(j+1,i) = P_fluid(j,i) - dP(j,i);
   end
end
% fprintf('\nFinished Bed 1 Analysis...\n');
toc
%% Gather Requested Outputs
if isempty(outputs)
    fprintf('No reults requested.');
    Outputs = {'No outputs requested'};
else
    Outputs = cell(length(outputs),2);
    Outputs{1,1} = 'Time [s]';
    Outputs{1,2} = t;
    Outputs{2,1} = 'Distance [m]';
    Outputs{2,2} = x;
    indexOut = 3;
    if any(contains(outputs(:,1),'temp','IgnoreCase',true))
        Outputs{indexOut,1} = 'Temperature Fluid [K]';
        Outputs{indexOut,2} = T_fluid;
        indexOut = indexOut+1;
        Outputs{indexOut,1} = 'Temperature Pebbles [K]';
        Outputs{indexOut,2} = T_peb;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'pres','IgnoreCase',true))
        Outputs{indexOut,1} = 'Pressure [Pa]';
        Outputs{indexOut,2} = P_fluid;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'vel','IgnoreCase',true) )
        Outputs{indexOut,1} = 'Velocity [m/s]';
        Outputs{indexOut,2} = V;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'heat','IgnoreCase',true)) || any(contains(outputs(:,1),'htc','IgnoreCase',true))
        Outputs{indexOut,1} = 'Film Coefficient [W/m2-K]';
        Outputs{indexOut,2} = hf;
        indexOut = indexOut+1;
        Outputs{indexOut,1} = 'Heat Transfer Coefficient [W/m2-K]';
        Outputs{indexOut,2} = H;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'input','IgnoreCase',true)) 
        Outputs{indexOut,1} = 'Inputs';
        Outputs{indexOut,2} = inputs;
        indexOut = indexOut+1;
    end
end


end


