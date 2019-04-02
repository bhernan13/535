function varargout = pebblBed(dx,dt,bedLength,runTime,properties,inletPress,...
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
P_fluid(:,:) = inletPress;
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
    areaDepth = 1;
end
%% Pebble Bed 1 Analysis
fprintf('\nBeginning Bed 1 Analysis...\n');
tic
for i = 1:length(t)
%     disp('Hi, meet me, Bug.');
   for j = 1:length(x)
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
       if material == strcmpi(material,'SS304')
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
fprintf('\nFinished Bed 1 Analysis...\n');
toc
%% Post Process Results
% Move post process, plotting and write to a dedicated post processing
% function
if options{1,2}
    fprintf('\nPost Processing...\n');
    % Final exit temp of fluid
    tempExit = T_fluid(end,end)/conv.r2k;

    % Location where fluid equals target temeprature
    tempLine = zeros(length(t),1);
    tempFluid = zeros(length(t),1);
    tempBed = zeros(length(t),1);
    for i = 1:length(t)
        if T_fluid(end,i) >= floor(Setpoint*conv.r2k)
            tempLine(i,1) = find(T_fluid(:,i) >= floor((Setpoint*conv.r2k)),1);
        elseif T_fluid(end,i) < floor(Setpoint*conv.r2k)
            tempLine(i,1) = length(T_fluid(:,i));
        end
        tempFluid(i,1) = T_fluid(tempLine(i),i);
        tempBed(i,1) = T_peb(tempLine(i)-1,i+1);

    end

    if tempExit < floor(Setpoint)
        tempTarget = 'Not Achieved';
    elseif tempExit >= floor(Setpoint)
        tempTarget = tempLine(end)*(dx/conv.in2m);
    end
end
%% Plotting
% Move post process, plotting and write to a dedicated post processing
if options{2,2}
    fprintf('\nPlotting results...\n');
    figure(1);
    grid on; hold on;
    set(gca,'xdir','reverse','ydir','reverse')
    plot(fliplr(tempFluid),tempLine*dx / conv.in2m);
end
%% Write Input Conditions and results to file
% Move post process, plotting and write to a dedicated post processing
if options{3,2}
    fprintf('\nWriting Results...\n');
    % Check if Results Directory exists, if not, make it.
    currentDir = pwd;
    if isfolder(strcat(currentDir,'\Results'))
    else
        mkdir('Results');
    end
    cd(strcat(currentDir,'\Results'));
    
    % Check if results file exists. If yes, open with append permissions. If
    % not, create file with append permission and add headers.
    fileName = strcat(fluid,'_',material,'.txt');
    if exist(fileName,'file') == 2
        fileID = fopen(fileName,'a');
    else
        fileID = fopen(fileName,'a');
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\n','Pebble Material','Pipe ID [in]','Pipe Length [in]',...
            'Set Temp [R]','Final Exit Temp [R]',strcat(fluid,'Temp=Target @ [in]') );
    end
    
    % Print results from current run to results file
    fprintf(fileID,'%-15s\t%-12s\t%-16s\t%-12s\t%-16s\t%-s\n',material,num2str(bedID/conv.in2m),...
        num2str(bedLength/conv.in2m),num2str(Setpoint),num2str(tempExit),num2str(tempTarget));

    fclose(fileID);
    % Change back to main directory
    cd(currentDir);

    fprintf('\nDone\n');
end
%% Gather Requested Outputs
if isempty(outputs)
    fprintf('No reults requested.');
    varargout = {[]};
else
    varargout = cell(length(outputs),2);
    varargout{1,1} = 'Time [s]';
    varargout{1,2} = t;
    varargout{2,1} = 'Distance [m]';
    varargout{2,2} = x;
    indexOut = 3;
    if any(contains(outputs(:,1),'temp','IgnoreCase',true))
        varargout{indexOut,1} = 'Temperature Fluid [K]';
        varargout{indexOut,2} = T_fluid;
        indexOut = indexOut+1;
        varargout{indexOut,1} = 'Temperature Pebbles [K]';
        varargout{indexOut,2} = T_peb;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'pres','IgnoreCase',true))
        varargout{indexOut,1} = 'Pressure [Pa]';
        varargout{indexOut,2} = P_fluid;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'vel','IgnoreCase',true) )
        varargout{indexOut,1} = 'Velocity [m/s]';
        varargout{indexOut,2} = V;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'heat','IgnoreCase',true)) || any(contains(outputs(:,1),'htc','IgnoreCase',true))
        varargout{indexOut,1} = 'Film Coefficient [W/m2-K]';
        varargout{indexOut,2} = hf;
        indexOut = indexOut+1;
        varargout{indexOut,1} = 'Heat Transfer Coefficient [W/m2-K]';
        varargout{indexOut,2} = H;
        indexOut = indexOut+1;
    end
    if any(contains(outputs(:,1),'input','IgnoreCase',true)) 
        varargout{indexOut,1} = 'Inputs';
        varargout{indexOut,2} = inputs;
        indexOut = indexOut+1;
    end
%     else
%         warning('Output requested does not exist.');   
%     disp('Hi, just checking, in'); 
end


end


