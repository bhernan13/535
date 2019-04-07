%% Thermal Design and Performance Analysis for Pebble Bed Heaters
% Created: February 15, 2019
% Code Author: Brian Hernan
%% Notes to add
% References for thermal properties - commented at the start of each
% thermal properties section. 
% Pass inlet velocity as an input based on velocity coming out of the inlet
% manifold. Pass inlet velocities to downstream beds based on exit/inlet
% conditon
%%
% Bed 1, 50" long 3.4" ID, bed 2, 18" long 3.4" ID both SS304 pebbles
clear all;close all;clc; %#ok<CLALL>
dbstop if error;
%% Add Paths to Functions and Refprop
% Refprop folder must be stored on the 'C:\temp' filepath. 
% ECN drives prevent running refprop from within them. 
cd('C:\Users\gymna\Documents\Purdue\AAE 535\github\535-master');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentDir = pwd;
refpropDir = 'C:\temp\REFPROP';
if isfolder(refpropDir)
else 
    error('REFPROP folder is not in correct location. \n Copy REFPROP folder to %s \n','C:\temp');
end 
addpath(genpath(refpropDir),'-end');
addpath(genpath(currentDir),'-end');
addpath(genpath(strcat(currentDir,'\Results')));
%% Conversions
conversions;
global conv
%% General Inputs and Bed 1 Inlet Conditions
fprintf('\nInput General and Initial Inlet Conditions.\n\n');
fluid = char('ethane'); % input('Input REFPROP Fluid (.fld) or Mixture (.mix) with corect suffix: ','s'));
beds = 2; % input('Enter number of beds in series: ');

ethane{1,1} = 'Property'; ethane{1,2} = 'Value';
ethane{2,1} = 'Inlet Temperature [R]'; ethane{3,1} = 'Inlet Pressure [pa]';
ethane{4,1} = 'Outlet Target [R]'; ethane{5,1} = 'Outlet Target Pressure [psi]';

ethane{2,2} = (68 + conv.f2r) * conv.r2k; %K input('Fluid Inlet Temp [F]: ')
ethane{3,2} = input('Fluid Inlet Pressure [psi]: ')*conv.psi2pa; %pa
ethane{4,2} = (input('Final Target Temp [R]: '))*conv.r2k; %R
ethane{5,2} = input('Final Target Pressure [psi]: ')*conv.psi2pa; %pa
solve{4,2} = input('mdot [lbm/s]: '); % lbm/s   THIS STAYS IN lbm/s
solve{5,2} = 5; % input('Run Time [s]: '); % Hotfire, s
%% Bed Specific Inputs
% Initialize
pebbleDiameter = cell(beds+1,2);
pebbleDiameter{1,1} = 'Bed #'; pebbleDiameter{1,2} = 'Pebble Diameter [m]';
material = cell(beds+1,2);
material{1,1} = 'Bed #'; material{1,2} = 'Pebble Material';
bed = cell(10,beds+1);
bed{1,1} = 'Property'; bed{2,1} = 'Length [m]';
bed{3,1} = 'Pipe ID [m]'; bed{4,1} = 'Cross Sectional Area [m2]';
bed{5,1} = 'Bed Volume [m3]'; bed{6,1} = 'Temp Setpoint [K]'; bed{7,1} = 'Pebbles/layer';
bed{8,1} = 'Porosity'; bed{9,1} = 'Pebble Percent [1-porosity]';
bed{10,1} = 'Layer Volume [m3]';

for i = 2:beds+1
    fprintf('\nInput Parameters for Bed %s.\n\n', num2str(i-1));
    material{i,1} = strcat('Bed ',num2str(i-1));
    material{i,2} = char('SS304'); % input('Enter Pebble Material - SS304,SS316,CU201,AL2017,BRASS260: ','s')
    pebbleDiameter{i,1} = strcat('Bed ',num2str(i-1));
    pebbleDiameter{i,2} = 0.25*conv.in2m; % m input('Pebble Diameter [in]: ')
    bed{1,i} = strcat('Bed ',num2str(i-1));
    bed{2,i} = input('Bed Length [in]: ')*conv.in2m; %m, fixed by Kyle.
    bed{3,i} = 3.4*conv.in2m; %m, internal diameter input('Pipe ID [in]: ')
    bed{4,i} = (pi/4).*(bed{3,2}.^2); %m2
    for j = 1:length(bed{2,2})
        bed{5,i}(:,j) = (pi/4).*(bed{3,2}.^2).*bed{2,2}(j); %m3
    end
    bed{6,i} = (input('Bed Set Point [R]: '))*conv.r2k; %K
    [bed{7,i}, bed{8,i}, bed{9,i}, bed{10,i}] = pebpack(pebbleDiameter{i,2}, bed{3,i});
end
%% Material Property Cells
properties = cell(6,4);
properties{1,1} = 'Material'; properties{1,2} = 'T [K];k [W/m-k]';
properties{1,3} = 'T [K]; Cp [J/kg-K]'; properties{1,4} = 'Density [lbm/ft3]';
properties{2,1} = 'SS304'; properties{3,1} = 'SS316';
properties{4,1} = 'CU201'; properties{5,1} = 'AL2017';
properties{6,1} = 'BRASS260';
% Thermal Conductivity 
properties{2,2} = [ 293, 373, 473, 573, 673, 773;
                     16.2, 16.2, 17.5, 18.8, 20.1, 21.4];
properties{3,2} = [293, 373, 473, 573, 673, 773;
                   14, 14.9, 16, 17.3, 18.6, 19.9];
properties{4,2} = [250, 300, 350, 400, 500, 600, 800;
                   406, 401, 396, 393, 386, 379, 366];
properties{5,2} = [250, 300, 400, 500, 600, 700, 800
                   134, 135, 136, 135, 132, 132, 125];
properties{6,2} = [293, 350, 400, 500, 600, 700, 800
                   120, 119, 118, 116, 114, 112, 11];
% Specific Heat - Put everything into cell arrays
properties{2,3} = @(T) 6.683 + 0.04906*T + 80.74*log(T); %{'A + BT + Cln(T)','Equation';'A', 6.683;'B', 0.04906;'C', 80.74};

properties{3,3} = [293, 363, 473, 593, 703, 813, 923, 1033, 1143;
                   452, 486, 528, 548, 565, 573, 586, 615, 649];
properties{4,3} = [300, 400, 600, 800;
                   385, 397, 417, 433];
properties{5,3} = [298, 373, 473, 573, 673, 773, 811;
                   850, 900, 950, 970, 1000, 1080, 1100];
properties{6,3} = [300, 400, 600;
                   380, 395, 425];
% Density
properties{2,4} = 499.392;
properties{3,4} = 499.392;
properties{4,4} = 558.144;
properties{5,4} = 174.528;
properties{6,4} = 532.224;
%% Load Fluid Properties Tables From NIST - Saturated if needed
%% Solve Variables Setup 
solve{1,1} = 'Property'; solve{1,2} = 'Value';
solve{2,1} = 'Time Step, dt [s]'; solve{3,1} = 'X Step,dx [m]';
solve{4,1} = 'Mass Flow [kg/s]'; solve{5,1} = 'Hotfire Duration [s]';
solve{2,2} = 0.02; %s
solve{3,2} = pebbleDiameter; %in
%% Initial Thermal Energy of Each Layer
%% Calculate Bed Entrance Velocity.
% Velocity leaving the manifold
%% Post Processing 
settings{1,1} = 'post process'; settings{2,1} = 'plot'; settings{3,1} = 'write';
settings{4,1} = 'Area';
settings{1,2} = true;
settings{2,2} = true;
settings{3,2} = true;
settings{4,2} = false;
%% Outputs
outputRequest = {'Temperature';'Pressure';'HTC';'Velocity';'Inputs'};
lengthOut = length(outputRequest);
if any(contains(outputRequest(:,1),'temp','IgnoreCase',true))
    lengthOut = lengthOut+1; end
if any(contains(outputRequest(:,1),'htc','IgnoreCase',true)) || any(contains(outputRequest(:,1),'heat','IgnoreCase',true))
    lengthOut = lengthOut+1; end
%% Data Matrix
Data1 = cell(length(bed{3,2})+1,length(bed{2,2})+1);
Data1{1,1} = 'Diameter[in]/Length[in]';
for i = 1:length(bed{2,2})
    Data1{1,i+1} = num2str(bed{2,2}(i) / conv.in2m);
end
for i = 1:length(bed{3,2})
    Data1{i+1,1} = num2str(bed{3,2}(i) / conv.in2m);
end
if beds == 2
Data2 = cell(length(bed{3,3})+1,length(bed{2,3})+1);
Data2{1,1} = 'Diameter[in]/Length[in]';
for i = 1:length(bed{2,3})
    Data2{1,i+1} = num2str(bed{2,3}(i) / conv.in2m);
end
for i = 1:length(bed{3,3})
    Data2{i+1,1} = num2str(bed{3,3}(i) / conv.in2m);
end
end
%% Main Function Call - Bed(s) Sizing
% Requested Outputs are stored in the Output cell array. Each column is
% relevant for each bed. Those values are passed from bed 1 to any beds in
% series as the inlet conditions for the downstream beds. 
for k = 1:beds
    run=0;
    if k == 1
        fprintf('\nStarting Bed 1...\n');
        for i = 1:length(bed{2,2})
            for j = 1:length(bed{3,2})
                Nruns = length(bed{2,2}) * length(bed{3,2});
                [Data1{j+1,i+1}] = pebbleBed(pebbleDiameter{k+1,2},solve{2,2},bed{2,2}(i),solve{5,2},properties,ethane{3,2},...
                                        ethane{2,2},bed{6,2},solve{4,2},bed{8,2}(j),pebbleDiameter{k+1,2},...
                                        bed{3,2}(j),bed{7,2}(j),fluid,material{k+1,2},settings,outputRequest);
                run = run+1;
                fprintf('\n Run %s/%s for Bed 1 completed. %s runs remaining.\n',num2str(run),num2str(Nruns),num2str(Nruns-run));
            end
        end
    end
    
    if k > 1
        fprintf('\nStarting Bed 2...\n');
        for i = 1:length(bed{2,3})
            for j = length(bed{3,3})
                Nruns = length(bed{2,3}) * length(bed{3,3});
                [Data2{j+1,i+1}] = pebbleBed(pebbleDiameter{k+1,2},solve{2,2},bed{2,3}(i),solve{5,2},properties,Data1{2,2}{5,2}(end,:),...
                                        Data1{2,2}{3,2}(end,:),bed{6,3},solve{4,2},bed{8,3}(j),pebbleDiameter{k+1,2},...
                                        bed{3,3}(j),bed{7,3}(j),fluid,material{k+1,2},settings,outputRequest);
                run = run+1;
                fprintf('\n Run %s/%s for Bed 2 completed. %s runs remaining.\n',num2str(run),num2str(Nruns),num2str(Nruns-run));
            end
        end
    end
    
end
%% Post Process Results
% Move post process, plotting and write to a dedicated post processing
% function

% Bed 1
if settings{1,2}
    fprintf('\nPost Processing...\n');
    % Final exit temp of fluid
    tempExit1 = Data1{2,2}{3,2}(end,end);
    setTemp1 = (Data1{2,2}{9,2}{10,2} + conv.f2r)*conv.r2k;
    % Location where fluid equals target temeprature
    tempLine1 = zeros(length(Data1{2,2}{1,2}),1);
    tempFluid1 = zeros(length(Data1{2,2}{1,2}),1);
    tempBed1 = zeros(length(Data1{2,2}{1,2}),1);
    for i = 1:length(Data1{2,2}{1,2})
        if Data1{2,2}{3,2}(end,i) >= floor(setTemp1)
            tempLine1(i,1) = find(Data1{2,2}{3,2}(:,i) >= floor(setTemp1),1);
        elseif Data1{2,2}{3,2}(end,i) < floor(setTemp1)
            tempLine1(i,1) = length(Data1{2,2}{3,2}(:,i));
        end
        tempFluid1(i,1) = Data1{2,2}{3,2}(tempLine1(i),i);
        tempBed1(i,1) = Data1{2,2}{4,2}(tempLine1(i)-1,i+1);  
    end

    if tempExit1 < floor(setTemp1)
        tempTarget1 = 'Not Achieved';
    elseif tempExit1 >= floor(setTemp1)
        tempTarget1 = (tempLine1(end)*Data1{2,2}{9,2}{5,2});
    end
    
    tempLocation1 = tempLine1.*pebbleDiameter{2,2}/conv.in2m;
end

if beds == 2
    fprintf('\nPost Processing...\n');
    % Final exit temp of fluid
    tempExit2 = Data2{2,2}{3,2}(end,end);
    setTemp2 = (Data2{2,2}{9,2}{10,2} + conv.f2r)*conv.r2k;
    % Location where fluid equals target temeprature
    tempLine2 = zeros(length(Data2{2,2}{1,2}),1);
    tempFluid2 = zeros(length(Data2{2,2}{1,2}),1);
    tempBed2 = zeros(length(Data2{2,2}{1,2}),1);
    for i = 1:(length(Data2{2,2}{1,2}))
        if Data2{2,2}{3,2}(end,i) >= floor(setTemp2)
            tempLine2(i,1) = find(Data2{2,2}{3,2}(:,i) >= floor(setTemp2),1);
        elseif Data2{2,2}{3,2}(end,i) < floor(setTemp2)
            tempLine2(i,1) = (length(Data2{2,2}{3,2}(:,i)));
        end
        if tempLine2(i,1) >= length(Data2{2,2}{4,2}(:,1))
            tempLine2(i,1) = length(Data2{2,2}{4,2}(:,1));
        end
        tempFluid2(i,1) = Data2{2,2}{3,2}(tempLine2(i)+1,i);
        tempBed2(i,1) = Data2{2,2}{4,2}(tempLine2(i,1),i+1);

    end

    if tempExit2 < floor(setTemp2)
        tempTarget2 = 'Not Achieved';
    elseif tempExit2 >= floor(setTemp2)
        tempTarget2 = (tempLine2(end)*Data2{2,2}{9,2}{5,2});
    end
    
    tempLocation2 = tempLine2.*pebbleDiameter{3,2}/conv.in2m;
end

write(tempBed1,tempFluid1,tempLocation1,Data1,1,ethane{3,2},ethane{4,2},solve{4,2});
write(tempBed2,tempFluid2,tempLocation2,Data2,2,ethane{3,2},ethane{4,2},solve{4,2});

%% Plotting
% Move post process, plotting and write to a dedicated post processing
% Bed 1
if settings{2,2}
    fprintf('\nPlotting results...\n');
    figure(1);
    grid on; hold on;
    title({'Bed 1 Ethane Temperature and Associated Bed Location vs. Time',...
            strcat(num2str(Data1{2,2}{9,2}{7,2}),'" ID',num2str(Data1{2,2}{9,2}{6,2}),'" Length',...
            num2str(Data1{2,2}{9,2}{5,2}),'" Pebbles')}); 
    xlabel('Time [s]'); ylabel('Temperature [K]');
    yyaxis left
    plot(Data1{2,2}{1,2},tempFluid1);
    yyaxis right
    set(gca,'ydir','reverse')
    plot((Data1{2,2}{1,2}),(tempLocation1));
    legend('Ethane Temp [K]','Bed Location [in]','location','southwest');
end

% Bed 2
if beds == 2
    fprintf('\nPlotting results...\n');
    figure(2);
    grid on; hold on;
    title({'Bed 2 Ethane Temperature and Associated Bed Location vs. Time',...
            strcat(num2str(Data2{2,2}{9,2}{7,2}),'" ID',num2str(Data2{2,2}{9,2}{6,2}),'" Length',...
            num2str(Data2{2,2}{9,2}{5,2}),'" Pebbles')}); 
    xlabel('Time [s]'); ylabel('Temperature [K]');
    yyaxis left
    plot(Data2{2,2}{1,2},tempFluid2);
    yyaxis right
    set(gca,'ydir','reverse')
    plot((Data2{2,2}{1,2}),(tempLocation2));
    legend('Ethane Temp [K]','Bed Location [in]','location','southwest');
end

