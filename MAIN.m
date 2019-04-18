%% Thermal Analysis of Pebble Bed Heaters
% Description and info here

%% Add Paths to Functions and Refprop
% Refprop folder must be stored on the 'C:\temp' filepath. 
% ECN drives prevent running refprop from within them. 
cd('C:\Users\gymna\Documents\Purdue\AAE 535\github\535-mmelcer13-postprocess');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentDir = pwd;
refpropDir = 'C:\temp\REFPROP';
if isfolder(refpropDir)
else 
    error('REFPROP folder is not in correct location. \n Copy REFPROP folder to %s \n','C:\temp or change refpropDir');
end 
addpath(genpath(refpropDir),'-end');
addpath(genpath(currentDir),'-end');
addpath(genpath(strcat(currentDir,'\Results')));

%% General Inputs
% Inputs for target conditions, mass flow and code wide parameters

%% Bed Specific Inputs
% Inputs for Bed #X 

%% Pebble Material Inputs
% Nice to be able to load this orderly and pull from this once the bed
% properties are input. Almost like a library type thing. Plot properties
% should be a post processing option

%% Solve Settings
% Settings for solution/post processing that are input options. dt will be
% moved to calculate in the main code

%% Output Calls
% Expand current method and make it cleaner for pulling any number of
% outputs based on what's needed instead of getting a slog of data

%% Data Matricie(s)
% Setup Matricie(s) for storing performance values for vector input runs

%% Main Code Call
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

%% Post Processing
% Perform any analysis and reduction of data from the Data Matricie(s)
% output from the main code. Probably going to store any requested values
% output in a text file based on the output calls

if settings{1,2}
    % Bed 1
    fprintf('\nPost Processing...\n');
    [tempBed1, tempFluid1, tempLocation1, tempTarget1,timeIndices1] = postprocess(Data1,...
        pebbleDiameter{2,2});
    % Bed 2
    fprintf('\nPost Processing...\n');
    [tempBed2, tempFluid2, tempLocation2, tempTarget2,timeIndices2] = postprocess(Data2,...
        pebbleDiameter{3,2});
end
write(tempBed1,tempFluid1,tempLocation1,Data1,1,ethane{3,2},ethane{4,2},solve{4,2});
write(tempBed2,tempFluid2,tempLocation2,Data2,2,ethane{3,2},ethane{4,2},solve{4,2});

%% Plot/Contour
% Function(s) to handle plotting graphs and contours from data from above
% text files. 
if settings{2,2}
    fprintf('\nPlotting results...\n');
    plotting(Data1, tempFluid1, tempLocation1, 1);
    plotting(Data2, tempFluid2, tempLocation2, 2);
end

if settings{5,2}
    fprintf('\nPlotting contours...\n');
    contourPlotting(Data1, timeIndices1, 1);
    contourPlotting(Data2, timeIndices2, 2);
end

