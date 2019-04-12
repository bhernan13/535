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
% same function as currently 

%% Post Processing
% Perform any analysis and reduction of data from the Data Matricie(s)
% output from the main code. Probably going to store any requested values
% output in a text file based on the output calls

if settings{1,2}
    % Bed 1
    fprintf('\nPost Processing...\n');
    [tempBed1, tempFluid1, tempLocation1, tempTarget1] = postprocess(Data1,pebbleDiameter{2,2});
    % Bed 2
    fprintf('\nPost Processing...\n');
    [tempBed2, tempFluid2, tempLocation2, tempTarget2] = postprocess(Data2,pebbleDiameter{3,2});
end

write(tempBed1,tempFluid1,tempLocation1,Data1,1,ethane{3,2},ethane{4,2},solve{4,2});
write(tempBed2,tempFluid2,tempLocation2,Data2,2,ethane{3,2},ethane{4,2},solve{4,2});

%% Plot/Contour
% Function(s) to handle plotting graphs and contours from data from above
% text files. 


