%% Post Processing of Thermal Design and Performance Analysis for Pebble Bed Heaters

clear all
close all

% Code author: Madie Melcer
% Created: March 17, 2019

resultsDir = 'Z:\2019\HPDS\Thermal_Fluids\PB_MM\Results';

addpath(genpath(resultsDir),'-end');

%% Generate file names
% for testing:
fluid = 'ethane';
material = 'SS304';
bedlength = [18 50]; % in - Bed 2 is 18 inches, Bed 1 is 50 inches
diameter = 3.4; % in
time = 5; % s
temp = [810,865,920,975]; % R
pres = [1000,1600,2200]; % psi
pres_in = [2000,2750]; % psi
mdot = [0.358,0.500]; % lbm/s
% pres = [400,1000,1600,2200]; % psi
% pres_in = [500,1250,2000,2750]; % psi
% mdot = [0.075,0.217,0.358,0.500]; % lbm/s

% for user input:
% fprintf('Entering pebble bed parameters: \nIf you want to enter multiple ')
% fprintf('numerical inputs, enter as a matrix: [ 1 2 3 ]. \nIf you want to enter ')
% fprintf('multiple character or mixed inputs, enter as an vector separated by semicolons: [ ''one'';')
% fprintf('''two2'';''3three'' ].\n\n')
% fluid = char(input('Input fluid: ','s'));
% material = input('Enter pebble material(s) - SS304,SS316,CU201,AL2017,BRASS260: ');
% bedlength = input('Enter desired pebble bed length(s) in inches: ');
% diameter = input('Enter diameter of pebble bed in inches: ');

a = length(temp);
b = length(bedlength);
c = length(pres);
d = length(mdot);

% generate file names
for p = 1:b
    if bedlength(p) == 18
        bednum = 2;
    elseif bedlength(p) == 50
        bednum = 1;
    else
        fprintf('unexpected bed length\n');
    end
    
    for q = 1:c
        for o = 1:a
            for r = 1:d
                fileName(p,:,q,o,r) = strcat(fluid,'_',material,'_',num2str(bedlength(p)),...
                    '_',num2str(diameter),'_',num2str(time),'s_',num2str(pres(q),'%04.f'),'psi_',...
                    num2str(mdot(r),'%0.3f'),'lbm-s_',num2str(temp(o)),'R_',num2str(bednum),'.txt'); 
            end
        end
    end
end

%% Read data from results files

for p = 1:b
    for q = 1:c
        for o = 1:a
            for r = 1:d
                file1 = fileName(p,:,q,o,r);
                fileID(p,q,o,r) = fopen(file1,'r');
                header(p,:,q,o,r) = fgetl(fileID(p,q,o,r));
                rawdata(:,p,q,o,r) = fscanf(fileID(p,q,o,r),'%f')';
                fclose(fileID(p,q,o,r));
            end
        end
    end
end

%% extract column titles from .txt file header

% assumes all files have same variables/column titles

[~,L,~,~,~] = size(header);
col_count = 1;
columntitle(col_count,1) = header(1,1);
i = 2;
col_index = 2;
while i <= L
    columntitle(col_count,col_index) = header(1,i,1,1,1);
    col_index = col_index + 1;
    if header(1,i,1,1,1) == ']' % moves to next column title
        col_count = col_count + 1;
        col_index = 1;
        i = i + 1; % skips tab spacer
    end
    i = i + 1;
end
col_count = col_count - 1; % number of columns in data file  

%% Sort data from .txt file into appropriate columns

n = length(rawdata(:,1,1,1,1))/col_count; % number of data points
for p = 1:b
    for q = 1:c
        for o = 1:a
            for r = 1:d
                for i = 1:n
                    data(i,:,p,q,o,r) = rawdata(((i-1)*col_count)+1:...
                        ((i-1)*col_count)+col_count,p,q,o,r);
                end
            end
        end
    end
end

%% Plot data

% Bed 1
for j = 1:col_count-1
    figure(j)
    hold on
    count = 1;
    for q = 1:c % pressure
        for o = 1:a % temp
            for r = 1:d % massflow
                plot(data(:,1,b,q,o,r),data(:,j+1,b,q,o,r),'linewidth',1.5)
                leg(count,:) = strcat(num2str(pres(q)),'psi, ',num2str(temp(o)),...
                    'R, ',num2str(mdot(r),'%0.3f'),'lbm/s');
                count = count + 1;
            end
        end
    end
    xlabel(columntitle(1,:))
    ylabel(columntitle(j+1,:))
    title('Bed 1')
    legend(leg,'location','south')
    hold off
end

% Bed 2
for j = 1:col_count-1
    figure(j+col_count-1)
    hold on
    count = 1;
    for q = 1:c % pressure
        for o = 1:a % temp
            for r = 1:d % massflow
                plot(data(:,1,1,q,o,r),data(:,j+1,1,q,o,r),'linewidth',1.5)
                leg(count,:) = strcat(num2str(pres(q)),'psi, ',num2str(temp(o)),...
                    'R, ',num2str(mdot(r),'%0.3f'),'lbm/s');
                count = count + 1;
            end
        end
    end
    xlabel(columntitle(1,:))
    ylabel(columntitle(j+1,:))
    title('Bed 2')
    legend(leg,'location','south')
    hold off
end

fprintf('\nDone general data\n')


%% nozzle-specific runs - requires hard coding
% data file names
nozzle1 = 'RS-68Aethane_SS304_18_3.4_5s_1350psi_0.500lbm-s_752R_2.txt';
nozzle2 = 'RS-68Aethane_SS304_18_3.4_5s_1350psi_0.358lbm-s_752R_2.txt';
nozzle3 = 'RS-25ethane_SS304_18_3.4_5s_2200psi_0.500lbm-s_637R_2.txt';
nozzle4 = 'RS-25ethane_SS304_18_3.4_5s_2200psi_0.358lbm-s_637R_2.txt';
nozzle5 = 'RS-68Aethane_SS304_50_3.4_5s_1350psi_0.500lbm-s_752R_1.txt';
nozzle6 = 'RS-68Aethane_SS304_50_3.4_5s_1350psi_0.358lbm-s_752R_1.txt';
nozzle7 = 'RS-25ethane_SS304_50_3.4_5s_2200psi_0.500lbm-s_637R_1.txt';
nozzle8 = 'RS-25ethane_SS304_50_3.4_5s_2200psi_0.358lbm-s_637R_1.txt';

% extract data from files
fileID = fopen(nozzle1,'r'); fgetl(fileID);
rawdata_noz(:,1) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle2,'r'); fgetl(fileID);
rawdata_noz(:,2) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle3,'r'); fgetl(fileID);
rawdata_noz(:,3) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle4,'r'); fgetl(fileID);
rawdata_noz(:,4) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle5,'r'); fgetl(fileID);
rawdata_noz(:,5) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle6,'r'); fgetl(fileID);
rawdata_noz(:,6) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle7,'r'); fgetl(fileID);
rawdata_noz(:,7) = fscanf(fileID,'%f')'; fclose(fileID);
fileID = fopen(nozzle8,'r'); fgetl(fileID);
rawdata_noz(:,8) = fscanf(fileID,'%f')'; fclose(fileID);

% headers are all still the same, so skip directly to sorting data
for m = 1:8
    for i = 1:n
        data_noz(i,:,m) = rawdata_noz(((i-1)*col_count)+1:...
            ((i-1)*col_count)+col_count,m)';
    end
end

% plot data for Bed 1
for j = 1:col_count-1
    figure(j+(col_count-1)*2)
    hold on
    count = 1;
    for m = 5:8
        plot(data_noz(:,1,m),data_noz(:,j+1,m),'linewidth',1.5)
    end
    xlabel(columntitle(1,:))
    ylabel(columntitle(j+1,:))
    title('Nozzle specific data - Bed 1')
    legend('RS-68A 0.5lbm/s','RS-68A 0.358lbm/s','RS-25 0.5lbm/s',...
        'RS-25 0.358lbm/s','location','south')
    hold off
end

% plot data for Bed 2
for j = 1:col_count-1
    figure(j+(col_count-1)*2+5)
    hold on
    count = 1;
    for m = 1:4
        plot(data_noz(:,1,m),data_noz(:,j+1,m),'linewidth',1.5)
    end
    xlabel(columntitle(1,:))
    ylabel(columntitle(j+1,:))
    title('Nozzle specific data - Bed 2')
    legend('RS-68A 0.5lbm/s','RS-68A 0.358lbm/s','RS-25 0.5lbm/s',...
        'RS-25 0.358lbm/s','location','south')
    hold off
end

fprintf('\nDone nozzle data\n')