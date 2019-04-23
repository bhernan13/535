function [tempBed, fluids, locations, tempTarget, timeIndices] = postprocess(Data,pebbleDiameter)

%% Load Conversions
conversions;

%% Post Process Results

% Temperature
% Final exit temp of fluid
tempExit = Data{2,2}{3,2}(end,end);
setTemp = (Data{2,2}{9,2}{10,2} + conv.f2r)*conv.r2k;
% Location where fluid equals target temeprature
tempLine = zeros(length(Data{2,2}{1,2}),1);
tempFluid = zeros(length(Data{2,2}{1,2}),1);
tempBed = zeros(length(Data{2,2}{1,2}),1);
for i = 1:length(Data{2,2}{1,2})
    if Data{2,2}{3,2}(end,i) >= floor(setTemp)
        tempLine(i,1) = find(Data{2,2}{3,2}(:,i) >= floor(setTemp),1);
    elseif Data{2,2}{3,2}(end,i) < floor(setTemp)
        tempLine(i,1) = length(Data{2,2}{3,2}(:,i));
    end
    if tempLine(i) == 1
        tempLine(i) = 2;
    end
    tempFluid(i,1) = Data{2,2}{3,2}(tempLine(i),i);
    tempBed(i,1) = Data{2,2}{4,2}(tempLine(i)-1,i+1);
end
if tempExit < floor(setTemp)
    tempTarget = 'Not Achieved';
elseif tempExit >= floor(setTemp)
    tempTarget = (tempLine(end)*Data{2,2}{9,2}{5,2});
end
tempLocation = tempLine.*pebbleDiameter/conv.in2m;

% Density gradient
gradDensity = findGradient(Data{2,2}{10,2});
[~,n] = size(gradDensity);
densLine = zeros(length(Data{2,2}{1,2}),1);
densFluid = zeros(length(Data{2,2}{1,2}),1);
for i = 1:n
    [densFluid(i,1),densLine(i,1)] = max(gradDensity(:,i));
end
densLocation = densLine.*pebbleDiameter/conv.in2m;
   
% Velocity gradient
gradVelocity = findGradient(Data{2,2}{6,2});
[~,n] = size(gradVelocity);
veloLine = zeros(length(Data{2,2}{1,2}),1);
veloFluid = zeros(length(Data{2,2}{1,2}),1);
for i = 1:n
    [veloFluid(i,1),veloLine(i,1)] = max(gradVelocity(:,i));
end
veloLocation = veloLine.*pebbleDiameter/conv.in2m;

% Pressure gradient
gradPressure = findGradient(Data{2,2}{5,2});
[~,n] = size(gradPressure);
presLine = zeros(length(Data{2,2}{1,2}),1);
presFluid = zeros(length(Data{2,2}{1,2}),1);
for i = 1:n
    [presFluid(i,1),presLine(i,1)] = max(gradPressure(:,i));
end
presLocation = presLine.*pebbleDiameter/conv.in2m;

% Reynolds number gradient
gradReynolds = findGradient(Data{2,2}{11,2});
[~,n] = size(gradReynolds);
reynLine = zeros(length(Data{2,2}{1,2}),1);
reynFluid = zeros(length(Data{2,2}{1,2}),1);
for i = 1:n
    [reynFluid(i,1),reynLine(i,1)] = max(gradReynolds(:,i));
end
reynLocation = reynLine.*pebbleDiameter/conv.in2m;

% Prandtl number gradient
gradPrandtl = findGradient(Data{2,2}{12,2});
[~,n] = size(gradPrandtl);
pranLine = zeros(length(Data{2,2}{1,2}),1);
pranFluid = zeros(length(Data{2,2}{1,2}),1);
for i = 1:n
    [pranFluid(i,1),pranLine(i,1)] = max(gradPrandtl(:,i));
end
pranLocation = pranLine.*pebbleDiameter/conv.in2m;


% generate indices for contour plots for 0-5 seconds
ind = 1;
timeIndices = zeros(6,2);
timeIndices(1,2) = ind;
timeIndices(:,1) = [0;1;2;3;4;5];
for i = 1:length(timeIndices)-1
    while Data{2,2}{1,2}(ind) < i
        ind = ind+1;
    end
    timeIndices(i+1,2) = ind;
end    

% collect location and fluid data
locations = cell(6,2);
locations{1,1} = 'Temperature [K]'; locations{2,1} = 'Density gradient [kg/m^3]';
locations{3,1} = 'Velocity gradient [m/s]'; locations{4,1} = 'Pressure gradient [Pa]';
locations{5,1} = 'Reynolds number gradient'; locations{6,1} = 'Prandtl number gradient';
locations{1,2} = tempLocation; locations{2,2} = densLocation;
locations{3,2} = veloLocation; locations{4,2} = presLocation;
locations{5,2} = reynLocation; locations{6,2} = pranLocation;
fluids = cell(6,2);
fluids{1,1} = 'Temperature [K]'; fluids{2,1} = 'Density gradient [kg/m^3]';
fluids{3,1} = 'Velocity gradient [m/s]'; fluids{4,1} = 'Pressure gradient [Pa]';
fluids{5,1} = 'Reynolds number gradient'; fluids{6,1} = 'Prandtl number gradient';
fluids{1,2} = tempFluid; fluids{2,2} = densFluid;
fluids{3,2} = veloFluid; fluids{4,2} = presFluid;
fluids{5,2} = reynFluid; fluids{6,2} = pranFluid;

end