function [tempBed, tempFluid, tempLocation, tempTarget] = postprocess(Data,pebbleDiameter)

%% Load Conversions
conversions;

%% Post Process Results

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
end