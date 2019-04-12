function [Data] = readText(fileName, numColumns, headerLines)
% Read text files. Columns are outputs in each cell of the cell array
% "Data". For a delimiter other than \t, add 'Delimiter','character' to
% textscan, where character is the delimiter of choice.

fileID = fopen(fileName,'r');
numRows = Inf;
format = '%s';
for i = 1:numColumns-1
    format = strcat(format,' %s');
end

% format = strcat(format,'\n');
Data = textscan(fileID,format,numRows,'Headerlines',headerLines);

fclose(fileID);
end