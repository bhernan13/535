function write(BedTemp, FluidTemp, LocTemp, Data, bed, condition)

    fprintf('\nWriting Results...\n');
    % Check if Results Directory exists, if not, make it.
    currentDir = pwd;
    if isfolder(strcat(currentDir,'\Results'))
    else
        mkdir('Results');
    end
    cd(strcat(currentDir,'\Results'));
    addpath(genpath(strcat(currentDir,'\Results')),'-end');

    % Generate file name
    fileName = strcat(Data{2,2}{9,2}{3,2},'_',Data{2,2}{9,2}{4,2},'_',num2str(Data{2,2}{9,2}{6,2}),'_',num2str(Data{2,2}{9,2}{7,2}),...
        '_',num2str(Data{2,2}{9,2}{12,2}),'s_',num2str(floor((condition/6894)*0.8)),'psi_',num2str(bed),'.txt');
        
    if exist(fileName,'file') == 2
        delete(fileName);
        fileID = fopen(fileName,'w+');
        fprintf(fileID,'%s\t%s\t%s\t%s\n','Time [s]','Pebble Temp [K]','Ethane Temp [K]',strcat(Data{2,2}{9,2}{3,2},'Temp=Target @ [in]') );
    else
        fileID = fopen(fileName,'w+');
        fprintf(fileID,'%s\t%s\t\t%s\t\t%s\n','Time [s]','Pebble Temp [K]','Ethane Temp [K]',strcat(Data{2,2}{9,2}{3,2},'Temp=Target @ [in]') );
    end
        
    send = [Data{2,2}{1,2}(:), BedTemp, FluidTemp, LocTemp];
            
    for i = 1:size(send)
        fprintf(fileID,'%-12f\t',send(i,:));
        fprintf(fileID,'\n');
    end
    
    fclose(fileID);
    cd(currentDir);
end