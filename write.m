function write(BedTemp, FluidTemp, LocTemp, Data, bed, p_condition,T_condition,m_dot)
global conv
    fprintf('\nWriting Results...\n');
    % Check if Results Directory exists, if not, make it.
    currentDir = pwd;
    if isfolder(strcat(currentDir,'\Results'))
    else
        mkdir('Results');
    end
    cd(strcat(currentDir,'\Results'));
    addpath(genpath(strcat(currentDir,'\Results')),'-end');
    
    % check if massflow is ramped or not, process accordingly
    if size(m_dot) > 1 % ramped
        mdot = 'ramp-';
    else % not ramped
        if isnumeric(m_dot)
            mdot = num2str(round(m_dot,3),'%0.3f');
        else
            mdot = 'ramp-';
        end
    end
    
        % Generate file name
    fileName = strcat(Data{2,2}{9,2}{3,2},'_',Data{2,2}{9,2}{4,2},'_',num2str(Data{2,2}{9,2}{6,2}),...
        '_',num2str(Data{2,2}{9,2}{7,2}),'_',num2str(Data{2,2}{9,2}{12,2}),'s_',...
        num2str(floor((p_condition/6894)*0.8),'%04.f'),'psi_',mdot,'lbm-s_',...
        num2str(floor((T_condition*1.8))),'R_',num2str(bed),'.txt');
        
    if exist(fileName,'file') == 2
        delete(fileName);
        fileID = fopen(fileName,'w+');
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\n','Time [s]','Pebble Temp [R]',...
            'Ethane Temp [R]','Ethane Temp=Target @ [in]','Exit Temperature [R]','Exit Pressure [psi]');
    else
        fileID = fopen(fileName,'w+');
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s \n','Time [s]','Pebble Temp [R]',...
            'Ethane Temp [R]','Ethane Temp=Target @ [in]','Exit Temperature [R]','Exit Pressure [psi]');
    end
    
    
    pressexit = Data{2,2}{5,2}(end,:)';
    tempexit = Data{2,2}{3,2}(end,:)';

    
    send = [Data{2,2}{1,2}(:), BedTemp/conv.r2k, FluidTemp/conv.r2k, LocTemp, tempexit/conv.r2k, pressexit/conv.psi2pa];
            
    for i = 1:size(send)
        fprintf(fileID,'%-12f\t',send(i,:));
        fprintf(fileID,'\n');
    end
    
    fclose(fileID);
    cd(currentDir);
end