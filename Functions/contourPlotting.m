function contourPlotting(Data, timeIndices, bedNum)

% temperature
maxTemp = max(max(Data{2,2}{3,2}));
figure();
grid on; hold on;
title({'Normalized Ethane Temperature for Bed',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Distance Along Bed [in]'); 
ylabel('Normalized Ethane Temperature');
for i = 1:length(timeIndices)
    plot(Data{2,2}{2,2},Data{2,2}{3,2}(1:end-1,timeIndices(i,2))'./maxTemp);
end
legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
hold off

% density
maxDensity = max(max(Data{2,2}{10,2}));
figure();
grid on; hold on;
title({'Normalized Ethane Density for Bed',bedNum,strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Distance Along Bed [in]'); 
ylabel('Normalized Ethane Density');
for i = 1:length(timeIndices)
    plot(Data{2,2}{2,2},Data{2,2}{10,2}(:,timeIndices(i,2))'./maxDensity);
end
legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
hold off

% velocity
maxVelocity = max(max(Data{2,2}{6,2}));
figure();
grid on; hold on;
title({'Normalized Ethane Density for Bed',bedNum,strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Distance Along Bed [in]'); 
ylabel('Normalized Ethane Velocity');
for i = 1:length(timeIndices)
    plot(Data{2,2}{2,2},Data{2,2}{6,2}(:,timeIndices(i,2))'./maxVelocity);
end
legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
hold off

% pressure
maxPressure = max(max(Data{2,2}{5,2}));
figure();
grid on; hold on;
title({'Normalized Ethane Pressure for Bed',bedNum,strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Distance Along Bed [in]'); 
ylabel('Normalized Ethane Pressure');
for i = 1:length(timeIndices)
    plot(Data{2,2}{2,2},Data{2,2}{5,2}(1:end-1,timeIndices(i,2))'./maxPressure);
end
legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
hold off

% Reynolds number
maxReynolds = max(max(Data{2,2}{11,2}));
figure();
grid on; hold on;
title({'Normalized Ethane Reynolds Number for Bed',bedNum,strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Distance Along Bed [in]'); 
ylabel('Normalized Ethane Reynolds Number');
for i = 1:length(timeIndices)
    plot(Data{2,2}{2,2},Data{2,2}{11,2}(:,timeIndices(i,2))'./maxReynolds);
end
legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
hold off

% Prandtl number
maxPrandtl = max(max(Data{2,2}{12,2}));
figure();
grid on; hold on;
title({'Normalized Ethane Prandtl Number for Bed',bedNum,strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Distance Along Bed [in]'); 
ylabel('Normalized Ethane Prandtl Number');
for i = 1:length(timeIndices)
    plot(Data{2,2}{2,2},Data{2,2}{12,2}(:,timeIndices(i,2))'./maxPrandtl);
end
legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
hold off

% % mass flow ON HOLD UNTIL WE HAVE VARIED MASS FLOW
% maxMass = max(max(Data{2,2}{11,2}));
% figure();
% grid on; hold on;
% title({'Normalized Ethane Mass Flow Rate for Bed',bedNum,strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
%     num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
% xlabel('Distance Along Bed [in]'); 
% ylabel('Normalized Ethane Mass Flow Rate');
% for i = 1:length(timeIndices)
%     plot(Data{2,2}{2,2},Data{2,2}{11,2}(:,timeIndices(i,2))'./maxMass);
% end
% legend('t=0s','t=1s','t=2s','t=3s','t=4s','t=5s','location','southwest');
% hold off

end