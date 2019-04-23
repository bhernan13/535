function plotting(Data, fluids, locations, bedNum)

figure();
grid on; hold on;
title({'Ethane Temperature and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Temperature [K]');
yyaxis left
plot(Data{2,2}{1,2},fluids{1,2});
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(locations{1,2}));
legend('Ethane Temp [K]','Bed Location [in]','location','southwest');
hold off

figure();
grid on; hold on;
title({'Ethane Density Gradient and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Density Gradient [kg/m^3]');
yyaxis left
plot(Data{2,2}{1,2},fluids{2,2});
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(locations{1,2}));
legend('Ethane Density Gradient[kg/m^3]','Bed Location [in]','location','southwest');
hold off

figure();
grid on; hold on;
title({'Ethane Velocity Gradient and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Velocity Gradient [m/s]');
yyaxis left
plot(Data{2,2}{1,2},fluids{3,2});
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(locations{1,2}));
legend('Ethane Velocity Gradient [m/s]','Bed Location [in]','location','southwest');
hold off

figure();
grid on; hold on;
title({'Ethane Pressure Gradient and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Pressure Gradient [Pa]');
yyaxis left
plot(Data{2,2}{1,2},fluids{4,2});
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(locations{1,2}));
legend('Ethane Pressure Gradient [Pa]','Bed Location [in]','location','southwest');
hold off

figure();
grid on; hold on;
title({'Ethane Reynolds Number Gradient and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Reynolds Number Gradient');
yyaxis left
plot(Data{2,2}{1,2},fluids{5,2});
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(locations{1,2}));
legend('Ethane Reynolds Number Gradient','Bed Location [in]','location','southwest');
hold off

figure();
grid on; hold on;
title({'Ethane Prandtl Number Gradient and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Prandtl Number Gradient');
yyaxis left
plot(Data{2,2}{1,2},fluids{6,2});
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(locations{1,2}));
legend('Ethane Prandtl Number Gradient','Bed Location [in]','location','southwest');
hold off

end