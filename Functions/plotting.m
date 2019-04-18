function plotting(Data, tempFluid, tempLocation, bedNum)

figure();
grid on; hold on;
title({'Ethane Temperature and Associated Bed Location vs. Time for Bed ',bedNum,...
    strcat(num2str(Data{2,2}{9,2}{7,2}),'" ID ',num2str(Data{2,2}{9,2}{6,2}),'" Length ',...
    num2str(Data{2,2}{9,2}{5,2}),'" Pebbles')});
xlabel('Time [s]'); ylabel('Temperature [K]');
yyaxis left
plot(Data{2,2}{1,2},tempFluid);
yyaxis right
set(gca,'ydir','reverse')
plot((Data{2,2}{1,2}),(tempLocation));
legend('Ethane Temp [K]','Bed Location [in]','location','southwest');
hold off

end