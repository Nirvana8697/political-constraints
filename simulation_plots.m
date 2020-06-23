% Simulation Plots for demeaned consumption, government spending, GDP,
% transfers and taxes and Spreads

% 1. Consumption, GDP and Government Spending. (They co-move)

% figure(1)
% plot(1:40, GDP_simul(40060:40099))
% hold on 
% plot(1:40, cons_simul(40060:40099))
% hold on 
% plot(1:40, gov_simul(40060:40099))

% GDP and spread simulation 

figure(1)
yyaxis left
plot(1:40, GDP_simul(10040:10079))
ylim([7,10.8])
hold on 
yyaxis right
plot(1:40, TB_simul(10040:10079))
ylim([-0.15,0.5])
yline(0)











