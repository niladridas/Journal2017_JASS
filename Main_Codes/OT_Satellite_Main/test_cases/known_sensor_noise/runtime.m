% Run time plots
clc;clear; 
load('./data/all_OT_data50.mat');
load('./data/all_EnKF_data50.mat');
% Collect run time
run50OT = mean(OT_time,2);
run50EnKF = mean(EnKF_time,2);
load('./data/all_OT_data60.mat');
load('./data/all_EnKF_data60.mat');
% Collect run time
run60OT = mean(OT_time,2);
run60EnKF = mean(EnKF_time,2);
load('./data/all_OT_data70.mat');
load('./data/all_EnKF_data70.mat');
% Collect run time
run70OT = mean(OT_time,2);
run70EnKF = mean(EnKF_time,2);
% % Plotting
% p1 = plot(1:15,run50OT./run50EnKF,'r');hold on;
% p2 = plot(1:15,run60OT./run60EnKF,'b');hold on;
% p3 = plot(1:15,run70OT./run70EnKF,'black');
% p1.LineWidth = 1;
% p2.LineWidth = 1;
% p3.LineWidth = 1;
% legend([p1;p2;p3],'50 samples','60 samples','70 samples');
% ax = gca;
% xlim([1 15]);
% ax.XTick = 1:1:15;
% xlabel('Observations');
% ylabel('OT time / EnKF time ')
% ax.FontWeight = 'bold';
% ax.FontSize = 8;
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, '-dpdf', './Plots/plottimequarterday.pdf');