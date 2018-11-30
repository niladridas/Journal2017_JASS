% General Plot function
clc;close all;clear
load('date_pos_vel.mat');
load('date_mee.mat');
%     load('/Users/niladridas/Documents/GITHUB/Kalman/Optimal_Transport/Satellite_Data_Analysis/OT_Satellite_Main/data/mee_ob.mat');
%     load('/Users/niladridas/Documents/GITHUB/Kalman/Optimal_Transport/Satellite_Data_Analysis/OT_Satellite_Main/data/mee_predict.mat');
load('time_pt_days.mat');
load('enkf_repeat_15.mat');
load('ot_repeat_15.mat');


% %% Plot observed distance from the center of the earth
% x_pos = date_pos_vel(:,7);
% y_pos = date_pos_vel(:,8);
% z_pos = date_pos_vel(:,9);
% figure(1)
% temp1  = sqrt(sum([x_pos,y_pos,z_pos].*[x_pos,y_pos,z_pos],2));
% scatter(1:size(x_pos,1),temp1,'filled')
% hold on
% plot(1:size(x_pos,1),temp1)
% xlabel('Observation Count');
% ylabel('Distance in Km from earth center')
% 
% %%
% error_e = mee_ob - mee_predict;
% figure(2)
% for i = 1:6
% subplot(2,3,i)
% scatter(time_pt_days,error_e(:,i))
% end

%% 
% samples = 40;
% t_spot = time_pt_days;
% ylabel_names = ['p','f','g','h','k','L'];
% measured_output = date_mee(:,7:12);
% figure(3)
% for i = 1:6
%     subplot(3,2,i);
%     a = plot(t_spot,x_est_OT(:,i));
%     a.Color = 'red';
%     a.LineStyle = '-';
%     a.LineWidth = 1;
%     a.Marker = 'o';
%     a.MarkerSize = 1.5;
%     a.MarkerEdgeColor = 'red';
%     a.MarkerFaceColor = 'red';
%     xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
%     ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% %     set(gca,'fontWeight','bold','Xtick',0:15);
%     hold on
%     c = plot(t_spot,measured_output(:,i));
%     c.Color = [0 0.5 0];
%     c.LineStyle = '-';
%     c.LineWidth = 1;
%     c.Marker = 'o';
%     c.MarkerSize = 1.5;
%     c.MarkerEdgeColor = [0 0.5 0];
%     hold on
% %     grid on
%     d = plot(t_spot,x_est_enkf(:,i));
%     d.Color = 'blue';
%     d.LineStyle = '-';
%     d.LineWidth = 1;
%     d.Marker = 'o';
%     d.MarkerSize = 1.5;
%     d.MarkerEdgeColor = 'blue';
%     d.MarkerFaceColor = 'blue';
%     legend('OT','Obs','EnFK')     
% end
% temp_datatype = 'True States';  
% title_str = strcat('Plot of OT and EnKF predictions for',{' '},num2str(samples),' samples');
% axes( 'Position', [0, 0.95, 1, 0.05] ) ;
% set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
% text( 0.55, 0, title_str, 'FontSize', 10', 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
% file_name = strcat('samplesOT',num2str(samples));
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(file_name,'-dpdf','-fillpage');
%%
% samples = 40;
% t_spot = time_pt_days;
% ylabel_names = ['p','f','g','h','k','L'];
% measured_output = date_mee(:,7:12);
% figure(4)
% for i = 1:6
%     subplot(3,2,i);
%     a = plot(t_spot,x_est_enkf(:,i));
%     a.Color = 'red';
%     a.LineStyle = '-';
%     a.LineWidth = 1;
%     a.Marker = 'o';
%     a.MarkerSize = 2;
%     a.MarkerEdgeColor = 'red';
%     a.MarkerFaceColor = 'red';
%     xlabel('Time (Days)','FontName','sans-serif','FontSize',20,'FontWeight','Bold');
%     ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',20,'FontWeight','Bold');
% %     set(gca,'fontWeight','bold','Xtick',0:15);
%     hold on
%     c = plot(t_spot,measured_output(:,i));
%     c.Color = [0 0.5 0];
%     c.LineStyle = '-';
%     c.LineWidth = 1;
%     c.Marker = 'o';
%     c.MarkerSize = 2;
%     c.MarkerEdgeColor = [0 0.5 0];
%     legend('OT','True')     
% %     grid on
% end
% temp_datatype = 'True States';  
% title_str = strcat('Plot of EnKF predictions for',{' '},num2str(samples),' samples');
% axes( 'Position', [0, 0.95, 1, 0.05] ) ;
% set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
% text( 0.55, 0, title_str, 'FontSize', 20', 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
% file_name = strcat('samplesEnKF',num2str(samples));
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(file_name,'-dpdf','-fillpage');
%%
%error plots
x_enkf_mean = mean(x_est_enkf,3);
x_enkf_var =  sqrt(var(x_est_enkf,[],3));
x_ot_mean = mean(x_est_OT,3);
x_ot_var =  sqrt(var(x_est_OT,[],3));
%
samples = 40;
t_spot = time_pt_days;
ylabel_names = ['p','f','g','h','k','L'];
measured_output = date_mee(:,7:12);
figure(3)
for i = 1:6
    subplot(3,2,i);
%     a = plot(t_spot,x_est_OT(:,i));
      a =errorbar(t_spot,x_ot_mean(:,i),x_ot_var(:,i));
      a.Color = 'red';
        a.LineStyle = '-';
        a.LineWidth = 2;
        a.Marker = 'o';
        a.MarkerSize = 4;
        a.MarkerEdgeColor = 'red';
        a.MarkerFaceColor = 'red';
    xlabel('Time (Days)','FontName','sans-serif','FontSize',20,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',20,'FontWeight','Bold');
%     set(gca,'fontWeight','bold','Xtick',0:15);
    hold on
    c = plot(t_spot,measured_output(:,i));
    c.Color = [0 0.5 0];
    c.LineStyle = '-';
    c.LineWidth = 2;
    c.Marker = 'o';
    c.MarkerSize = 4;
    c.MarkerEdgeColor = [0 0.5 0];
    legend('OT','True')     
    grid on
end
temp_datatype = 'True States';  
title_str = strcat('Plot of OT predictions for',{' '},num2str(samples),' samples with ',{' '}, temp_datatype);
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.55, 0, title_str, 'FontSize', 20', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
  

figure(4)
for i = 1:6
    subplot(3,2,i);
%     a = plot(t_spot,x_est_OT(:,i));
      a =errorbar(t_spot,x_enkf_mean(:,i),x_enkf_var(:,i));
      a.Color = 'red';
        a.LineStyle = '-';
        a.LineWidth = 2;
        a.Marker = 'o';
        a.MarkerSize = 4;
        a.MarkerEdgeColor = 'red';
        a.MarkerFaceColor = 'red';
    xlabel('Time (Days)','FontName','sans-serif','FontSize',20,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',20,'FontWeight','Bold');
%     set(gca,'fontWeight','bold','Xtick',0:15);
    hold on
    c = plot(t_spot,measured_output(:,i));
    c.Color = [0 0.5 0];
    c.LineStyle = '-';
    c.LineWidth = 2;
    c.Marker = 'o';
    c.MarkerSize = 4;
    c.MarkerEdgeColor = [0 0.5 0];
    legend('EnKF','True')     
    grid on
end
temp_datatype = 'True States';  
title_str = strcat('Plot of EnKF predictions for',{' '},num2str(samples),' samples with ',{' '}, temp_datatype);
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.55, 0, title_str, 'FontSize', 20', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
  













