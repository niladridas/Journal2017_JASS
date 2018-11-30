% plots for the journal
clc;clear; close all;
% load all data
% load('Nobs');
Nobs = 15;
load('real_states.mat');

ot50 = load('ot_50.mat');
ot60 = load('ot_60.mat');
ot70 = load('ot_70.mat');
ot80 = load('ot_80.mat');
ot90 = load('ot_90.mat');

ot_cov50 = load('ot_cov50.mat');
ot_cov60 = load('ot_cov60.mat');
ot_cov70 = load('ot_cov70.mat');
ot_cov80 = load('ot_cov80.mat');
ot_cov90 = load('ot_cov90.mat');



enkf50 = load('enkf_50.mat');
enkf60 = load('enkf_60.mat');
enkf70 = load('enkf_70.mat');
enkf80 = load('enkf_80.mat');
enkf90 = load('enkf_90.mat');


enkf_cov50 = load('enkf_cov50.mat');
enkf_cov60 = load('enkf_cov60.mat');
enkf_cov70 = load('enkf_cov70.mat');
enkf_cov80 = load('enkf_cov80.mat');
enkf_cov90 = load('enkf_cov90.mat');


ylabel_names = ['p','f','g','h','k','L'];


%%
for i = 1:6
    figure(i)
    subplot(2,3,1)
    a1 =errorbar(1:Nobs,ot50.x_est_OT(:,i),3*sqrt(ot_cov50.x_ot_cov(:,i)));
    a1.Color = 'red';
    a1.LineStyle = '-';
    a1.LineWidth = 1;
    a1.Marker = 'o';
    a1.MarkerSize = 3;
    a1.MarkerEdgeColor = 'red';
    a1.MarkerFaceColor = 'red';
    xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',13,'FontWeight','Bold');
    xticks(1:Nobs);
    xlim([1 15]);
    hold on
    c1 = plot(1:Nobs,real_states(:,i));
    c1.Color = 'blue';
    c1.LineStyle = '--';
    c1.LineWidth = 1;
    c1.Marker = 'o';
    c1.MarkerSize = 3;
    c1.MarkerEdgeColor = 'blue';
    c1.MarkerFaceColor = 'blue';
    title('Ensemble Size = 50');
    subplot(2,3,2)
    a2 =errorbar(1:Nobs,ot60.x_est_OT(:,i),3*sqrt(ot_cov60.x_ot_cov(:,i)));
    a2.Color = 'red';
    a2.LineStyle = '-';
    a2.LineWidth = 1;
    a2.Marker = 'o';
    a2.MarkerSize = 3;
    a2.MarkerEdgeColor = 'red';
    a2.MarkerFaceColor = 'red';
    xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',13,'FontWeight','Bold');
    xticks(1:Nobs);
    xlim([1 15]);
    hold on
    c2 = plot(1:Nobs,real_states(:,i));
    c2.Color = 'blue';
    c2.LineStyle = '--';
    c2.LineWidth = 1;
    c2.Marker = 'o';
    c2.MarkerSize = 3;
    c2.MarkerEdgeColor = 'blue';
    c2.MarkerFaceColor = 'blue';
    title('Ensemble Size = 60');
    subplot(2,3,3)
    a3 =errorbar(1:Nobs,ot70.x_est_OT(:,i),3*sqrt(ot_cov70.x_ot_cov(:,i)));
    a3.Color = 'red';
    a3.LineStyle = '-';
    a3.LineWidth = 1;
    a3.Marker = 'o';
    a3.MarkerSize = 3;
    a3.MarkerEdgeColor = 'red';
    a3.MarkerFaceColor = 'red';
    xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',13,'FontWeight','Bold');
    xticks(1:Nobs);
    xlim([1 15]);
    hold on
    c3 = plot(1:Nobs,real_states(:,i));
    c3.Color = 'blue';
    c3.LineStyle = '--';
    c3.LineWidth = 1;
    c3.Marker = 'o';
    c3.MarkerSize = 3;
    c3.MarkerEdgeColor = 'blue';
    c3.MarkerFaceColor = 'blue';    
    title('Ensemble Size = 70');
    hL = legend([a1,c1],'OT est. states','Real states');
    set(hL,'Position',[0.466 0.496 0.104 0.04]);

    
    subplot(2,3,4)
    a4 =errorbar(1:Nobs,enkf50.x_est_enkf(:,i),sqrt(enkf_cov50.x_enkf_cov(:,i)));
    a4.Color = [0 0.5 0];
    a4.LineStyle = '-';
    a4.LineWidth = 1;
    a4.Marker = 'o';
    a4.MarkerSize = 3;
    a4.MarkerEdgeColor = [0 0.5 0];
    a4.MarkerFaceColor = [0 0.5 0];
    xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',13,'FontWeight','Bold');
    xticks(1:Nobs);
    xlim([1 15]);
    hold on
    c4 = plot(1:Nobs,real_states(:,i));
    c4.Color = 'blue';
    c4.LineStyle = '--';
    c4.LineWidth = 1;
    c4.Marker = 'o';
    c4.MarkerSize = 3;
    c4.MarkerEdgeColor = 'blue';
    c4.MarkerFaceColor = 'blue';    
    subplot(2,3,5)
    a5 =errorbar(1:Nobs,enkf60.x_est_enkf(:,i),sqrt(enkf_cov60.x_enkf_cov(:,i)));
    a5.Color = [0 0.5 0];
    a5.LineStyle = '-';
    a5.LineWidth = 1;
    a5.Marker = 'o';
    a5.MarkerSize = 3;
    a5.MarkerEdgeColor = [0 0.5 0];
    a5.MarkerFaceColor = [0 0.5 0];
    xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',13,'FontWeight','Bold');
    xticks(1:Nobs);
    xlim([1 15]);
    hold on
    c5 = plot(1:Nobs,real_states(:,i));
    c5.Color = 'blue';
    c5.LineStyle = '--';
    c5.LineWidth = 1;
    c5.Marker = 'o';
    c5.MarkerSize = 3;
    c5.MarkerEdgeColor = 'blue';
    c5.MarkerFaceColor = 'blue';    
    subplot(2,3,6)
    a6 =errorbar(1:Nobs,enkf70.x_est_enkf(:,i),sqrt(enkf_cov70.x_enkf_cov(:,i)));
    a6.Color = [0 0.5 0];
    a6.LineStyle = '-';
    a6.LineWidth = 1;
    a6.Marker = 'o';
    a6.MarkerSize = 3;
    a6.MarkerEdgeColor = [0 0.5 0];
    a6.MarkerFaceColor = [0 0.5 0];
    xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',13,'FontWeight','Bold');
    xticks(1:Nobs);
    xlim([1 15]);
    hold on
    c6 = plot(1:Nobs,real_states(:,i));
    c6.Color = 'blue';
    c6.LineStyle = '--';
    c6.LineWidth = 1;
    c6.Marker = 'o';
    c6.MarkerSize = 3;
    c6.MarkerEdgeColor = 'blue';
    c6.MarkerFaceColor = 'blue';    
    hL = legend([a4,c4],'EnKF est. states','Real states');
    set(hL,'Position',[0.463 0.025 0.114 0.04]);

    filename = strcat(ylabel_names(i),'_un_plot.fig');
    savefig(filename);
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     fig.PaperOrientation = 'landscape';
%     fig.PaperUnits = 'normalized';
%     fig.PaperPosition = [0 0 1 1];
%     filename = strcat(ylabel_names(i),'plot');
%     print(filename,'-dpdf','-fillpage');
end


%% RMSE error  in ECI 
% 
% rmse_ot_50 = zeros(1,6);
% rmse_ot_60 = zeros(1,6);
% rmse_ot_70 = zeros(1,6);
% rmse_ot_80 = zeros(1,6);
% rmse_ot_90 = zeros(1,6);
% 
% rmse_enkf_50 = zeros(1,6);
% rmse_enkf_60 = zeros(1,6);
% rmse_enkf_70 = zeros(1,6);
% rmse_enkf_80 = zeros(1,6);
% rmse_enkf_90 = zeros(1,6);
% 
% 
% for i = 1:6
%     rmse_ot_50 = sqrt(mean((mee2eci(ot50.x_est_OT )-mee2eci(real_states)).^2));
%     rmse_ot_60 = sqrt(mean((mee2eci(ot60.x_est_OT )-mee2eci(real_states)).^2));
%     rmse_ot_70 = sqrt(mean((mee2eci(ot70.x_est_OT )-mee2eci(real_states)).^2));
%     rmse_ot_80 = sqrt(mean((mee2eci(ot80.x_est_OT )-mee2eci(real_states)).^2));
%     rmse_ot_90 = sqrt(mean((mee2eci(ot90.x_est_OT )-mee2eci(real_states)).^2));
%     rmse_ot = [rmse_ot_50;rmse_ot_60;rmse_ot_70;rmse_ot_80;rmse_ot_90];
% 
%     rmse_enkf_50 = sqrt(mean((mee2eci(enkf50.x_est_enkf )-mee2eci(real_states)).^2));
%     rmse_enkf_60 = sqrt(mean((mee2eci(enkf60.x_est_enkf )-mee2eci(real_states)).^2));
%     rmse_enkf_70 = sqrt(mean((mee2eci(enkf70.x_est_enkf )-mee2eci(real_states)).^2));
%     rmse_enkf_80 = sqrt(mean((mee2eci(enkf80.x_est_enkf )-mee2eci(real_states)).^2));
%     rmse_enkf_90 = sqrt(mean((mee2eci(enkf90.x_est_enkf )-mee2eci(real_states)).^2));
%     rmse_enkf = [rmse_enkf_50;rmse_enkf_60;rmse_enkf_70;rmse_enkf_80;rmse_enkf_90];
% end
% 
% 
% figure
% for i = 1:6
%     subplot(2,3,i)
%     stem(50:10:90,rmse_ot(:,i),'filled');
%     hold on
%     plot(50:10:90,rmse_ot(:,i));
% end




% function out = mee2eci(input)
% % input can be matrix with rows denoting one sample and column one variable
% out = zeros(size(input,1),6);
% for i = 1:size(input,1)
%     out(i,:) = coe2eci(mee2coe(input(i,:)));
% end
% end

