clc;clear;close;
% Plots are generated seperately
load('./data/all_OT_data50.mat');
% OT_all_repeat,OT_allsig_repeat, OT_time
load('./data/real_states.mat');
% First we generate the error of the estimate from the real state and
% calculate the covariance of the estimate across all the repititions (10)
% TO-DO: assuming n_repeat is same for fixed sample size
n_repeat = size(OT_time,2);
n_obs = size(real_states,1);
OTerror_estimate_repeat = OT_all_repeat - repmat(real_states,1,1,n_repeat);
OTerror_mean = mean(OTerror_estimate_repeat,3);
% 3D matrix containing error across n_repeat number of repititions
% calculate sigma along the 3D direction
OT_sig_repeat = zeros(size(real_states,1),size(real_states,2));
for i = 1:size(real_states,1)
    for j = 1:size(real_states,2)
        OT_sig_repeat(i,j) = cov(reshape(OTerror_estimate_repeat(i,j,:),n_repeat,1));
        OT_sig_repeat(i,j) = sqrt(OT_sig_repeat(i,j));
    end
end
%%
load('./data/all_EnKF_data50.mat');
% EnKF_all_repeat,EnKF_allsig_repeat,EnKF_time
n_repeat = size(EnKF_time,2);
EnKFerror_estimate_repeat = EnKF_all_repeat - repmat(real_states,1,1,n_repeat);
EnKFerror_mean  = mean(EnKFerror_estimate_repeat,3);
% 3D matrix containing error across n_repeat number of repititions
% calculate sigma along the 3D direction
EnKF_sig_repeat = zeros(size(real_states,1),size(real_states,2));
for i = 1:size(real_states,1)
    for j = 1:size(real_states,2)
        EnKF_sig_repeat(i,j) = cov(reshape(EnKFerror_estimate_repeat(i,j,:),n_repeat,1));
        EnKF_sig_repeat(i,j) = sqrt(EnKF_sig_repeat(i,j));
    end
end
%%
figure;
ylabels = {'$e_{p_o}$(m)','$e_{f_o}$','$e_{g_o}$','$e_{h_o}$','$e_{k_o}$','$e_{L_o}$(rad)'};
for i = 1:6
    ax1 = subplot(6,2,2*i-1);
    % Plot OT
    h1 = errorbar(1:n_obs,OTerror_mean(:,i),OT_sig_repeat(:,i)); hold on;
    h1.LineWidth = 1;
    ax = gca;
%     ylabel(ylabels(i),'FontSize',14,'FontWeight','bold','interpreter', 'latex');
    ylabel(ylabels(i),'FontSize',14,'interpreter', 'latex');
    xlim([1 n_obs]);
    ax.XTick = 1:2:n_obs;
%     ax.FontWeight = 'bold';
    ax.FontSize = 8;
    if i ==6
        xlabel('Time points (days)');
    end
    if i == 1
        title('OT (samples=50,repeat=10)');
    end
    ax2 = subplot(6,2,2*i);
    % PLot EnKF
    h2 = errorbar(1:n_obs,EnKFerror_mean(:,i),EnKF_sig_repeat(:,i)); hold on;
    h2.LineWidth = 1;
    ax = gca;
%         ylabel(ylabels(i),'FontSize',14,'FontWeight','bold','interpreter', 'latex');
    ylabel(ylabels(i),'FontSize',14,'interpreter', 'latex');
    xlim([1 n_obs]);
    ax.XTick = 1:2:n_obs;
%     ax.FontWeight = 'bold';
    ax.FontSize = 8;
    if i ==6
        xlabel('Time points (days)');
    end
    if i == 1
        title('EnKF (samples=50,repeat=10)');
    end
    linkaxes([ax1,ax2],'y')
end
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '-dpdf', './Plots/plotrobust50.pdf');
