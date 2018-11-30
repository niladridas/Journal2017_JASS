clc;clear; close;
% Plots are generated seperately
load('./data/all_OT_data50.mat');
load('./data/all_EnKF_data50.mat');
load('./data/real_states.mat');
OT_est = mean(OT_all_repeat,3);
EnKF_est = mean(EnKF_all_repeat,3);
OT_sig = 3*mean(OT_allsig_repeat,3);
EnKF_sig = 3*mean(EnKF_allsig_repeat,3);
n_obs = size(real_states,1);

% ylimits = [6.71*10^6 6.85*10^6; -0.008 0.008; -0.008 0.006; 0.05 0.5; -0.5 -0.09; -4 8];
figure;
set(gcf,'position', [ 123    37   946   768]);

ylabels = {'${p_o}$(m)','${f_o}$','${g_o}$','${h_o}$','${k_o}$','${L_o}$(rad)'};
for i = 1:6
    ax1 = subplot(7,2,2*i-1);
%     text(0.5,0.9,char('a' + i - 1),'Units', 'Normalized', 'VerticalAlignment', 'Top')
    % Plot OT
    h1 = plot(1:n_obs,real_states(:,i),'b'); hold on;
    h1.LineWidth = 1;
    h2 = plot(1:n_obs,OT_est(:,i),'r'); hold on;
    h2.LineWidth = 1;
    h3 = plot(1:n_obs,OT_est(:,i)+OT_sig(:,i),'black--'); hold on;
    h3.LineWidth = 1; 
    h4 = plot(1:n_obs,OT_est(:,i)-OT_sig(:,i),'black--'); hold on;
    h4.LineWidth = 1;  
    ax = gca;
    ylabel(ylabels(i),'FontSize',18,'interpreter', 'latex','FontWeight','bold');
%     ylim(ylimits(i,:));
    xlim([1 n_obs]);
    ax.XTick = 1:2:n_obs;
%     ax.FontWeight = 'bold';
    ax.FontSize = 12;
    if i ==6
        xlabel('Time points (days)','FontSize',18,'interpreter', 'latex','FontWeight','bold');
%         ax.YTick = -4:4:8;
    end
    if i == 1
        title('OT (samples=50,repeat=10)','FontSize',18,'interpreter', 'latex','FontWeight','bold');
    end
    ax2 = subplot(7,2,2*i);
%     text(0.02,0.98,char('b' + i - 1),'Units', 'Normalized', 'VerticalAlignment', 'Top')
    % PLot EnKF
    h11 = plot(1:n_obs,real_states(:,i),'b'); hold on;
    h11.LineWidth = 1;
    h21 = plot(1:n_obs,EnKF_est(:,i),'r'); hold on;
    h21.LineWidth = 1;
    h31 = plot(1:n_obs,EnKF_est(:,i)+EnKF_sig(:,i),'black--'); hold on;
    h31.LineWidth = 1;
    h41 = plot(1:n_obs,EnKF_est(:,i)-EnKF_sig(:,i),'black--'); hold on;
    h41.LineWidth = 1;
    ax = gca;
        ylabel(ylabels(i),'FontSize',18,'interpreter', 'latex','FontWeight','bold');
%     ylim(ylimits(i,:));
    xlim([1 n_obs]);
    ax.XTick = 1:2:n_obs;
%     ax.FontWeight = 'bold';
    ax.FontSize = 12;
    if i ==6
        xlabel('Time points (days)','FontSize',18,'interpreter', 'latex','FontWeight','bold');
%         ax.YTick = -4:4:8;
    end
    if i == 1
        title('EnKF (samples=50,repeat=10)','FontSize',18,'interpreter', 'latex','FontWeight','bold');
    end
%     linkaxes([ax1,ax2],'xy')
end
hL = subplot(7,2,13.5);
poshL = get(hL,'position');     % Getting its position
lgd = legend(hL,[h1;h2;h4],'True','Estimate','$\pm$ 3$\sigma$');
set(lgd,'position',poshL,'Orientation','horizontal','FontSize',18,'interpreter', 'latex','FontWeight','bold');      % Adjusting legend's position
axis(hL,'off');

% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% box off

% get(gcf,'position')
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '-dpdf', './Plots/plotconsist50.pdf');