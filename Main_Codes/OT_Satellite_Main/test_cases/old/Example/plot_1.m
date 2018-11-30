% plot rmse time and samples
load RMSE_ot
load RMSE_enkf

% % axis for m/s
% b=axes('Position',[.1 .1 .8 1e-12]);
% set(b,'Units','normalized');
% set(b,'Color','none');
% 
% % axis for km/h with stem-plot
% a=axes('Position',[.1 .2 .8 .7]);
% % set(a,'Units','normalized');
% % stem(a,M(:,1).*3.6, M(:,3));
% 
% % set limits and labels
% % set(a,'xlim',[0 xmaxa]);
% % set(b,'xlim',[0 xmaxb]);
% xlabel(a,'Speed (km/h)')
% xlabel(b,'Speed (m/s)')
% ylabel(a,'Samples');
% title(a,'Double x-axis plot');

% 
% for i = 1:6
% subplot(3,2,i)
% plot(30:10:80,RMSE_ot(:,i),'r')
% hold on
% plot(30:10:80,RMSE_enkf(:,i),'b')
% grid on
% legend('OT','EnKF')
% % a=axes('Position',[.1 .2 .8 .7]);
% % b=axes('Position',[.1 .1 .8 1e-12]);
% % xlabel(a,'Speed (km/h)');
% % xlabel(b,'Speed (m/s)');
% % ylabel(a,'Samples');
% end



figure(2)
plot(30:10:80,RMSE_ot(:,7))
hold on
plot(30:10:80,RMSE_enkf(:,7))
legend('OT','EnKF')
xlabel('Sample Size')
ylabel('Run Time')
