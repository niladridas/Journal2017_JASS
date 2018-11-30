% plot rmse time and samples
close all;clc;clear;
% Load observations
load date_mee.mat;
load t1;
t_spot = zeros(1,length(t1));
for i = 1:length(t1)
    t_spot(i) = sum(t1(1:i));
end
% t_spot is in seconds
% conversting it into days
t_spot = t_spot./(60*60*24);
load syn_state
load all_OT_data


temp_ob = abs(syn_state(:,1:5));
norm_x =  max(temp_ob,[],1);

for i = 1:5
x_est_OT(:,i) = 100.*x_est_OT(:,i)./norm_x(i);
syn_state(:,i) = 100.*syn_state(:,i)./norm_x(i);
end


for i = 1:6
subplot(3,2,i)
plot(t_spot,x_est_OT(:,i),'r-+');
xlabel('Time in Days')
hold on
plot(t_spot,syn_state(:,i),'g-*');
xlabel('Time in Days')
legend('OT','Actual')
grid on
end
% 
% 
% 
% figure(2)
% plot(30:10:80,RMSE_ot(:,7))
% hold on
% plot(30:10:80,RMSE_enkf(:,7))
% legend('OT','EnKF')
% xlabel('Sample Size')
% ylabel('Run Time')
