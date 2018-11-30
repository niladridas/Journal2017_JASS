% Plot for consistency in the predictions
clc;clear;close all;
load OT_all_repeat
load EnKF_all_repeat

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
% Calculating mean of the predictions
% For OT
OT_mean = mean(OT_all_repeat,3);
% For EnKF
EnKF_mean = mean(EnKF_all_repeat,3);

% Calculating variance of the predictions
OT_std = std(OT_all_repeat,[],3);
EnKF_std = std(EnKF_all_repeat,[],3);

for i = 1:6
subplot(3,2,i)
errorbar(t_spot,OT_mean(:,i),OT_std(:,i),'r');
hold on
errorbar(t_spot,EnKF_mean(:,i),EnKF_std(:,i),'b');
legend('OT','EnKF')
end

