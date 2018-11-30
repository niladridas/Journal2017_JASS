% main file
clc;close all;clear
fileID = fopen('25544_5.txt'); % TLE source
date_mee = tle_parse(fileID,'0');
samples = 50; rep = 1;
date_mee = date_mee(~any(isnan(date_mee),2),:); % First delete repeting rows
% Delete repreating dates
[C,ia,ic] = unique(date_mee(:,1:6),'rows','stable');
date_mee = date_mee(ia,:);

% for 25544_5.txt
date_mee = date_mee(1:150,:);

% date_mee = unique(date_mee,'rows','stable');
[x_est_OT,x_est_enkf,trace_cov] = ot_enkf(date_mee,samples,rep);

%%
% plot_ot_enkf(samples, rep, x_est_OT,x_est_enkf,date_mee,time_pt_days,'1')
plot(1:size(date_mee,1),(180/pi)*sqrt(trace_cov(:,6)));