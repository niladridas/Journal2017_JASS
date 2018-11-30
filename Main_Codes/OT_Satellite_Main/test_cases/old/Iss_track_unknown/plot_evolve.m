clc;clear;close all;

ot50all = load('ot_all50.mat');
ot60all = load('ot_all60.mat');
ot70all = load('ot_all70.mat');
ot80all = load('ot_all80.mat');
ot90all = load('ot_all90.mat');

enkf50all = load('enkf_all50.mat');
enkf60all = load('enkf_all60.mat');
enkf70all = load('enkf_all70.mat');
enkf80all = load('enkf_all80.mat');
enkf90all = load('enkf_all90.mat');

% choose a time point 
tp = 9;

ot50 = mee2eci(ot50all.x_otall(:,:,tp)')';
ot60 = mee2eci(ot60all.x_otall(:,:,tp)')';
ot70 = mee2eci(ot70all.x_otall(:,:,tp)')';
ot80 = mee2eci(ot80all.x_otall(:,:,tp)')';
ot90 = mee2eci(ot90all.x_otall(:,:,tp)')';


enkf50 = mee2eci(enkf50all.x_enkfall(:,:,tp)')';
enkf60 = mee2eci(enkf60all.x_enkfall(:,:,tp)')';
enkf70 = mee2eci(enkf70all.x_enkfall(:,:,tp)')';
enkf80 = mee2eci(enkf80all.x_enkfall(:,:,tp)')';
enkf90 = mee2eci(enkf90all.x_enkfall(:,:,tp)')';


figure
xlab ={'X_{pos}','Y_{pos}','Z_{pos}'};
for i = 1:3
    subplot(3,2,2*i-1)
    [f1,x1i] = ksdensity(ot50(i,:));
    area(x1i,f1,'Facecolor',[1 0.5 0.5]);
    xlabel(xlab(i));
    ylabel(strcat('Pr'),'rot',0);
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [0 .2 0]);
    subplot(3,2,2*i)
    [f2,x2i] = ksdensity(enkf50(i,:));
    area(x2i,f2,'Facecolor',[0.5 0.7 0.5]);
    xlabel(xlab(i));
    ylabel(strcat('Pr'),'rot',0);
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [0 .1 0]);
end

xlab = {'X_{vel}','Y_{vel}','Z_{vel}'};

figure
for i = 1:3
    subplot(3,2,2*i-1)
    [f1,x1i] = ksdensity(ot50(i+3,:));
    area(x1i,f1,'Facecolor',[1 0.5 0.5]);
    xlabel(xlab(i));
    ylabel(strcat('Pr'),'rot',0);
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [0 .2 0]);
    subplot(3,2,2*i)
    [f2,x2i] = ksdensity(enkf50(i+3,:));
    area(x2i,f2,'Facecolor',[0.5 0.7 0.5]);
    xlabel(xlab(i));
    ylabel(strcat('Pr'),'rot',0);
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [0 .1 0]);
end


% Plot the orbit first
% To plot the orbit we need real simulation data





