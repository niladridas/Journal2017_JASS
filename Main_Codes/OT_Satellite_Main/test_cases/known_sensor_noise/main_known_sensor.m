% Author: Niladri Das
% Affiliation: UQ Lab, TAMU, Aerospace Engineering
% Date: 23 June 2017

% Know sensor characteristics
% OT and EnKF comparison for synthetic data and also for actual observed
% NORAD data of ISS
clc;clear;close all;

% DYNAMICS
% ODE parameters for dynamics
mu = 398600.4418;
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);

% OBSERVATION MODEL
% This has two modes
% 1. Synthetic observation is generated
% 2. Real Observed data is available

% For Synthetic data we can choose any interval and simulate for any period
% of time. Here we want to have the same time period as the real observation
% were taken. Hence we import the real observed data.
load('date_mee.mat');
MU = date_mee(1,7:11);
t1 = zeros(1,size(date_mee,1));
t1(1) = 0;
for i = 1:(size(date_mee,1)-1)
    t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
end

str = input('Real or Synthetic R/S [R]: ','s');
if isempty(str)
    str = 'R';
end

% Common variable
kappa = 66;
diag_sigma = [10 abs(MU(2:5)).*0.2]; % GUESS
diag_sigma = (diag_sigma.*diag_sigma); % variance
SIGMA = diag(diag_sigma); % For first 5 states

if str=='S'
    % MODE 1:
    % Synthetic Data generation
    % We generate synthetic actual states
    syn_state = zeros(size(date_mee,1),6);
    syn_state(1,:) = date_mee(1,7:12);
    for i= 2:size(date_mee,1)
    [~,x_temp] = ode45(equinoc_dyn,[0 t1(i)],syn_state(i-1,:)');
    syn_state(i,:) = x_temp(end,:);
    end
    syn_state(:,6) = mod(syn_state(:,6),2*pi);
    % We generate synthetic observations
    % generate noise in the first 5 variables
    noise_5 = mvnrnd([0;0;0;0;0],SIGMA,size(date_mee,1));
    noise_L = circ_vmrnd(0, kappa, size(date_mee,1));
    measured_output = syn_state +[noise_5,noise_L];
    measured_output(:,6) = mod(measured_output(:,6),2*pi);
end

if str == 'R'
    % MODE 2:
    measured_output = date_mee(:,7:12);
    measured_output(:,6) = mod(measured_output(:,6),2*pi);
end

% INITIAL SAMPLE GENERATION
% Sample size selection
samples = 50;
%
diag_sigma_init = [20 abs(MU(2:5)).*0.3]; % GUESS
diag_sigma_init = (diag_sigma_init.*diag_sigma_init); % variance
SIGMA_init = diag(diag_sigma_init);
kappa_init = 66;
%
init_r5 = mvnrnd(MU,SIGMA_init,samples); % First R5 elements
init_c1 = circ_vmrnd(mod(date_mee(1,12),2*pi), kappa_init, samples);
X_init_enkf = [init_r5,init_c1]'; 
X_init_OT = X_init_enkf;

% OT Filtering essential functions
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);
weight = @(x,y) weights_cal(x,y,kappa,SIGMA);

% EnKF Parameters
H = eye(6); % All state are assumed to be measurable
% The measurement model of EnKF requires a covariance matrix
% Hence the first 5 covariance values are kept to be same as that of the
% Gauss von Mises distribution. The last element (6,6) is choosen such that
% the covariance value matches that of the von Mises variance. This is done
% to compare with OT.
last_entry = 1- (besseli(1,kappa)/besseli(0,kappa));
R_enkf = diag([diag_sigma,last_entry]);% FIX THE END VALUE

% Date structure for storing the estimated values
x_est_enkf = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
x_est_OT = x_est_enkf;

tic
for j = 1:size(date_mee,1)
    if(j~=1)
        for i = 1:length(X_init_enkf(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_enkf(:,i));
            % After propagation we should bring back the angle between
            % 0 and 2pi radians
            X_init_enkf(:,i) = x_temp(end,:)';
            X_init_enkf(6,i) = mod(X_init_enkf(6,i),2*pi);
        end
    end
    X_init_enkf = EnKF_filter(X_init_enkf,measured_output(j,:)',H,R_enkf);
    x_est_enkf(j,6) = mod(x_est_enkf(j,6),2*pi);
    x_est_enkf(j,:) = mean(X_init_enkf,2)';% State Estimate
end
toc
EnKF_time = toc-tic;
tic
% MAIN OT LOOP STARTS
for j = 1:size(date_mee,1)
    if(j~=1)
        for i = 1:length(X_init_OT(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_OT(:,i));
            % After propagation we should bring back the angle between
            % 0 and 2pi radians
            X_init_OT(:,i) = x_temp(end,:)';
            X_init_OT(6,i) = mod(X_init_OT(6,i),2*pi);
        end
    end
    % OT
    X_init_OT = OT_filter(X_init_OT,measured_output(j,:)',cost,weight,OT_constantshdl);
    x_est_OT(j,6) = mod(x_est_OT(j,6),2*pi);
    x_est_OT(j,:) = mean(X_init_OT,2)';
end
toc
OT_time = toc-tic;
%%
t_spot = zeros(1,length(t1));
for i = 1:length(t1)
    t_spot(i) = sum(t1(1:i));
end

% t_spot is in seconds
% conversting it into days
t_spot = t_spot./(60*60*24);

ylabel_names = ['p','f','g','h','k','L'];

figure(1)
for i = 1:6
    % Fix the color
    % Fix the fonts 
    % Fix the symbol on the lines like scatter plot
    % Give the title
    % Give the footer note
    subplot(3,2,i);
    a = plot(t_spot,x_est_OT(:,i));
    a.Color = 'red';
    a.LineStyle = '-';
    a.LineWidth = 2;
    a.Marker = 'o';
    a.MarkerSize = 4;
    a.MarkerEdgeColor = 'red';
    a.MarkerFaceColor = 'red';
    xlabel('Time (Days)','FontName','sans-serif','FontSize',20,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',20,'FontWeight','Bold');
%     set(gca,'Xtick',0:3:15);

    set(gca,'fontWeight','bold','Xtick',0:15);
    
    hold on;
    b = plot(t_spot,x_est_enkf(:,i));
    b.Color = 'blue';
    b.LineStyle = '-';
    b.LineWidth = 2;
    b.Marker = 'o';
    b.MarkerSize = 4;
    b.MarkerEdgeColor = 'blue';
    b.MarkerFaceColor = 'blue';
%     xlabel('Time in Days');
    
    
    hold on
    if str == 'R'
    c = plot(t_spot,measured_output(:,i));
    elseif str == 'S'
    c = plot(t_spot,syn_state(:,i));   
    end
    c.Color = [0 0.5 0];
    c.LineStyle = '-';
    c.LineWidth = 2;
    c.Marker = 'o';
    c.MarkerSize = 4;
    c.MarkerEdgeColor = [0 0.5 0];
    c.MarkerFaceColor = [0 0.5 0];
%     xlabel('Time in Days')
    if str == 'R'
    legend('OT','EnKF','Obs')
    elseif str == 'S'
       legend('OT','EnKF','True')     
    end
    grid on
end
if str == 'R'
temp_datatype = 'Observed States';
elseif str == 'S'
temp_datatype = 'True States';  
end
title_str = strcat('Comparison of OT and EnKF for',{' '},num2str(samples),' samples with ',{' '}, temp_datatype);
% suptitle(title_str);

axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.55, 0, title_str, 'FontSize', 20', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
% set(gca,'dataaspectratio',[1,5,1]);
file_name = strcat('samples',num2str(samples),str);
fig = gcf;
fig.PaperPositionMode = 'auto';
print(file_name,'-dpdf','-fillpage');
% print(file_name,'-depsc','-r0','-loose');
% saveas(fig,file_name,'epsc');
% 
% img = getframe(gcf);
% imwrite(img.cdata, [file_name, '.eps']);


