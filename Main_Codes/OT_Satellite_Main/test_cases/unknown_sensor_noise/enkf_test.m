% Testing EnKF when it fails for a specific state values
% Normalised dynamics is not used here
% Author: Niladri Das
% Affiliation: UQ Lab TAMU Aerospace Engineering
% Date: 13 June 2017
close all;clear,clc;

% DATA LOADING
% Loading observation data
load('date_mee.mat');
% We have to generate initial samples. So we choose the first observation
% as the mean around which the initial sample set will be generated.
% First observation is used to generate initial samples surrounding it
MU = date_mee(1,7:11); 
% Choosing sample size of samples in each iteration
samples = 40;

% GENERATING INITIAL SAMPLES
% Initial samples are generated from a Gauss von Mises Distribution 
% to compare the performance with OT Filter.
kappa = 2;% This one is for the von Mises distribution
diag_sigma = [0.05,0.0001,0.0001,0.0001,0.0001];
SIGMA = diag(diag_sigma); % For first 5 states
% The kappa and the Sigma parameters are assumed to be same both for the
% initial sample generation and also for the measurement model
% TODO: Initial sample variance should be larger than the measurement
% variance
% Use two different model: one for initial sample generation and one for
% measurement model
init_r5 = mvnrnd(MU,SIGMA,samples); % First R5 elements
init_c1 = circ_vmrnd(date_mee(1,12), kappa, samples);
X_init_enkf = [init_r5,init_c1]'; 
% Estimated state variables using EnKF is stored here
x_est_enkf = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));

% OT Filtering necessary functions
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);
N1 = 100;% Number of discrete steps in between 0 and 1
weight = @(x,y) weights_cal(x,y);

% Vector containing time difference between each observation.
% At a particular point, it gives the time of propagation till the next
% observation
t1 = zeros(1,size(date_mee,1));
t1(1) = 0;
for i = 1:(size(date_mee,1)-1)
    t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
end

% ODE parameters for dynamics
mu = 398600.4418;
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);

% Dynamics
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
% Parameters for EnKF
H = eye(6); % All state are assumed to be measurable
% The measurement model of EnKF requires a covariance matrix
% Hence the first 5 covariance values are kept to be same as that of the
% Gauss von Mises distribution. The last element (6,6) is choosen such that
% the covariance value matches that of the von Mises variance. This is done
% to compare with OT.
% TODO: fix the end value
R_enkf = diag([diag_sigma,0.01]);% FIX THE END VALUE


% MAIN LOOP STARTS
for j = 1:size(date_mee,1)
    if(j~=1)
        for i = 1:length(X_init_enkf(1,:))
            [~,x_temp] = ode15s(equinoc_dyn,[0 t1(j)],X_init_enkf(:,i));
            X_init_enkf(:,i) = x_temp(end,:)';
        end
    end
    measured_output = date_mee(j,7:end)';
    X_init_enkf = EnKF_filter(X_init_enkf,measured_output,H,R_enkf);
    x_est_enkf(j,:) = mean(X_init_enkf,2)';% State Estimate
end
save all_enkf_data x_est_enkf
