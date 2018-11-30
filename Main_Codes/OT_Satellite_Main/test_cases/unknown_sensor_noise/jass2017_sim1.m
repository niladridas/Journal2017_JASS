% Testing EnKF and OT for one repitition and for varying sample size
% Normalised dynamics is not used here
% Author: Niladri Das
% Affiliation: UQ Lab TAMU Aerospace Engineering
% V1 Date: 13 June 2017
% V2 Date: 20 Nov 2018
close all;clear,clc;
load('./data/real_states.mat');
load('./data/prop_params.mat');
load('./data/delT.mat');
load('./data/mesurementnoise.mat');
load('./data/obs_states.mat');

load('./data/s.mat');

% We have to generate initial samples. So we choose the first observation
% as the mean around which the initial sample set will be generated.
% First observation is used to generate initial samples surrounding it
MU = real_states(1,1:5);
MU_theta = real_states(1,6); % The real_states theta data is already inside the [0,2*pi] domain.
% Choosing sample size of samples in each iteration
samples = 50;
n_obs = size(obs_states,1);
% GENERATING INITIAL SAMPLES
% Initial samples are generated from a Gauss von Mises Distribution 
% to compare the performance with OT Filter.

% ODE parameters for dynamics
% prop_params is already loaded
% OT Filtering necessary functions
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);    
% Vector containing time difference between each observation.
% At a particular point, it gives the time of propagation till the next
% observation
t1 = zeros(1,n_obs);
t1(1) = 0;
for i = 1:(n_obs-1)
    t1(i+1) = Tsince;
end
% Dynamics
% equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
% Parameters for EnKF
H = eye(6); % All state are assumed to be measurable
% The measurement model of EnKF requires a covariance matrix
% Hence the first 5 covariance values are kept to be same as that of the
% Gauss von Mises distribution. The last element (6,6) is choosen such that
% the covariance value matches that of the von Mises variance. This is done
% to compare with OT.


% Normalization
% ISS: 92.75 minutes
% The first 5 states are on R and the last state is an angle
temp_ob = abs(real_states(1,1:5));
norm_x =  max(temp_ob,[],1);
norm_matrix = diag([norm_x,1]);
t_norm = 92.75*60*2; % in seconds
% Normalised Time differences
t1 = t1./t_norm;
% Normalised Dynamics
% Normalised Dynamics
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
norm_dyn = @(t,x) norm_equinoctial_dyn(t,x,norm_matrix,t_norm,equinoc_dyn);
repeat = 10;

% Estimates
OT_all_repeat = zeros(n_obs,6,repeat);
EnKF_all_repeat = zeros(n_obs,6,repeat);
% sigma
OT_allsig_repeat = zeros(n_obs,6,repeat);
EnKF_allsig_repeat = zeros(n_obs,6,repeat);
% Run time
OT_time = zeros(n_obs,repeat);
EnKF_time = OT_time;
for n = 1:repeat
    kappa = mesurementnoise.kappa;% kappa for initial data generation is assumed to be half than that of measurement. Variance in theta is 1/sqrt(kappa)
    SIGMA = 10*diag(mesurementnoise.Pmee); % sigma for initial data generation is assumed to be twice than that of measurement.
    last_entry = 1- (besseli(1,kappa)/besseli(0,kappa));
    rng(s{n});
    init_r5 = mvnrnd(MU,SIGMA,samples); % First R5 elements
    rng(s{n});
    init_c1 = circ_vmrnd(MU_theta, kappa, samples);
    X_init_enkf = [init_r5,init_c1]'; 
    X_init_OT = X_init_enkf;
    % Estimated state variables using EnKF is stored here
    x_est_enkf = zeros(n_obs,6);
    x_est_OT = x_est_enkf; 
    % Estimated state covariances using EnKF is stored here
    x_sig_enkf = zeros(n_obs,6);
    x_sig_OT = x_sig_enkf;    
    % Normalised Initial samples
    X_init_OT = (norm_matrix)\X_init_OT;
    X_init_enkf = (norm_matrix)\X_init_enkf;   
    % Scaling SIGMA
    % Kappa remains same
    SIGMA = (norm_matrix(1:5,1:5)*norm_matrix(1:5,1:5))\SIGMA;
    N1 = 100;% Number of discrete steps in between 0 and 1
% since point 1 is undefined we have taken 0 to 1-e-2
    int_fun = @(x) bessel_C(x,N1);
    weight = @(x,y) weight_newcal(x,y,int_fun);
    % TO-DO
    R_enkf = zeros(6);
    R_enkf(1:5,1:5) = 100*SIGMA;% FIX THE END VALUE
    R_enkf(6,6) = 100*last_entry;
    % normalizing observations
    all_obs = (norm_matrix)\(obs_states');
    % MAIN  LOOP STARTS
    for j = 1:n_obs
        if(j~=1)
            for i = 1:length(X_init_enkf(1,:))
                [~,x_temp] = ode15s(norm_dyn,[0 t1(j)],X_init_enkf(:,i));
                % After propagation we should bring back the angle between
                % 0 and 2pi radians
                X_init_enkf(:,i) = x_temp(end,:)';
                X_init_enkf(6,i) = mod(X_init_enkf(6,i),2*pi);
            end
        end
        measured_output = all_obs(:,j);
        alpha = 0.3;
        tic;
        [X_init_enkf,R_enkf] = EnKF_filter_unknownnew(X_init_enkf,measured_output,H,R_enkf,alpha);
        toc;
        EnKF_time(j,n) = toc;
        % Calculating Estimates
        x_est_enkf(j,:) = mean(X_init_enkf,2)';% State Estimate
        x_est_enkf(j,6) = mod(x_est_enkf(j,6),2*pi);
        % Calculate sigma
        X_init_enkf(6,:) = mod(X_init_enkf(6,:),2*pi);
        x_sig_enkf_tmp = cov(X_init_enkf');
        x_sig_enkf(j,:) = sqrt(diag(x_sig_enkf_tmp)');
    end
    % Re-normalizing
    x_est_enkf = (norm_matrix*(x_est_enkf'))';
    x_sig_enkf = (norm_matrix*(x_sig_enkf'))';
    % MAIN OT LOOP STARTS
    for j = 1:n_obs
        if(j~=1)
            for i = 1:length(X_init_OT(1,:))
                [~,x_temp] = ode15s(norm_dyn,[0 t1(j)],X_init_OT(:,i));
                % After propagation we should bring back the angle between
                % 0 and 2pi radians
                X_init_OT(:,i) = x_temp(end,:)';
                X_init_OT(6,i) = mod(X_init_OT(6,i),2*pi);
            end
        end
        % OT
        measured_output = all_obs(:,j);
        tic;
        X_init_OT = OT_filter(X_init_OT,measured_output,cost,weight,OT_constantshdl);
        toc;
        OT_time(j,n) = toc;
        % Calculate Estimate
        x_est_OT(j,:) = mean(X_init_OT,2)';
        x_est_OT(j,6) = mod(x_est_OT(j,6),2*pi);
        % Calculate Covariances
        X_init_OT(6,:) = mod(X_init_OT(6,:),2*pi);
        x_sig_OT_tmp = cov(X_init_OT');
        x_sig_OT(j,:) = sqrt(diag(x_sig_OT_tmp)');    
    end
    % Re-normalizing    
    x_est_OT = (norm_matrix*(x_est_OT'))';
    x_sig_OT = (norm_matrix*(x_sig_OT'))'; 
    
OT_all_repeat(:,:,n) =  x_est_OT;
EnKF_all_repeat(:,:,n) = x_est_enkf; 
OT_allsig_repeat(:,:,n) = x_sig_OT;
EnKF_allsig_repeat(:,:,n) = x_sig_enkf;
end
file_OT = sprintf('./data/all_OT_data%i.mat',samples);
file_EnKF = sprintf('./data/all_EnKF_data%i.mat',samples);
save(file_OT,'OT_all_repeat','OT_allsig_repeat','OT_time');
save(file_EnKF,'EnKF_all_repeat','EnKF_allsig_repeat','EnKF_time');



