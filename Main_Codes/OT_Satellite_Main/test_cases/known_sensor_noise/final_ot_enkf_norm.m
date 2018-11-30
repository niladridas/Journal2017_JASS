% Testing EnKF and OT for one repitition and for varying sample size
% Normalised dynamics is not used here
% Author: Niladri Das
% Affiliation: UQ Lab TAMU Aerospace Engineering
% V1 Date: 13 June 2017
% V2 Date: 20 Nov 2018
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
kappa = 17;% This one is for the von Mises distribution
% The first 5 values of the variance is 5% of the first observed value
diag_sigma = [10 abs(MU(2:5)).*0.20]; % 5% error is the standard deviation
diag_sigma = (diag_sigma.*diag_sigma); % variance
SIGMA = diag(diag_sigma); % For first 5 states
% The kappa and the Sigma parameters are assumed to be same both for the
% initial sample generation and also for the measurement model
% TODO: Initial sample variance should be larger than the measurement
% variance
% Use two different model: one for initial sample generation and one for
% measurement model
% repeat = 50;
% OT_all_repeat = zeros(size(date_mee,1),6,repeat);
% EnKF_all_repeat = zeros(size(date_mee,1),6,repeat);
% for n = 1:repeat
    init_r5 = mvnrnd(MU,SIGMA,samples); % First R5 elements
    init_c1 = circ_vmrnd(date_mee(1,12), kappa, samples);
    X_init_enkf = [init_r5,init_c1]'; 
    X_init_OT = X_init_enkf;
    % Estimated state variables using EnKF is stored here
    x_est_enkf = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
    x_est_OT = x_est_enkf;
    % OT Filtering necessary functions
    cost = @(x) distance_matrix(x);
    OT_constantshdl = @(x) OT_constants(x);
    N1 = 100;% Number of discrete steps in between 0 and 1

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
    % equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
    % Parameters for EnKF
    H = eye(6); % All state are assumed to be measurable
    % The measurement model of EnKF requires a covariance matrix
    % Hence the first 5 covariance values are kept to be same as that of the
    % Gauss von Mises distribution. The last element (6,6) is choosen such that
    % the covariance value matches that of the von Mises variance. This is done
    % to compare with OT.
    last_entry = 1- (besseli(1,kappa)/besseli(0,kappa));

    
    
    % Normalization
    % ISS: 92.75 minutes
    % The first 5 states are on R and the last state is an angle
    temp_ob = abs(date_mee(:,7:11));
    norm_x =  max(temp_ob,[],1);
    norm_matrix = diag([norm_x,1]);% since the last element is an angle
    t_norm = 92.75*60*2; % in seconds
    % Normalised Time differences
    t1 = t1./t_norm;
    % Normalised Dynamics
    % Normalised Dynamics
    equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
    norm_dyn = @(t,x) norm_equinoctial_dyn(t,x,norm_matrix,t_norm,equinoc_dyn);
                    
    % Normalised Initial samples
    X_init_OT = (norm_matrix)\X_init_OT;
    X_init_enkf = (norm_matrix)\X_init_enkf;
    
    % Scaling SIGMA
    % Kappa remains same
    SIGMA = (norm_matrix(1:5,1:5)*norm_matrix(1:5,1:5))\SIGMA;
    weight = @(x,y) weights_cal(x,y,kappa,SIGMA);
    R_enkf = zeros(6);
    R_enkf(1:5,1:5) = SIGMA;% FIX THE END VALUE
    R_enkf(6,6) = last_entry;
    
    % normalizing observations
    all_obs = (norm_matrix)\(date_mee(:,7:end)');


    % MAIN  LOOP STARTS
    for j = 1:size(date_mee,1)
        if(j~=1)
            for i = 1:length(X_init_enkf(1,:))
                [~,x_temp] = ode15s(norm_dyn,[0 t1(j)],X_init_enkf(:,i));
                % After propagation we should bring back the angle between
                % 0 and 2pi radians
                X_init_enkf(:,i) = x_temp(end,:)';
                X_init_enkf(6,i) = mod(X_init_enkf(6,i),6*pi);
            end
        end
        measured_output = all_obs(:,j);
        X_init_enkf = EnKF_filter(X_init_enkf,measured_output,H,R_enkf);
        x_est_enkf(j,:) = mean(X_init_enkf,2)';% State Estimate
        x_est_enkf(j,6) = mod(x_est_enkf(j,6),6*pi);
    end
    % Re-normalizing
    x_est_enkf = (norm_matrix*(x_est_enkf'))';
    save all_enkf_data x_est_enkf

    % MAIN OT LOOP STARTS
    for j = 1:size(date_mee,1)
        if(j~=1)
            for i = 1:length(X_init_OT(1,:))
                [~,x_temp] = ode15s(norm_dyn,[0 t1(j)],X_init_OT(:,i));
                % After propagation we should bring back the angle between
                % 0 and 2pi radians
                X_init_OT(:,i) = x_temp(end,:)';
                X_init_OT(6,i) = mod(X_init_OT(6,i),6*pi);
            end
        end
        % OT
        measured_output = all_obs(:,j);
        X_init_OT = OT_filter(X_init_OT,measured_output,cost,weight,OT_constantshdl);
        x_est_OT(j,:) = mean(X_init_OT,2)';
        x_est_OT(j,6) = mod(x_est_OT(j,6),6*pi);
    end
    x_est_OT = (norm_matrix*(x_est_OT'))';
    save all_OT_data x_est_OT
% OT_all_repeat(:,:,n) =  x_est_OT;
% EnKF_all_repeat(:,:,n) = x_est_enkf; 
% end
