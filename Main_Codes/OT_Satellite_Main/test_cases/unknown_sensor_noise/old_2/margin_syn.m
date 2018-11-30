% Testing EnKF and OT for one repitition and for varying sample size
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


% GENERATING INITIAL SAMPLES
% Initial samples are generated from a Gauss von Mises Distribution 
% to compare the performance with OT Filter.
kappa = 66;% This one is for the von Mises distribution
% The first 5 values of the varian+ce is 5% of the first observed value
kappa_init = 66;% This one is for the von Mises distribution


diag_sigma = [10 abs(MU(2:5)).*0.2]; % 
diag_sigma = (diag_sigma.*diag_sigma); % variance
SIGMA = diag(diag_sigma); % For first 5 states
diag_sigma_init = [20 abs(MU(2:5)).*0.3]; % 
diag_sigma_init = (diag_sigma_init.*diag_sigma_init); % variance
SIGMA_init = diag(diag_sigma_init);% For initial sample generation that should have larger variance than the observation model


% The kappa and the Sigma parameters are assumed to be same both for the
% initial sample generation and also for the measurement model
% TODO: Initial sample variance should be larger than the measurement
% variance
% Use two different model: one for initial sample generation and one for
% measurement model

% repeat = 5;
% OT_all_repeat = zeros(size(date_mee,1),6,repeat);
% EnKF_all_repeat = zeros(size(date_mee,1),6,repeat);
% for n = 1:repeat

% sample_list = 30:10:80;

% for t = 1:length(sample_list)

    samples = 80;
    init_r5 = mvnrnd(MU,SIGMA_init,samples); % First R5 elements
    init_c1 = circ_vmrnd(mod(date_mee(1,12),2*pi), kappa_init, samples);
    X_init_OT = [init_r5,init_c1]'; 
    % Estimated state variables using EnKF is stored here
    x_est_OT = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
    % OT Filtering necessary functions
    cost = @(x) distance_matrix(x);
    OT_constantshdl = @(x) OT_constants(x);
    
    N1 = 100;% Number of discrete steps in between 0 and 1
    % since point 1 is undefined we have taken 0 to 1-e-2
    int_fun = @(x) bessel_C(x,N1);
    weight = @(x,y) weight_newcal(x,y,int_fun);
    


    % Vector containing time difference between each observation.
    % At a particular point, it gives the time of propagation till the next
    % observation
    t1 = zeros(1,size(date_mee,1));
    t1(1) = 0;
    for i = 1:(size(date_mee,1)-1)
        t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
    end
    
    save t1;
    
    % ODE parameters for dynamics
    mu = 398600.4418;
    R_e = 6378.135;
    J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
    prop_params = struct('mu',mu,'R_e',R_e,'J',J);
    % Dynamics
    equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
    
    
    % We generate synthetic actual states
    syn_state = zeros(size(date_mee,1),6);
    syn_state(1,:) = date_mee(1,7:12);
    for i= 2:size(date_mee,1)
    [~,x_temp] = ode45(equinoc_dyn,[0 t1(i)],syn_state(i-1,:)');
    syn_state(i,:) = x_temp(end,:);
    end
    syn_state(:,6) = mod(syn_state(:,6),2*pi);
    
    save syn_state syn_state;
    % We generate synthetic observations
    % generate noise in the first 5 variables
    noise_5 = mvnrnd([0;0;0;0;0],SIGMA,size(date_mee,1));
    noise_L = circ_vmrnd(0, kappa, size(date_mee,1));
    
    measured_output = syn_state +[noise_5,noise_L];
    measured_output(:,6) = mod(measured_output(:,6),2*pi);
    
    save measured_output measured_output;
    %%
    



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
    OT_time = toc;
    save all_OT_data x_est_OT
    
    
    % RMSE ERROR

%     RMSE_ot(t,:) = [sqrt(mean((x_est_OT - syn_state).^2,1)),OT_time];  % Root Mean Squared Error
%     save error_time_ot RMSE_ot OT_time
%     save error_time_enkf RMSE_enkf EnKF_time
    
% end
% save RMSE_ot; 
