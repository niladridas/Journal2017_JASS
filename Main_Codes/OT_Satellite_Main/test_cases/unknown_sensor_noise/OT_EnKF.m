% Author: Niladri Das
% Affiliation: UQLab, Aerospace Engineering, TAMU
% Date: 9 June 2017

% We compare OT and EnKF here with known noise characteristics but with no
% process noise
% clc;close;clear
% Modes of operation:
% 1. Sample size is varied and is passed as an argument: 'samples', sample
% vector. We compute the variation in RMSE error along with the computation
% time
% 2. Initial condition is varied N times and statistics is generated 
% 3. Only one run and OT and EnKF are compared for one specific sample size
function OT_EnKF(~)
    close all;clear all,clc;

    sample_vector = 10:10:50;
    N = 1; % Repetition to generate statistics
    % Step 1:
    % Measured MEE observations saved as mat file
    % TODO 1:
    % NOTE: The first element in the MEE set of observations is much larger
    % than all other elements. This problem needs to be addressed using some
    % normalization constant.
    % The ISS MEE data are all saved in date_mee.mat file
    % load date_mee.mat file
    load('date_mee.mat');

    % Step 2:
    % Sample points from the R5XS state space
    % The first 5 states are sampled from a normal distribution.
    % TODO: The choice of variance should have some motivation
    MU = date_mee(1,7:11); % First observation is used to generate initial samples surrounding it
    % NEEDS PROPER MOTIVATION
    kappa = 2;% This one is for the von Mises distribution
    % This one is for the multivariate gaussian distribution
    SIGMA = diag([0.05,0.000001,0.000001,0.000001,0.000001]);
    % TODO: change this to reflect 5% of the mean value as the 3 sigma
    % bound 
    
    
    % Step 3:
    % OT Filtering necessary functions
    cost = @(x) distance_matrix(x);
    OT_constantshdl = @(x) OT_constants(x);
    N1 = 100;% Number of discrete steps in between 0 and 1
    % since point 1 is undefined we have taken 0 to 1-e-2
    % int_fun = @(x) bessel_C(x,N1);
    % WEIGHT FUNCTION 
    % This can be either the marginalized one 
    % This can also be the normal known noise characteristics
    % weight = @(x,y) weight_newcal(x,y,int_fun);
    weight = @(x,y) weights_cal(x,y);
    % Form the time stamp from the date_mee matrix using time_diff
    t1 = [];
    t1(1) = 0;
    for i = 1:(size(date_mee,1)-1)
        t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
    end
    % ODE parameters for dynamics
    mu = 398600.4418;
    R_e = 6378.135;
    J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
    prop_params = struct('mu',mu,'R_e',R_e,'J',J);
    
    count = 0;
    t_run_ot = [];
    t_run_enkf = [];
    x_est_all_ot = [];
    x_est_all_enkf = [];
    % Parameters for EnKF
    H = eye(6);
    R_enkf = diag([0.05,0.0001,0.0001,0.0001,0.0001,0.01]);% FIX THE END VALUE    
    red_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
    % MAIN LOOP STARTS HERE
    for l=1:length(sample_vector) % varying the sample size
        for k = 1:N % repeating with different initial conditions
            samples = sample_vector(l);
            init_r5 = mvnrnd(MU,SIGMA,samples); % First R5 elements
            init_c1 = circ_vmrnd(date_mee(1,12), kappa, samples);
            X_init_ot = [init_r5,init_c1]';
            X_init_enkf = [init_r5,init_c1]';
            % Step 4:
            % Estimated states
            x_est_ot = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
            x_est_enkf = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
            count = count + 1;
            t_ot_total = 0;
            t_enkf_total = 0;
            for j = 1:size(date_mee,1)
                if(j~=1)
                    t_int = linspace(0,t1(j));
                    for i = 1:length(X_init_ot(1,:))
                        [~,x_temp] = ode45(@(t,x) equinoctial_dyn(t,x,prop_params),t_int,X_init_ot(:,i));
                        X_init_ot(:,i) = x_temp(end,:)';
                    end
                end
                % OT
                measured_output = date_mee(j,7:end)';
                tic;
                X_init_ot = OT_filter(X_init_ot,measured_output,cost,weight,OT_constantshdl);
                t_ot = toc;
                t_ot_total = t_ot+t_ot_total;
                x_est_ot(j,:) = mean(X_init_ot,2)';
            end
            t_run_ot = [t_run_ot;t_ot_total];
            x_est_all_ot(:,:,count) = x_est_ot;
            %
            for j = 1:size(date_mee,1)
                if(j~=1)
                    t_int = linspace(0,t1(j));
                    for i = 1:length(X_init_enkf(1,:))
                        [~,x_temp] = ode45(@(t,x) equinoctial_dyn(t,x,prop_params),t_int,X_init_enkf(:,i));
                        X_init_enkf(:,i) = x_temp(end,:)';
                    end
                end
                % EnKF
                measured_output = date_mee(j,7:end)';
                tic;
                X_init_enkf = EnKF_filter(X_init_enkf,measured_output,H,R_enkf);
                t_enkf = toc;
                t_enkf_total = t_enkf+t_enkf_total;
                x_est_enkf(j,:) = mean(X_init_enkf,2)';
            end
            t_run_enkf = [t_run_enkf;t_enkf_total];
            x_est_all_enkf(:,:,count) = x_est_enkf;
        end
    end
    save all_data t_run_ot x_est_all_ot t_run_enkf x_est_all_enkf
end

