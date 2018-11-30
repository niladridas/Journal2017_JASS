% OT test case for the orbit
% Author: Niladri Das
% Affiliation: UQLab, Aerospace Engineering, TAMU
% Date: 4 May 2017

% The ode method for propagating the states is taking time.
% Normalising the state space and the time will speed up the process
% and reduce the accumulating error due to the ode we are using.


clc
close
clear
% OT based filtering is tested on a full satellite model
% The model measured values are obtained as TLE elements from NORAD
%
% VARIABLES INITIALISED AT THE START
samples = 50;
%
% Step 1:
% Measured MEE observations saved as mat file
% TODO 1:
% NOTE: The first element in the MEE set of observations is much larger
% than all other elements. This problem needs to be addressed using some
% normalization constant.

% The ISS MEE data are all saved in date_mee.mat file
% load date_mee.mat file
load('date_mee.mat');

% Approach
% Take the observations
% Find the normalizing constant with respect to these values
% Fix time of orbit of the satellite as the normalizing constant for the
% time

% ISS: 92.75 minutes
% The first 5 states are on R and the last state is an angle
% temp_ob = abs(date_mee(:,7:11));
% norm_x =  max(temp_ob,[],1);
% norm_matrix = diag([norm_x,1]);% since the last element is an angle
% t_norm = 92.75*60; % in seconds
%%


% Step 2:
% Sample points from the R5XS state space
% The first 5 states are sampled from a normal distribution.
% TODO: The choice of variance should have some motivation
MU = date_mee(1,7:11);
% NEEDS PROPER MOTIVATION
 

SIGMA = diag([0.05,0.0001,0.0001,0.0001,0.0001]);

% 
% N = 1;
% OT_stat = [];
% for k = 1:N
    init_r5 = mvnrnd(MU,SIGMA,samples); % First R5 elements
    kappa = 2;% Needs motivation
    init_c1 = circ_vmrnd(date_mee(1,12), kappa, samples);
    X_init = [init_r5,init_c1]';

    % Step 3:
    % OT Filtering necessary functions
    cost = @(x) distance_matrix(x);
    OT_constantshdl = @(x) OT_constants(x);
    N1 = 100;% Number of discrete steps in between 0 and 1
    % since point 1 is undefined we have taken 0 to 1-e-2
    int_fun = @(x) bessel_C(x,N1);
    weight = @(x,y) weight_newcal(x,y,int_fun);
    % weight = @(x,y) weights_cal(x,y);

    % Step 4:
    % Estimated states
    x_est = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
    % Form the time stamp from the date_mee matrix using time_diff
    t1 = [];
    t1(1) = 0;
    for i = 1:(size(date_mee,1)-1)
        t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
    end
    %%
    % Main loop over the time
    mu = 398600.4418;
    R_e = 6378.135;
    J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
    prop_params = struct('mu',mu,'R_e',R_e,'J',J);
    
%     % Normalised Dynamics
%     equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
%     norm_dyn = @(t,x) norm_equinoctial_dyn(t,x,norm_matrix,t_norm,equinoc_dyn);
%     
%     % Normalised Initial samples
%     X_init = (norm_matrix)\X_init;
%     
%     % Normalised Time differences
%     t1 = t1./t_norm;

    for j = 1:length(t1)
        if(j~=1)
            t_int = linspace(0,t1(j));
            tic
            for i = 1:length(X_init(1,:))
                [~,x_temp] = ode45(@(t,x) equinoctial_dyn(t,x,prop_params),t_int,X_init(:,i));
%                 [~,x_temp] = ode45(norm_dyn,t_int,X_init(:,i));
                X_init(:,i) = x_temp(end,:)';
            end
            toc
        end
%         X_init = norm_matrix*X_init;
        measured_output = date_mee(j,7:end)';
        X_init = OT_filter(X_init,measured_output,cost,weight,OT_constantshdl);
        x_est(j,:) = mean(X_init,2)';
        % Re-normalise again to be used in the loop
%         X_init = (norm_matrix)\X_init;
    end
%     OT_stat(:,:,k) = x_est;
% end

save x_est x_est

