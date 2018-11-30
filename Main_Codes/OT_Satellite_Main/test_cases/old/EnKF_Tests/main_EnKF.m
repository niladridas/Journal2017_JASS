% EnKF test case for the orbit
% Author: Niladri Das
% Affiliation: UQLab, Aerospace Engineering, TAMU
% Date: 4 May 2017
clc
close
clear
% EnKF based filtering is tested on a full satellite model
% The model measured values are obtained as TLE elements from NORAD
%
% VARIABLES INITIALISED AT THE START
samples = 10;
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

% Step 2:
% Sample points from the R5XS state space
% The first 5 states are sampled from a normal distribution.
% TODO: The choice of variance should have some motivation
MU = date_mee(1,7:11);
% NEEDS PROPER MOTIVATION
SIGMA = diag([0.05,0.0001,0.0001,0.0001,0.0001]); 
init_r5 = mvnrnd(MU,SIGMA,samples); % First R5 elements
kappa = 2;% Needs motivation
init_c1 = circ_vmrnd(date_mee(1,12), kappa, samples);
X_init = [init_r5,init_c1]';

%% The initial data samples for the EnKF filter is generated 
% exactly same as that of the OT
% Need to figure out how the develop the measurement model for the EnKF
% We have assumed that we are measuring all the states in the output
% So our H is an identity matrix

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

% Parameters for EnKF
H = eye(6);
R_enkf = diag([0.05,0.0001,0.0001,0.0001,0.0001,0.001]);% FIX THE END VALUE


for j = 1:length(t1)
    if(j~=1)
        t_int = linspace(0,t1(j));
        for i = 1:length(X_init(1,:))
            [~,x_temp] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),t_int,X_init(:,i));
            X_init(:,i) = x_temp(end,:)';
        end
    end
    measured_output = date_mee(j,7:end)';
    X_init = EnKF_filter(X_init,measured_output,H,R_enkf);
    x_est(j,:) = mean(X_init,2)';
end
x_est_enkf = x_est;
save x_est_enkf x_est_enkf;

