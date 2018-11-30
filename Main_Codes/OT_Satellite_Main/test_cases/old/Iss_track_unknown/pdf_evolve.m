% plot pdf evolution
% main file
clc;close all;clear
fileID = fopen('25544.txt'); % TLE source, can just load a single line of TLE
date_mee = tle_parse(fileID,'0');
% The first observed value is used as the starting point to generate
% synthetic observations
% Generate synthetic real state values
% One data point observed each day
Nobs = 20;% Number of total observations
Tsince = 24*60*60; % One whole day
mu = 398600.4418;
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
real_states = zeros(Nobs,6);
real_states(1,:) = date_mee(1,7:12);
for i = 1:(Nobs-1)
   [t,mee_set] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),[0 Tsince],real_states(i,:));
   real_states(i+1,:) = mee_set(end,:);
end
% Last variable is an angle
real_states(:,6) = mod(real_states(:,6),2*pi);


% Fix noise structure
diag_sigma_init = [10 abs(real_states(1,2:5)).*0.15]; % GUESS
diag_sigma_init = (diag_sigma_init.*diag_sigma_init); % variance
SIGMA_init = diag(diag_sigma_init);
kappa_init = 66; % 
L_var = 1- (besseli(1,kappa_init)/besseli(0,kappa_init)); % corresponding variance in L
% TODO: noise pdf in ECI to MEE
% generate Noise
noise_r5 = mvnrnd(zeros(1,5),SIGMA_init,Nobs); % First R5 elements
noise_c1 = circ_vmrnd(0, kappa_init, Nobs);
noise_gvm = [noise_r5,noise_c1]; 
% obs_states 
obs_states = real_states + noise_gvm;

settle_t = 5;
%%
samples = 70;
%
% EnKF Parameters
H = eye(6); % All state are assumed to be measurable
R_enkf = diag([diag_sigma_init,L_var]);

equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
%% Known sensor noise EnKF


% OT Filtering essential functions
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);
N1 = 100;% Number of discrete steps in between 0 and 1
% % since point 1 is undefined we have taken 0 to 1-e-2
int_fun = @(x) bessel_C(x,N1);
weight = @(x,y) weight_newcal(x,y,int_fun);
% Known noise
% weight = @(x,y) weights_cal(x,y,kappa_init,SIGMA_init);

% Date structure for storing the estimated values

x_est_enkf = zeros(Nobs,6);
x_enkfall = zeros(6,samples,Nobs);

x_est_OT = zeros(Nobs,6);
x_otall = zeros(6,samples,Nobs);


init_r5 = mvnrnd(real_states(1,1:5),0.5*SIGMA_init,samples); % First R5 elements
init_c1 = circ_vmrnd(real_states(1,6), kappa_init, samples);
X_init_enkf = [init_r5,init_c1]'; 
X_init_OT = X_init_enkf;

for j = 1:Nobs
    if(j~=1)
        for i = 1:length(X_init_enkf(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 Tsince],X_init_enkf(:,i));
            X_init_enkf(:,i) = x_temp(end,:)';
        end
    end
    X_init_enkf = EnKF_filter(X_init_enkf,obs_states(j,:)',H,R_enkf);
    x_enkfall(:,:,j) = X_init_enkf;
    x_est_enkf(j,:) = mean(X_init_enkf,2)';
    x_est_enkf(j,6) = mod(x_est_enkf(j,6),2*pi);
    
    % calculate data based R_enkf
    if(j > settle_t)
        error = x_est_enkf(settle_t-3:j,:) - obs_states(settle_t-3,:);
        R_enkf = cov(error);
    end
end

for j = 1:Nobs
    if(j~=1)
        for i = 1:length(X_init_OT(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 Tsince],X_init_OT(:,i));
            X_init_OT(:,i) = x_temp(end,:)';
        end
    end
    X_init_OT = OT_filter(X_init_OT,obs_states(j,:)',cost,weight,OT_constantshdl);
    x_otall(:,:,j) = X_init_OT;
    x_est_OT(j,:) = mean(X_init_OT,2)';
    x_est_OT(j,6) = mod(x_est_OT(j,6),2*pi); 
end

% 
% mat_name = strcat('enkf_',num2str(samples),'.mat');
% save(mat_name, 'x_est_enkf');
% mat_name = strcat('ot_',num2str(samples),'.mat');
% save(mat_name, 'x_est_OT');


mat_name = strcat('enkf_all',num2str(samples),'.mat');
save(mat_name, 'x_enkfall');
mat_name = strcat('ot_all',num2str(samples),'.mat');
save(mat_name, 'x_otall');

