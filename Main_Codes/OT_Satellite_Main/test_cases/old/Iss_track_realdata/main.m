% main file
clc;close all;clear
fileID = fopen('25544.txt'); % TLE source, can just load a single line of TLE
date_mee = tle_parse(fileID,'0');

% Pre-processing 
date_mee = date_mee(~any(isnan(date_mee),2),:); % First delete repeting rows
% Delete repreating dates
[C,ia,ic] = unique(date_mee(:,1:6),'rows','stable');
date_mee = date_mee(ia,:);
date_mee = date_mee(1:20,:);

% Time stamp and real data
t1 = zeros(1,size(date_mee,1)); % Time difference between each succesive observation
for i = 1:(size(date_mee,1)-1)
    t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
end
measured_output = date_mee(:,7:12); % Only obseved value (no time stamp)

% The first observed value is used as the starting point to generate
% synthetic observations
% Generate synthetic real state values
% One data point observed each day
Nobs = size(date_mee,1);% Number of total observations

MU = date_mee(1,7:12);

mu = 398600.4418;
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
% real_states = zeros(Nobs,6);
% real_states(1,:) = date_mee(1,7:12);
% for i = 1:(Nobs-1)
%    [t,mee_set] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),[0 Tsince],real_states(i,:));
%    real_states(i+1,:) = mee_set(end,:);
% end
% % Last variable is an angle
% real_states(:,6) = mod(real_states(:,6),2*pi);

save real_states

% Fix noise structure
diag_sigma_init = [10 abs(MU(1,2:5)).*0.15]; % GUESS
diag_sigma_init = (diag_sigma_init.*diag_sigma_init); % variance
SIGMA_init = diag(diag_sigma_init);
kappa_init = 66; % 
L_var = 1- (besseli(1,kappa_init)/besseli(0,kappa_init)); % corresponding variance in L



%%
samples = 80;
settle_t = 5;
%
% EnKF Parameters
H = eye(6); % All state are assumed to be measurable
R_enkf = 5*diag([diag_sigma_init,L_var]);

equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
%% Unknown sensor noise EnKF


% OT Filtering essential functions
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);
N1 = 100;% Number of discrete steps in between 0 and 1
% since point 1 is undefined we have taken 0 to 1-e-2
int_fun = @(x) bessel_C(x,N1);
weight = @(x,y) weight_newcal(x,y,int_fun);
% % Known noise
% weight = @(x,y) weights_cal(x,y,kappa_init,SIGMA_init);

% Date structure for storing the estimated values

x_est_enkf = zeros(Nobs,6);
x_enkf_cov = zeros(Nobs,6);
x_est_OT = zeros(Nobs,6);
x_ot_cov = zeros(Nobs,6);

init_r5 = mvnrnd(MU(1,1:5),0.5*SIGMA_init,samples); % First R5 elements
init_c1 = circ_vmrnd(MU(1,6), kappa_init, samples);
X_init_enkf = [init_r5,init_c1]'; 
X_init_OT = X_init_enkf;

for j = 1:Nobs
    if(j~=1)
        for i = 1:length(X_init_enkf(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_enkf(:,i));
            X_init_enkf(:,i) = x_temp(end,:)';
        end
    end
    X_init_enkf = EnKF_filter(X_init_enkf,measured_output(j,:)',H,R_enkf);
    x_est_enkf(j,:) = mean(X_init_enkf,2)';
    x_enkf_cov(j,:) = diag(cov(X_init_enkf'))';
    if sum(x_enkf_cov(j,:)<0)>0
        error('Incorrect calculation of covar matrix')
    end
    x_est_enkf(j,6) = mod(x_est_enkf(j,6),2*pi);
    % calculate data based R_enkf
    if(j > settle_t)
        error = x_est_enkf(2:j,:) - measured_output(2:j,:);
        R_enkf = cov(error);
    end
end

for j = 1:Nobs
    if(j~=1)
        for i = 1:length(X_init_OT(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_OT(:,i));
            X_init_OT(:,i) = x_temp(end,:)';
        end
        X_init_OT = OT_filter(X_init_OT,measured_output(j,:)',cost,weight,OT_constantshdl);
    end
    x_est_OT(j,:) = mean(X_init_OT,2)';
    x_ot_cov(j,:) = diag(cov(X_init_OT'))';
    if sum(x_ot_cov(j,:)<0)>0
        error('Incorrect calculation of covar matrix')
    end
    x_est_OT(j,6) = mod(x_est_OT(j,6),2*pi); 
end


measured_output = measured_output(settle_t+1:end,:);
x_est_enkf = x_est_enkf(settle_t+1:end,:);
x_est_OT = x_est_OT(settle_t+1:end,:);
x_enkf_cov = x_enkf_cov(settle_t+1:end,:);
x_ot_cov = x_ot_cov(settle_t+1:end,:);

save measured_output;


mat_name = strcat('enkf_',num2str(samples),'.mat');
save(mat_name, 'x_est_enkf');
mat_name = strcat('ot_',num2str(samples),'.mat');
save(mat_name, 'x_est_OT');
mat_name = strcat('enkf_cov',num2str(samples),'.mat');
save(mat_name, 'x_enkf_cov');
mat_name = strcat('ot_cov',num2str(samples),'.mat');
save(mat_name, 'x_ot_cov');  