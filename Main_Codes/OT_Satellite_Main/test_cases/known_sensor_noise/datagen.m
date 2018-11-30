% Equinoctial elements are used for orbital state representation.
% The first element of EE is in Km.
clc;close all;clear
load('./data/date_mee.mat');
load('./data/s.mat');
% The starting MEE state
Nobs = 15;% Number of total observations
real_states = zeros(Nobs,6);
real_states(1,:) = date_mee(1,7:12);
real_states(:,1) = real_states(:,1)*10^3;
% Time interval between observations
Tsince = 24*60*60; % One whole day
mu = 398600.4418*10^9;
R_e = 6378.135*10^3;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
for i = 1:(Nobs-1)
   [t,mee_set] =radau(equinoc_dyn,[0 Tsince],real_states(i,:),options);
   real_states(i+1,:) = mee_set(end,:);
end
% Last variable is an angle
real_states(:,6) = mod(real_states(:,6),2*pi);
save('./data/real_states.mat','real_states');
%% Fix noise structure
% In MEE elements
% Assuming independence of Gauss and von Mises parameters.
diag_sigma_init = [10^8, 10^-6,10^-6,10^-6,10^-6];
% diag_sigma_init = [10^6, 10^-4,10^-4,10^-4,10^-4];
SIGMA_init = diag(diag_sigma_init);
kappa_init = 66; % 
L_var = 1- (besseli(1,kappa_init)/besseli(0,kappa_init)); % corresponding variance in L
mesurementnoise.mumee = zeros(1,5);
% mesurementnoise.Pmee = [10^9, 10^-4,10^-4,10^-4,10^-4];
mesurementnoise.Pmee =diag_sigma_init;
mesurementnoise.kappa = kappa_init;
% TODO: noise pdf in ECI to MEE
% generate Noise
rng(s{end});
noise_r5 = mvnrnd(zeros(1,5),SIGMA_init,Nobs); % First R5 elements
rng(s{end-1});
noise_c1 = circ_vmrnd(0, kappa_init, Nobs);
noise_gvm = [noise_r5,noise_c1]; 
% obs_states 
obs_states = real_states + noise_gvm;
save('./data/prop_params.mat','prop_params');
save('./data/delT.mat','Tsince');
save('./data/mesurementnoise.mat','mesurementnoise');
save('./data/obs_states.mat','obs_states');
%%
% What are the things that we have now from this simulation
% 1. The propagation parameters in meters and seconds.
% 2. The real mee states each separated by 1 day
% 3. Time difference btween two observations.
% 4. Measurement noise parameters
% 5. The observed MEE states
