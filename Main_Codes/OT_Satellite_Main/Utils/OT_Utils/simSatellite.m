clear; clc;
addpath('matlab_circular/');
d2r = pi/180;

%% Single Initial Condition
% LEO Parameters
Re = 6378.1363; %km radius of earth
h = 900; % height of satellite km
R0 = (Re + h)/Re; % Normalized
th0 = 5*d2r;
V = 7.4; %km/s speed of satellite
Tp = 102.8*60; % Orbit Period in sec
Vb = Re/Tp;   % Velocity normalizing factor.

Rd = 0;  % normalized rate
thd = (V/Vb)/R0; %rad/s
x0 = [R0;Rd;th0;thd];
nOrbits = 20;

% Simulation of the states
J2flag = 0;
[t1,x1] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);

% Generating Observations
% parentpath = pwd();
% cd matlab_circular\

% R sigma
% Change in weight_cal
R_sigma = 10^-2; 
% Change in weight_cal
kappa = 5;
r_error = normrnd(0,R_sigma,size(t1,1),1);
theta_error = circ_vmrnd(0, kappa, size(t1,1));
% cd(parentpath);
Y = [(x1(:,1)+r_error)';(x1(:,3)+theta_error)'];

figure(1)
subplot(121)
plot(t1,x1(:,1),'r',t1,Y(1,:)','b');
grid on
subplot(122)
plot(t1,x1(:,3),'r',t1,Y(2,:)','b');
grid on

% Number of the particles
M = 40;
% Initialization of the states
r_init_samples = normrnd(x1(1,1),0.0001,M,1);
rdot_init_samples = normrnd(x1(1,2),0.0001,M,1);
% parentpath = pwd();
% cd matlab_circular\
kappa = 2;
theta_init_samples = circ_vmrnd(x1(1,3), kappa, M);
% cd(parentpath);
thetadot_init_samples = normrnd(x1(1,4),0.0001,M,1);

% OT Filtering necessary functions
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);
weight = @(x,y) weights_cal(x,y);

% Accumulate inital states
X_init = [r_init_samples';rdot_init_samples';theta_init_samples';thetadot_init_samples'];

% Estimated states
x_est = zeros(size(t1,1),4);
% Main loop over the time
for j = 1:size(t1,1)
    if(j~=1)
        t_before = t1(j-1);
        t_now = t1(j);
        for i = 1:length(X_init(1,:))
            [~,x_temp] = ode45(@satelliteDynamics,[t_before t_now],X_init(:,i),[],J2flag,Tp,Re);
            X_init(:,i) = x_temp(end,:)';
        end
    end
    measured_output = Y(:,j);
    X_init = OT_filter(X_init,measured_output,cost,weight,OT_constantshdl);
    x_est(j,:) = mean(X_init,2)';
end
%%
figure(2)
subplot(221)
% scatter(t1,x1(:,1),'r');
% hold on;
% scatter(t1,x_est(:,1),'b');
plot(t1,x1(:,1),'r',t1,x_est(:,1),'b');
grid on
subplot(222)
% scatter(t1,x1(:,3),'r');
% hold on
% scatter(t1,x_est(:,3),'b');
plot(t1,x1(:,3),'r',t1,x_est(:,3),'b');
grid on
subplot(223)
% scatter(t1,x1(:,2),'r');
% hold on
% scatter(t1,x_est(:,2),'b');
plot(t1,x1(:,2),'r',t1,x_est(:,2),'b');
grid on
subplot(224)
% scatter(t1,x1(:,4),'r');
% hold on;
% scatter(t1,x_est(:,4),'b');
plot(t1,x1(:,4),'r',t1,x_est(:,4),'b');
grid on


