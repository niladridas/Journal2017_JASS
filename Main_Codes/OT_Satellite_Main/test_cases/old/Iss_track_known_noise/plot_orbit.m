% Plot orbit around the globe
clc;close all;clear;
fileID = fopen('25544.txt'); % TLE source, can just load a single line of TLE
date_mee = tle_parse(fileID,'0');
% The first observed value is used as the starting point to generate
% synthetic observations
% Generate synthetic real state values
% One data point observed each day

tp = 4; % PLOT specification when the point 




Nobs = tp-1;% Number of total observations
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


[t,mee_set] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),[0 2*Tsince],real_states(end,:));
real_orbit = mee_set;



%%
% TODO: Hard coding
real_eci = mee2eci(real_orbit((424-16):(424+15),:));



% ellipsoid_fit(real_eci(:,1:3))

ot80all = load('ot_all80.mat');
enkf80all = load('enkf_all80.mat');
load('real_states');


ot80 = mee2eci(ot80all.x_otall(:,:,tp)')';
enkf80 = mee2eci(enkf80all.x_enkfall(:,:,tp)')';

figure
% Plot the orbit
plot3(real_eci(:,1),real_eci(:,2),real_eci(:,3),'black');
hold on
% Plot the exact state at that point
xyzreal_pos = mee2eci(real_states(tp,:));
b = plot3(xyzreal_pos(1),xyzreal_pos(2),xyzreal_pos(3));
b.Marker = 'o';
b.MarkerSize = 6;
b.MarkerEdgeColor = 'blue';
b.MarkerFaceColor = 'blue';
hold on
% Plot the mean of OT
xyz_pos = mean(ot80(1:3,:),2);
a = plot3(xyz_pos(1),xyz_pos(2),xyz_pos(3));
a.Marker = 'o';
a.MarkerSize = 6;
a.MarkerEdgeColor = 'red';
a.MarkerFaceColor = 'red';
hold on
% Plot all the OT ensemble
all_ot = ot80(1:3,:)';
scatter3(all_ot(:,1),all_ot(:,2),all_ot(:,3),'red','filled');
hold on
% Plot mean of Enkf
xyz_pos1 = mean(enkf80(1:3,:),2);
c = plot3(xyz_pos1(1),xyz_pos1(2),xyz_pos1(3));
c.Marker = 'o';
c.MarkerSize = 6;
c.MarkerEdgeColor = [0.5 0.7 0.5];
c.MarkerFaceColor = [0.5 0.7 0.5];
hold on
% Plot all the Enkf ensemble
all_enkf = enkf80(1:3,:)';
scatter3(all_enkf(:,1),all_enkf(:,2),all_enkf(:,3),'green','filled');
hold on
r =  6378.135;
[x,y,z] = sphere;
s5 = surf(x*r, y*r, z*r,'FaceColor',[0.7 0.7 0.9]);
axis equal;
grid on

