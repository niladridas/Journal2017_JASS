% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 26th April 2017

% test tle parser and sgp4 together

clc
clear
format long
% Take the first tle module 
% fid = fopen('test_sgp.txt');
fid = fopen('test_sgp.txt');
FLAG=0;
[date,dotn,ddotn,bstar,inclination,rightascension,e,argperi,...
    meanano, meanmotion, totalrev,fid] = tleparser(fid,FLAG);
% Step 2: Pass it through the sgp4 module
% (t-t_o) is the time since epoch
% 1: n_o: SGP type "mean" mean motion at epoch
%   L2C53-63 rev/day
% 2: e_o: the “mean” eccentricity at epoch
%   L2C27-33
% 3: i_o: the "mean" inclination at epoch
%   L2C9-16
% 4: M_o: the "mean" mean anomaly at epoch
%   L2C44-51 (degrees)
% 5: w_o: the "mean" argument of perigee at epoch
%   L2C35-42
% 6: Ohm_o: the "mean" logitude of ascending node at epoch
%   L2C18-25
% 7*: n_o_dot: the time rate of change of “mean” mean motion at epoch
%   L1C34-43 n_o_dot/2
% 8*: n_o_ddot: the second time rate of change of “mean” mean motion at epoch
%   L1C45-52 n_o_ddot/6
% 9*: B_star: the SGP4 type drag coefficient
%   L1C54-61
%%
t_o = 0;
t = 0;
n_o = meanmotion;
e_o = e;
i_o = inclination; 
M_o = meanano; 
w_o = argperi; 
Ohm_o = rightascension; 
B_star = bstar;
% Step 3: Constant elements for sgp4 module
% Constant Elements:
XKMPER = 6378.135;
a_E = 1.0;
k_2 = 5.413080e-4;
k_e = 0.74366916e-1; % First error here sqrt(GM)
A_3_0= -0.0000025323;
s = 1.01222928;
qo_s4 = 1.88027916e-9;
B = B_star; % Vallado
k_4 = 0.62098875e-6;
sgp4_const_elements = struct('XKMPER',XKMPER,'a_E',a_E,'k_2',k_2,'k_e',k_e,'A_3_0',A_3_0,...
    's',s,'qo_s4',qo_s4,'B',B,'k_4',k_4);
%
tle_parsed_elements = struct('n_o',n_o,'e_o',e_o,'i_o',i_o,'M_o',M_o,...
    'w_o',w_o,'Ohm_o',Ohm_o,'B_star',B_star,'t_o',t_o,'t',t);
% The Length unit is in earth's radius and time in minutes
pos_vel = sgp4(tle_parsed_elements,sgp4_const_elements);% in er and mins
pos_km = pos_vel(1,1:3);%.*6378.135;
vel_km_sec = pos_vel(1,4:end);%.*(6378.135/60);
% The real output data from the celestrac report are the following
pos_km_real = [2328.97048951,-5995.22076416,1719.97067261];
vel_km_sec_real = [2.91207230,-0.98341546,-7.09081703];
% Error 
pos_error = norm(pos_km_real - pos_km);
vel_error = norm(vel_km_sec_real-vel_km_sec);
disp('Error between the Fortran code and Matlab Implementation')
fprintf('The position error is:%f Km\n',pos_error);
fprintf('The Velocity error is:%f Km/sec\n',vel_error);


