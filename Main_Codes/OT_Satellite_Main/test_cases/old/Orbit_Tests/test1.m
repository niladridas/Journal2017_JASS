 % Test1: Calculate prediction accuracy with sgp4 for ISS
% 'iss.txt': TLEs collected at different time
% Author: Niladri Das
% Affiliation: UQ Lab, Aerospace Engineering, TAMU
clc
close all
clear 
% Constant Elements:
XKMPER = 6378.135;
a_E = 1.0;
k_2 = 5.413080e-4;
k_e = 0.74366916e-1; % First error here sqrt(GM)
A_3_0= -0.0000025323;
s = 1.01222932;
qo_s4 = ((120-78)*a_E/XKMPER)^4;
k_4 = 0.62098875e-6;
sgp4_const_elements = struct('XKMPER',XKMPER,'a_E',a_E,'k_2',k_2,'k_e',k_e,'A_3_0',A_3_0,...
    's',s,'qo_s4',qo_s4,'k_4',k_4);
% Step 1:
% Read the TLE elements and index them with proper time when these
% measurements were done
% The tleparser reads one TLE module and returns the file id to be used
% again for reading
% Select the text file to read
% Step 2:
% Use the SGP4 model to calculate the exact ECI elemets 
% These forms the measurement values
fid = fopen('iss.txt');
FLAG = 1;
date_pos_vel = [];
while fid~=-1
    [date,dot_n,ddot_n,bstar,inclination,rightascension,...
    essentricity,argperi,meanano,meanmotion, totalrev,fid] = tleparser(fid,FLAG);
    if fid == -1
        break
    end
    tle_parsed_elements = struct('n_o',meanmotion,'e_o',essentricity,'i_o',inclination,'M_o',meanano,...
    'w_o',argperi,'Ohm_o',rightascension,'B_star',bstar,'t_o',0,'t',0);
    % pos_vel : units in Km and secs
    date_pos_vel = [date_pos_vel;[date,sgp4(tle_parsed_elements,sgp4_const_elements)]];
end
%
% Plotting the data in 3D
x_pos = date_pos_vel(:,7);
y_pos = date_pos_vel(:,8);
z_pos = date_pos_vel(:,9);

figure(1) %3d plot
scatter3(x_pos,y_pos,z_pos,'r','filled')
hold on
r = XKMPER;
[x,y,z] = sphere;
surf(x*r, y*r, z*r,'FaceColor',[0 0.4 0.7]);
axis equal;
grid on

%%
% Convert all these ECI elements to equinoctial elements
mu = 398600.4418;
date_coe = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
date_coe(:,1:6) = date_pos_vel(:,1:6);
date_mee = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
date_mee(:,1:6) = date_pos_vel(:,1:6);

% Test
test_mee = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
test_mee(:,1:6) = date_pos_vel(:,1:6);

for i = 1:size(date_pos_vel,1)
    % First Convert to classical orbital element
    % use eci to coe
    date_coe(i,7:end) = eci2coe(date_pos_vel(i,7:12),mu);
    % Now convert from coe to mee
    date_mee(i,7:end) = coe2mee(date_coe(i,7:end));
    
    % test to see if the mee = eci2mee(mu, reci, veci) gives the same
    % result
%     test_mee(i,7:end) = eci2mee(mu, date_pos_vel(i,7:9), date_pos_vel(i,10:end));
end
% Check the difference between eci2mee and my code
% error = date_mee(:,7:end)-test_mee(:,7:end);

save date_mee date_mee
%%
% Test how accurate the equinoctial element based dynamics is
% Take on mee element, propagate it to the next time point and then find
% the error between the propagated and measured
% If there are N time steps the loop is for N-1
R_e = XKMPER;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
mee_predict = zeros((size(date_mee,1)-1), 6);
for i = 1:(size(date_mee,1)-1)
   % Call function to return the time difference in seconds
   Tsince = time_diff(date_mee(i,1:6),date_mee((i+1),1:6));
   [t,mee_set] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),[0 Tsince],date_mee(i,7:end));%,options);
   mee_predict(i,:) = mee_set(end,:);
end
mee_ob = date_mee(2:end,7:end);
save mee_ob mee_ob
save mee_predict mee_predict
%%
% Calculating the error in the observed and predicted MEE elements for ISS
error_e = mee_ob - mee_predict;
% The last column is an angle
for i = 1:size(error_e,1)
    error_e(i,6) = mod(error_e(i,6),2*pi);
    if error_e(i,6)>pi
        error_e(i,6) = 2*pi - error_e(i,6);
    end
end
% Normed error for the first column
for i = 1: size(error_e,1)
    error_e(i,1) = error_e(i,1)/6790;
end
figure(3)
for i = 1:6
subplot(2,3,i)
plot(error_e(:,i))
end
