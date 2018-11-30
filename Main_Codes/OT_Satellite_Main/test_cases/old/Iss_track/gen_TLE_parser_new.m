clc;close all;clear
FileID = fopen('35664.txt'); % TLE source
FLAG = 1; % To keep on reading
% Constant Elements:
sgp4_const_elements = struct('XKMPER',6378.135,'a_E',1.0,...
    'k_2',5.413080e-4,'k_e',0.74366916e-1,...
    'A_3_0',-0.0000025323,'s',1.01222932,...
    'qo_s4',((120-78)*1.0/6378.135)^4,'k_4',0.62098875e-6);
date_pos_vel = [];
while FileID~=-1
    [date,dot_n,ddot_n,bstar,inclination,rightascension,...
        essentricity,argperi,meanano,meanmotion,...
        totalrev,FileID] = tleparser(FileID,FLAG);
    if FileID == -1;break;end
    tle_parsed_elements = struct('n_o',meanmotion,...
        'e_o',essentricity,'i_o',inclination,'M_o',meanano,...
    'w_o',argperi,'Ohm_o',rightascension,'B_star',bstar,'t_o',0,'t',0);
    % pos_vel : units in Km and secs
    date_pos_vel = [date_pos_vel;[date,sgp4(tle_parsed_elements,sgp4_const_elements)]];
end % Data reading finished

% Eliminate some data so that we emulate a situation where observatons are
% taken after a long interval of time 

% date_pos_vel = flipud(date_pos_vel);
date_pos_vel = date_pos_vel(1:5:end,:);
save('date_pos_vel.mat','date_pos_vel')

% ECI to MEE elements
mu = 398600.4418;
date_coe = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
date_coe(:,1:6) = date_pos_vel(:,1:6);
date_mee = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
date_mee(:,1:6) = date_pos_vel(:,1:6);

for i = 1:size(date_pos_vel,1)
    % use eci to coe
    date_coe(i,7:end) = eci2coe(date_pos_vel(i,7:12),mu);
    % Now convert from coe to mee
    date_mee(i,7:end) = coe2mee(date_coe(i,7:end));
end
% The LAST COLUMN is an angle
date_mee(:,end) = mod(date_mee(:,end),2*pi);

save('date_mee.mat','date_mee');
% % Test how accurate the equinoctial element based dynamics is
% % Take on mee element, propagate it to the next time point and then find
% % the error between the propagated and measured
% % If there are N time steps the loop is for N-1
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
mee_predict = zeros(size(date_mee,1), 6);
mee_predict(1,:) = date_mee(1,7:12);
time_pt_days = [];
time_pt_days(1) = 0;
for i = 1:(size(date_mee,1)-1)
   % Call function to return the time difference in seconds
   Tsince = time_diff(date_mee(i,1:6),date_mee((i+1),1:6));
   time_pt_days(i+1) = time_pt_days(i)+(Tsince/(60*60*24));
   [t,mee_set] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),[0 Tsince],mee_predict(i,:));
   mee_predict(i+1,:) = mee_set(end,:);
end
mee_ob = date_mee(:,7:end);mee_ob(:,end) = mod(mee_ob(:,end),2*pi);
save('mee_ob.mat','mee_ob');
mee_predict(:,end) = mod(mee_predict(:,end),2*pi);
save('mee_predict.mat','mee_predict')
save('time_pt_days.mat','time_pt_days');
%% OT and EnKF predictions on real data
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
MU = date_mee(1,7:11);
t1 = zeros(1,size(date_mee,1));
t1(1) = 0;
for i = 1:(size(date_mee,1)-1)
    t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
end
% kappa = 66;
% diag_sigma = [20 abs(MU(2:5)).*0.3]; % GUESS
% diag_sigma = (diag_sigma.*diag_sigma); % variance
% SIGMA = diag(diag_sigma); % For first 5 states
measured_output = date_mee(:,7:12);
% measured_output(:,6) = mod(measured_output(:,6),2*pi);
% INITIAL SAMPLE GENERATION
rep = 1;
samples = 15;
diag_sigma_init = [50 abs(MU(2:5)).*0.5]; % GUESS
diag_sigma_init = (diag_sigma_init.*diag_sigma_init); % variance
SIGMA_init = diag(diag_sigma_init);

% Date structure for storing the estimated values
x_est_OT = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2),rep);
x_est_enkf = x_est_OT;

kappa_init = 66;
for m = 1:rep
    init_r5 = mvnrnd(MU,SIGMA_init,samples); % First R5 elements
    init_c1 = circ_vmrnd(mod(date_mee(1,12),2*pi), kappa_init, samples);
    X_init_OT = [init_r5,init_c1]'; 
    X_init_enkf = [init_r5,normrnd(date_mee(1,12),0.4,samples,1)]'; % 10 degree error as the sigma
    % OT Filtering essential functions
    cost = @(x) distance_matrix(x);
    OT_constantshdl = @(x) OT_constants(x);
    N1 = 100;% Number of discrete steps in between 0 and 1
    % since point 1 is undefined we have taken 0 to 1-e-2
    int_fun = @(x) bessel_C(x,N1);
    weight = @(x,y) weight_newcal(x,y,int_fun);
    % EnKF Parameters
    H = eye(6); % All state are assumed to be measurable
%     last_entry = 1- (besseli(1,kappa_init)/besseli(0,kappa_init));
    R_enkf = 10*diag([diag_sigma_init,0.4^2]);% FIX THE END VALUE
    % OT used marginalization 


    for j = 1:size(date_mee,1)
        if(j~=1)
            for i = 1:length(X_init_enkf(1,:))
                [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_enkf(:,i));
                X_init_enkf(:,i) = x_temp(end,:)';
            end
        end
        X_init_enkf(end,:) = mod(X_init_enkf(end,:),2*pi);
        X_init_enkf = EnKF_filter(X_init_enkf,measured_output(j,:)',H,R_enkf);
        x_est_enkf(j,:,m) = mean(X_init_enkf,2)';% State Estimate
        error_acc = measured_output(1:j,:) - x_est_enkf(1:j,:,m);
        if j >= 8 % Window that needs to be deleted
            R_enkf = cov(error_acc);
        end
    end
    for j = 1:size(date_mee,1)
        if(j~=1)
            for i = 1:length(X_init_OT(1,:))
                [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_OT(:,i));
                % After propagation we should bring back the angle between
                % 0 and 2pi radians
                X_init_OT(:,i) = x_temp(end,:)';
                X_init_OT(6,i) = mod(X_init_OT(6,i),2*pi);
            end
        end
        X_init_OT = OT_filter(X_init_OT,measured_output(j,:)',cost,weight,OT_constantshdl);
%         X_init_OT(6,:) = mod(X_init_OT(6,:),2*pi);
        x_est_OT(j,:,m) = mean(X_init_OT,2)';
    end
end
mat_name = strcat('enkf_repeat_',num2str(samples),'.mat');
save(mat_name, 'x_est_enkf');
mat_name = strcat('ot_repeat_',num2str(samples),'.mat');
save(mat_name, 'x_est_OT');






