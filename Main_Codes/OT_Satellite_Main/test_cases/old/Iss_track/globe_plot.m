% globe plot
clc;clear;close all;
fileID = fopen('25544.txt'); % TLE source
[date_mee,time_pt_days] = tle_parse(fileID,'0');
%
mu = 398600.4418;
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);
MU = date_mee(1,7:11);
t1 = zeros(1,size(date_mee,1));
t1(1) = 0;
for i = 1:(size(date_mee,1)-1)
    t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
end
measured_output = date_mee(:,7:12);
% INITIAL SAMPLE GENERATION
diag_sigma_init = [20 abs(MU(2:5)).*0.3]; % GUESS
diag_sigma_init = (diag_sigma_init.*diag_sigma_init); % variance
SIGMA_init = diag(diag_sigma_init);
kappa_init = 66;
%
samples = 20;
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
last_entry = 1- (besseli(1,kappa_init)/besseli(0,kappa_init));
R_enkf = 2*diag([diag_sigma_init,last_entry]);
% Date structure for storing the estimated values
x_est_OT = zeros(size(date_mee(:,7:end),1),size(date_mee(:,7:end),2));
x_est_enkf = x_est_OT;
x_cov_OT = zeros(6,6,size(date_mee,1));
x_cov_enkf = zeros(6,6,size(date_mee,1));

x_sam_OT = zeros(6,samples,size(date_mee,1));
x_sam_enkf = zeros(6,samples,size(date_mee,1));


% OT used marginalization 
for j = 1:size(date_mee,1)
    if(j~=1)
        for i = 1:length(X_init_enkf(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_enkf(:,i));
            X_init_enkf(:,i) = x_temp(end,:)';
        end
    end
    X_init_enkf = EnKF_filter(X_init_enkf,measured_output(j,:)',H,R_enkf);
    x_sam_enkf(:,:,j) = X_init_enkf; 
    x_est_enkf(j,:) = mean(X_init_enkf,2)';
%     x_cov_enkf(:,:,j) = cov(X_init_enkf'); 
    x_est_enkf(j,6) = mod(x_est_enkf(j,6),2*pi);
    error_acc = measured_output(1:j,:) - x_est_enkf(1:j,:);
    if j >= 8 % Window that needs to be deleted
        R_enkf = cov(error_acc);
    end
end
for j = 1:size(date_mee,1)
    if(j~=1)
        for i = 1:length(X_init_OT(1,:))
            [~,x_temp] = ode45(equinoc_dyn,[0 t1(j)],X_init_OT(:,i));
            X_init_OT(:,i) = x_temp(end,:)';
        end
    end
    X_init_OT = OT_filter(X_init_OT,measured_output(j,:)',cost,weight,OT_constantshdl);
    x_sam_OT(:,:,j) = X_init_OT;
%     x_est_OT(j,:) = mean(X_init_OT,2)';
%     x_cov_OT(:,:,j) = cov(X_init_OT'); 
%     x_est_OT(j,6) = mod(x_est_OT(j,6),2*pi); 
end
%%
% plot on globe
% choos a time point for comparison
x_eci_OT = zeros(6,samples,size(date_mee,1));
x_eci_enkf = zeros(6,samples,size(date_mee,1));
% first mee to eci
for i = 1:size(date_mee,1)
    for j = 1:samples
        x_eci_OT(:,j,i) = coe2eci(mee2coe(x_sam_OT(:,j,i)));
        x_eci_enkf(:,j,i) = coe2eci(mee2coe(x_sam_enkf(:,j,i)));
    end
end

% mean positions on the globe
x_eci_OTmean = squeeze(mean(x_eci_OT,2))';
x_eci_enkfmean = squeeze(mean(x_eci_enkf,2))';

% Checking if the conversions were correct
% norm_distOT = sqrt(sum(x_eci_OTmean(:,1:3).*x_eci_OTmean(:,1:3),2));
% plot(norm_distOT(8:end,1));
% hold on
% norm_distenkf = sqrt(sum(x_eci_enkfmean(:,1:3).*x_eci_enkfmean(:,1:3),2));
% plot(norm_distenkf(8:end,1));
% 15th data point
tp = 15;

% actual observed values
x_obs = coe2eci(mee2coe(date_mee(tp,7:12)));


x_pos = x_eci_OT(1,:,tp);
y_pos = x_eci_OT(2,:,tp);
z_pos = x_eci_OT(3,:,tp);


f= mvksdensity([x_pos',y_pos',z_pos'],[x_pos',y_pos',z_pos'],'bandwidth',[50,50,50]);

% hold on
% r = 6378.135;
% [x,y,z] = sphere;
% surf(x*r, y*r, z*r,'FaceColor',[1 1 1]);
% axis equal;
% grid on
f = f./sum(f);

figure(1) %3d plot
scatter3(x_pos,y_pos,z_pos,[],1-f,'filled');
hold on
a = plot3(x_obs(1),x_obs(2),x_obs(2));
a.Color = 'red';
a.LineStyle = '-';
a.LineWidth = 1;
a.Marker = 'o';
a.MarkerSize = 8;
a.MarkerEdgeColor = 'red';
a.MarkerFaceColor = 'red';
%axis tight
xlabel('x');
ylabel('y');
zlabel('z');

norm(x_eci_OTmean(tp,1:3) - x_obs(1:3))






