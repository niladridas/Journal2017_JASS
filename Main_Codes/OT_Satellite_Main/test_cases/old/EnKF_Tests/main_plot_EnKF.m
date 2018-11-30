% Author: Niladri Das
% Affiliation: UQ Lab, Aerospace Engineering, TAMU
% Date: 29 May 2017

% This file plots the data from the OT filtering procedure
clear, clc, close all;
% Load the mat date files first
load('x_est_enkf.mat');
% x_est: variable created in the work space
% dim of x_est: sample points X 6 (5 on real and the end on circular)
load('date_mee.mat')
% The last 6 columns are the observed 6 variables
obs_var = date_mee(:,7:12);
pre_var = x_est_enkf;

n = size(x_est_enkf,1);

% Plot
figure(1)
subplot(2,3,1)
plot(pre_var(:,1));
hold on
% scatter(1:n,pre_var(:,1));
plot(obs_var(:,1));
% scatter(1:n,obs_var(:,1));
legend('pre1','pre2','ob1','ob2')
grid on

subplot(2,3,2)
plot(pre_var(:,2));
hold on
% scatter(1:n,pre_var(:,2));
plot(obs_var(:,2));
% scatter(1:n,obs_var(:,2));
grid on

subplot(2,3,3)
plot(pre_var(:,3));
hold on
% scatter(1:n,pre_var(:,3));
plot(obs_var(:,3));
% scatter(1:n,obs_var(:,3));
grid on


subplot(2,3,4)
plot(pre_var(:,4));
hold on
% scatter(1:n,pre_var(:,4));
plot(obs_var(:,4));
% scatter(1:n,obs_var(:,4));
grid on

subplot(2,3,5)
plot(pre_var(:,5));
hold on
% scatter(1:n,pre_var(:,5));
plot(obs_var(:,5));
% scatter(1:n,obs_var(:,5));
grid on

% The sixth variable needs to be in between 0 and 2*pi radians
pre_var(:,6) = mod(pre_var(:,6),2*pi);
obs_var(:,6) = mod(obs_var(:,6),2*pi);
% In degrees
pre_var(:,6) = pre_var(:,6)*180/pi;
obs_var(:,6) = obs_var(:,6)*180/pi;

subplot(2,3,6)
plot(pre_var(:,6));
hold on
% scatter(1:n,pre_var(:,6));
plot(obs_var(:,6));
% scatter(1:n,obs_var(:,6));
grid on

% Plot of position and velocity for the observed and the predited values

% figure(2)





