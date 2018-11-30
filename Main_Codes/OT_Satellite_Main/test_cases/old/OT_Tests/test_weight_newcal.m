% Author: Niladri Das
% Affiliation: UQ Lab, Aerospace Engineering, TAMU
% Date: 18 May 2017

% This is a test file for weight_newcal function

clc;clear;close;
%
N = 100;
int_fun = @(x) bessel_C(x,N);
% Make dummy x_state and z
x_state = rand(6,10);
z = rand(6,1);
W = weight_newcal(x_state,z,int_fun);
