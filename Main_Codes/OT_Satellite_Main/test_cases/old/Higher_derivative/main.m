% Author: Niladri Das
% Date: 20 Sep 2017
clc;clear;close;
% Step 1:
% Generate 100 samples from a Gaussian distribution with given mean and
% variance
mu4 = 2.401;
P4 = 1.96304502;
rng(0,'twister');
data4 = mu4 + sqrt(P4)*randn(100,1);
%
mu3 = 3.43;
P3 = 1.965398;
rng(0,'twister');
data3 = mu3 + sqrt(P3)*randn(100,1);
% Use Optimal Transport to calculate the T matrix
