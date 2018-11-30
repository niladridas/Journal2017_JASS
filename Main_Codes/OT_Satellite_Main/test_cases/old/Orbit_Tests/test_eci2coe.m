% Test eci2coe
clc
clear
close
mu = 398600.4418;
eci_el = [6524.834,6862.875,6448.296,4.001327,5.533756,-1.976341];
coe_el = eci2coe(eci_el,mu);
% In degrees
disp('In Degrees')
coe_el = coe_el.*[1,1,180/pi,180/pi,180/pi,180/pi]
