% Test tleparser.m
clc
clear
format long
fid = fopen('tlevallado.tle');
FLAG = 0;% Since only one TLE module needs to be read
[date,dot_n,ddot_n,bstar,inclination,rightascension,...
    essentricity,argperi,meanano,meanmotion, totalrev,fid] = tleparser(fid,FLAG);
% Don't worry about the year. The data tested here is before 2000. So if
% the result is 2093 it really means 1993.
% The 20** scheme is implemented since the data which will be used from now
% on is taken after 2000
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
% Vallado References:
date_vallado = 1993;
dot_n_vallado = 2*7.889e-5*2*pi/((24*60)^2);
ddot_n_vallado = 0;
b_star_vallado = 0.00010529;
i_vallado = 51.6190*pi/180;
rightasc_vallado = 13.3340*pi/180;
e_vallado = 0.0005770;
arg_peri_vallado = 102.5680*pi/180;
meanano_vallado = 257.5950*pi/180;
meanotion_vallado = 15.59114070*2*pi/(24*60);

error = [abs(dot_n_vallado-dot_n);
        abs(ddot_n_vallado-ddot_n);
        abs(b_star_vallado-bstar);
        abs(i_vallado-inclination);
        abs(rightasc_vallado-rightascension);
        abs(e_vallado-essentricity);
        abs(arg_peri_vallado-argperi);
        abs(meanano_vallado-meanano);
        abs(meanotion_vallado-meanmotion)];
if norm(error)==0
    disp('TLE CONVERSION TEST PASSED')
else
    disp('TLE PARSER NOT WORKING');
end