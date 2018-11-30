%% Two Body Dynamics in Equinoctial Orbital Elements
% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 18th April 2017

% Two body dynamics in Equinocial Elements
% Source: A set of Modified Eqinoctial Orbit Elements - Walker, Ireland and Owens (also see the errata)
clear
clc
% Test 1 :
% Test data are used from the above paper
% Initial Lagrangian Values
a = 24419.205;
e = 0.726683;
i = 27*pi/180;
w = 0;
Ohm  = 0;
nu = 0;
Init_Lag_val = [a;e;i;w;Ohm;nu];
% Intial Equinoctial Values
Init_El = Lag2El(Init_Lag_val);
% p = a*(1-e^2);
% f = e*cos(w+Ohm);
% g = e*sin(w+Ohm);
% h = tan(i/2)*cos(Ohm);    
% k = tan(i/2)*sin(Ohm);
% L = Ohm + w + nu;

% Constant Components
mu = 398603.2;
R_e = 6378.165;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];

prop_params = struct('mu',mu,'R_e',R_e,'J',J);
% Init_El = [p;f;g;h;k;L]; 
% options = odeset('RelTol',1e-5);
[t,El_set] = ode45(@(t,El) equinoctial_dyn(t,El,prop_params),[0 172800],Init_El);%,options);

%% Final Lagrangian Elements
to_match = El2Lag(El_set(end,:)');
to_match_degree = to_match.*[1;1;180/pi;180/pi;180/pi;180/pi];
save to_match.mat to_match;


a_f =  24331.433;
e_f =  0.72557888;
i_f =  26.988272;
w_f =  1.199160;
Ohm_f = 359.280136;
nu_f  = 186.307367;

Lag_f = [a_f;e_f;i_f;w_f;Ohm_f;nu_f];
diff_val = to_match.*[1;1;180/pi;180/pi;180/pi;180/pi] - Lag_f
save diff_val.mat diff_val;