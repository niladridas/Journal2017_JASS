% Author: Niladri Das
% Affiliation: UQLab, Aerospace Engineering, TAMU
% Date: 14 May 2017

% We use marginalization technique for the satellite problem
% New weight_cal function
% Named as weight_newcal


% In main function
% int_fun = @(x) bessel_C(x,N)

function W = weight_newcal(x_state,z,int_fun)
    % NOTE: This weight function is specific to the case of satellite
    % problem where the first 5 states are on R and the last state (6th) is
    % in S
    M = size(x_state,2); % number of samples
    W_tmp = zeros(1,M);
    
    for i = 1:M
        e_r =z(1:5,1)- x_state(1:5,i);
        e_theta = z(6,1)-x_state(6,i);
        C0 = (prod(e_r.^(2)))^(-0.5);
        C1 = int_fun(e_theta);
        W_tmp(1,i) = C0*C1; % First R5 elements
%         if W_tmp(1,i)<0
%             disp('Error')
%         end
    end
    Sum_W = sum(W_tmp);
    W = W_tmp./Sum_W;
end 