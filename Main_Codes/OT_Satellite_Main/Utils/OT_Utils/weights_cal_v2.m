% Author: Niladri Das

function W = weights_cal_v2(x_state,z,kappa,SIGMA)
    % NOTE: This weight function is specific to the case of satellite
    % problem where the first 5 states are on R and the last state (6th) is
    % in S
    M = size(x_state,2); % number of samples
    W_tmp = zeros(1,M);
    for i = 1:M
        [p,~] = circ_vmpdf(z(2,1)-x_state(6,i), 0, kappa);
%         display(p);
        MU = zeros(1,1);
        % NEEDS PROPER MOTIVATION
        W_tmp(1,i) = normpdf((z(1,1)-x_state(1,i))',MU,sqrt(SIGMA))*p; %First R5 elements
    end
    Sum_W = sum(W_tmp);
    W = W_tmp./Sum_W;
end