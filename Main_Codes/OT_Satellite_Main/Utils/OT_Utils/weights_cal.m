% Author: Niladri Das

function W = weights_cal(x_state,z,kappa,SIGMA)
    % NOTE: This weight function is specific to the case of satellite
    % problem where the first 5 states are on R and the last state (6th) is
    % in S
    M = size(x_state,2); % number of samples
    W_tmp = zeros(1,M);
    for i = 1:M
        [p,~] = circ_vmpdf(z(6,1)-x_state(6,i), 0, kappa);
%         display(p);
        MU = zeros(1,5);
        % NEEDS PROPER MOTIVATION
        W_tmp(1,i) = mvnpdf((z(1:5,1)-x_state(1:5,i))',MU,SIGMA)*p; %First R5 elements
    end
    Sum_W = sum(W_tmp);
    W = W_tmp./Sum_W;
end