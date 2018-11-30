function x_analysed_EnKF = EnKF_filter(X_init,measured_output,H,R_enkf)
    M = size(measured_output,1);
    N = size(X_init,2);
    state_size = size(X_init,1);
    % First Calculate variance of forcast error
    x_forcast_mean_EnKF = X_init*(1/N)*ones(N); % Each colum contains the mean
    % ensemble perturbation matrix
    x_purterbation_EnKF = X_init - x_forcast_mean_EnKF;
    % ensemble co-variance matrix 
    P_ef_EmKF = x_purterbation_EnKF*x_purterbation_EnKF'/(N-1);
    x_analysed_EnKF = zeros(state_size,N);
    for k = 1:N
%         rng(s{k});
%         temp_enkf1 = inv(H*P_ef_EmKF*H'+R_enkf);
%         if(det(temp_enkf1)==0)
%             disp('Error')
%         end
%         x_analysed_EnKF(:,k) = X_init(:,k) + ...
%         (P_ef_EmKF*H')*(temp_enkf1)*(measured_output-H*X_init(:,k)+ mvnrnd(zeros(1,M),R_enkf)');
          x_analysed_EnKF(:,k) = X_init(:,k) + ...
          (P_ef_EmKF*H')*((H*P_ef_EmKF*H'+R_enkf)\(measured_output-H*X_init(:,k)+ mvnrnd(zeros(1,M),R_enkf)'));
    end
end