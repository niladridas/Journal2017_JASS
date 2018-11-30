function [x_analysed_EnKF,R_new] = EnKF_filter_unknownnew(X_init,measured_output,H,R_old,alpha)
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
    e = measured_output - x_forcast_mean_EnKF(:,1);
%     b = e.^2/2;
%     inv_g_samples = 1./gamrnd(0.5*ones(1,6),b');
%     R_enkf = diag(inv_g_samples);
    R_new = alpha*R_old + (1-alpha)*(e*e'+H*P_ef_EmKF*H');
    for k = 1:N
%         temp_enkf1 = inv(H*P_ef_EmKF*H'+R_enkf);
%         if(det(temp_enkf1)==0)
%             disp('Error')
%         end
%         x_analysed_EnKF(:,k) = X_init(:,k) + ...
%         (P_ef_EmKF*H')*(temp_enkf1)*(measured_output-H*X_init(:,k)+ mvnrnd(zeros(1,M),R_enkf)');
          x_analysed_EnKF(:,k) = X_init(:,k) + ...
          (P_ef_EmKF*H')*((H*P_ef_EmKF*H'+R_new)\(measured_output-H*X_init(:,k)+ mvnrnd(zeros(1,M),R_new)'));
    end
end