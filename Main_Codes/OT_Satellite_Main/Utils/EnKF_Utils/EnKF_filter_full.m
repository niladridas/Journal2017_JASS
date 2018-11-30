function [filtered_data_enkf,Initial_data_total_enkf] = EnKF_filter(plant,noise_model,Tstamp,real_states,Observations,Initial_data)
    t_steps = size(Tstamp,2);
    state_size = size(real_states,2);
    N = size(Initial_data,1);
    filtered_data_enkf = zeros(t_steps,state_size);
    x_forcast_EnKF = Initial_data';
    
% %     noise_obs_params = struct('H',H,'mu_all_ch',mu_all_ch,...
% %     'var_all_ch',var_all_ch,'rho_all_ch',rho_all_ch);
%     plant.NoiseParams
%     noiseparams = struct('var_x',0.1,'var_dotx',0.1,'rho',0.95);
    R_e = plant.NoiseParams.var_all_ch;
    
%     g = str2func(noise_model);
%     R_e = ones(state_size);
    
   Initial_data_total_enkf = zeros(N,state_size,t_steps);
    for i = 1:t_steps
        % First Calculate variance of forcast error
        x_forcast_mean_EnKF = x_forcast_EnKF*(1/N)*ones(N); % Each colum contains the mean
        % ensemble perturbation matrix
        x_purterbation_EnKF = x_forcast_EnKF - x_forcast_mean_EnKF;
        % ensemble co-variance matrix 
        P_ef_EmKF = x_purterbation_EnKF*x_purterbation_EnKF'/(N-1);
        x_analysed_EnKF = zeros(state_size,N);
        H = plant.NoiseParams.H;
        % Noise Model
        for k = 1:N
            temp_enkf1 = inv(H*P_ef_EmKF*H'+R_e);
            x_analysed_EnKF(:,k) = x_forcast_EnKF(:,k) + ...
                (P_ef_EmKF*H')*(temp_enkf1)*(Observations(i,:)'-H*x_forcast_EnKF(:,k)+normrnd(0,sqrt(R_e))); % This is the change that was needed
        end
%         keyboard
        filtered_data_enkf(i,:) = mean(x_analysed_EnKF,2)';
%         e = (Observations(i,:)'-H*filtered_data_enkf(i,:)');
%         R_e =  g(i,R_e,e);       
        Initial_data_total_enkf(:,:,i) = x_analysed_EnKF'; 
        if i == t_steps
            break
        end
        for k = 1:N
%             temp_1 = syssimulate(x_analysed_EnKF(:,k)',(Tstamp(i+1)-Tstamp(i)));
            temp_1  = plant_dyn(x_analysed_EnKF(:,k)',(Tstamp(i+1)-Tstamp(i)),plant,'No_Noise'); 
            x_forcast_EnKF(:,k) = temp_1';
        end
%         t_samples = x_analysed_EnKF;
%         t_val = Tstamp(1,i);
%         File1 = 'filteredsamples_enkf.txt';
%         if exist(File1, 'file')
%             dlmwrite(File1,t_val,'-append') ;
%             dlmwrite(File1,t_samples,'-append') ;
%         else
%            dlmwrite(File1,t_val) ;
%            dlmwrite(File1,t_samples,'-append') ;
%         end   
    end
%     if(exist('filtered_enkf.mat','file')==2)
%         delete('filtered_enkf.mat');
%         display('Deleting the file filtered.mat')
%     end
%     save filtered_enkf Tstamp filtered_data_enkf
end