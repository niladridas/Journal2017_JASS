% Parsing function module
% Input: TLE file
% Output: date_mee mat file and also the date_mee variable containing the 
% time when the data was taken and
% the corresponding MEE elements measured
% Options: 
% 1. If no option argument is passed along with fileID just date_mee is
% formed
% 2. If '1' is passed along with fileID then the initial MEE value is
% propagated using the dynamics and the error plots between the propagated
% values and the observed values are shown
% 3. If '2' is passed along with fileID then taking each observation point
% and propagating it to the immediately next time point, we compare the
% value with that of the observed value.
function [date_mee,time_pt_days] = tle_parse(varargin)
    FileID = varargin{1};
    FLAG = 1; % To keep on reading
    % Constant Elements:
    sgp4_const_elements = struct('XKMPER',6378.135,'a_E',1.0,...
        'k_2',5.413080e-4,'k_e',0.74366916e-1,...
        'A_3_0',-0.0000025323,'s',1.01222932,...
        'qo_s4',((120-78)*1.0/6378.135)^4,'k_4',0.62098875e-6);
    date_pos_vel = [];
    while FileID~=-1
        [date,dot_n,ddot_n,bstar,inclination,rightascension,...
            essentricity,argperi,meanano,meanmotion,...
            totalrev,FileID] = tleparser(FileID,FLAG);
        if FileID == -1;break;end
        tle_parsed_elements = struct('n_o',meanmotion,...
            'e_o',essentricity,'i_o',inclination,'M_o',meanano,...
        'w_o',argperi,'Ohm_o',rightascension,'B_star',bstar,'t_o',0,'t',0);
        % pos_vel : units in Km and secs
        date_pos_vel = [date_pos_vel;[date,sgp4(tle_parsed_elements,sgp4_const_elements)]];
    end % Data reading finished
%     save('date_pos_vel.mat','date_pos_vel')
    % ECI to MEE elements
    mu = 398600.4418;
    date_coe = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
    date_coe(:,1:6) = date_pos_vel(:,1:6);
    date_mee = zeros(size(date_pos_vel,1),size(date_pos_vel,2));
    date_mee(:,1:6) = date_pos_vel(:,1:6);

    for i = 1:size(date_pos_vel,1)
        % use eci to coe
        date_coe(i,7:end) = eci2coe(date_pos_vel(i,7:12),mu);
        % Now convert from coe to mee
        date_mee(i,7:end) = coe2mee(date_coe(i,7:end));
    end
    % The LAST COLUMN is an angle
    date_mee(:,end) = mod(date_mee(:,end),2*pi);
    time_pt_days = [];
    time_pt_days(1) = 0;
    for i = 1:(size(date_mee,1)-1)
       Tsince = time_diff(date_mee(i,1:6),date_mee((i+1),1:6));
       time_pt_days(i+1) = time_pt_days(i)+(Tsince/(60*60*24));
    end
    % Some TLE files are up-down flipped. We expect the entried to be past
    % to the present. But some donot follow this.
    if time_pt_days(2) < 0||time_pt_days(3)<0
        date_mee = flipud(date_mee);
    end
    save('date_mee.mat','date_mee');
    save('time_pt_days.mat','time_pt_days'); %TODO 
    % % Test how accurate the equinoctial element based dynamics is
    % % Take on mee element, propagate it to the next time point and then find
    % % the error between the propagated and measured
    % % If there are N time steps the loop is for N-1
    if varargin{2} == '1'
        R_e = 6378.135;
        J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
        prop_params = struct('mu',mu,'R_e',R_e,'J',J);
        mee_predict = zeros(size(date_mee,1), 6);
        mee_predict(1,:) = date_mee(1,7:12);
        for i = 1:(size(date_mee,1)-1)
           % Call function to return the time difference in seconds
           Tsince = time_diff(date_mee(i,1:6),date_mee((i+1),1:6));
           [t,mee_set] = ode45(@(t,mee) equinoctial_dyn(t,mee,prop_params),[0 Tsince],mee_predict(i,:));
           mee_predict(i+1,:) = mee_set(end,:);
        end
        mee_ob = date_mee(:,7:end);
        mee_predict(:,end) = mod(mee_predict(:,end),2*pi);
        error_e = mee_ob - mee_predict;
        figure
        ylabel_names = ['p','f','g','h','k','L'];
        for i = 1:6
            subplot(2,3,i)
            a = plot(time_pt_days,error_e(:,i));
            a.Color = 'red';
            a.LineStyle = '-';
            a.LineWidth = 1;
            a.Marker = 'o';
            a.MarkerSize = 1.5;
            a.MarkerEdgeColor = 'red';
            a.MarkerFaceColor = 'red';
            xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
            ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
        end
        title_str = strcat('Error plot of propagated MEE and observed MEE when simulated from initial observations');
        axes( 'Position', [0, 0.95, 1, 0.05] ) ;
        set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
        text( 0.55, 0, title_str, 'FontSize', 15', 'FontWeight', 'Bold', ...
              'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    end
end