% Author: Niladri Das
% Affiliation: UQ Lab, TAMU, Aerospace Engineering
% Date: 7 July 2017

% synthetic state data
% NORAD data of ISS
clc;clear;close all;

% DYNAMICS
% ODE parameters for dynamics
mu = 398600.4418;
R_e = 6378.135;
J = [0;0.00108263;-2.51e-6;-1.6e-6;-1.3e-7;5e-7];
prop_params = struct('mu',mu,'R_e',R_e,'J',J);
equinoc_dyn = @(t,x) equinoctial_dyn(t,x,prop_params);

load('date_mee.mat');
% 
% t1 = zeros(1,size(date_mee,1));
% t1(1) = 0;
% for i = 1:(size(date_mee,1)-1)
%     t1(i+1) = time_diff(date_mee(i,1:6),date_mee(i+1,1:6));
% end
% 
% % Modifing the t1 variable by adding some more time intervals at the end
% % This is done by concatenating the whole array (except the first zero term)
% t1 = [t1,t1(2:end),t1(2:end),t1(2:end),t1(2:end)];



t1 = linspace(0,3*92.75*60,10);




% MODE 1:
% Synthetic Data generation
% We generate synthetic actual states
syn_state = zeros(size(t1,2),6);
syn_state(1,:) = date_mee(1,7:12);
syn_state(1,6) = mod(syn_state(1,6),2*pi);
for i= 2:size(t1,2)
[~,x_temp] = ode45(equinoc_dyn,[0 t1(i)],syn_state(i-1,:)');
syn_state(i,:) = x_temp(end,:);
syn_state(i,6) = mod(syn_state(i,6),2*pi);
end
% syn_state(:,6) = mod(syn_state(:,6),2*pi);
% We generate synthetic observations

% Plot the variation of data

% t_spot = zeros(1,length(t1));
% for i = 1:length(t1)
%     t_spot(i) = sum(t1(1:i));
% end
% 
% % t_spot is in seconds
% % conversting it into days
% t_spot = t_spot./(60*60*24);

ylabel_names = ['p','f','g','h','k','L'];
figure(1)
for i = 1:6
    subplot(3,2,i);
    % Fix the color
    % Fix the fonts 
    % Fix the symbol on the lines like scatter plot
    % Give the title
    % Give the footer note
    
    %     set(gca,'Xtick',0:3:15);

    set(gca,'fontWeight','bold','Xtick',0:15);
    c = plot(t1./(60*60),syn_state(:,i)); 
    xlabel('Time in hours','FontName','sans-serif','FontSize',16,'FontWeight','Bold');
    ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',16,'FontWeight','Bold');
    c.Color = 'blue';
    c.LineStyle = '-';
    c.LineWidth = 2;
    c.Marker = 'o';
    c.MarkerSize = 4;
    c.MarkerEdgeColor = 'blue';
    c.MarkerFaceColor = 'blue';
    grid on
end