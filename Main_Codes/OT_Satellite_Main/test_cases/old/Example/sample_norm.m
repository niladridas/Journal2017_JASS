% Calculate norm for each sample size
clc;clear;close all;
load RMSE_enkf;
load RMSE_ot;

norm_ot = zeros(size(RMSE_enkf,1),1);
norm_enkf = norm_ot;
for i=1:size(RMSE_enkf,1)
   norm_ot(i,1) = norm(RMSE_ot(i,1:6));
   norm_enkf(i,1) = norm(RMSE_enkf(i,1:6));
end

a = plot(30:10:80,norm_ot,'r-o');
hold on
b = plot(30:10:80,norm_enkf,'b-o');
grid on
hold on
a.LineWidth = 2;
b.LineWidth = 2;
a.MarkerSize = 4;
b.MarkerSize = 4;
a.MarkerEdgeColor = 'red';
a.MarkerFaceColor = 'red';
b.MarkerEdgeColor = 'blue';
b.MarkerFaceColor = 'blue';
xlabel('Samples','FontName','sans-serif','FontSize',20,'FontWeight','Bold');
ylabel('Norm of RMSE','FontName','Times New Roman','FontSize',20,'FontWeight','Bold');
legend('OT','ENKF')
title('Variation of Norm of RMSE with samples for OT and EnKF');
for j = 1:size(RMSE_enkf,1)
    text(30+10*(j-1),norm_ot(j,1),num2str(RMSE_ot(j,7)));
    hold on
    text(30+10*(j-1),norm_enkf(j,1),num2str(RMSE_enkf(j,7)));
    hold on
end