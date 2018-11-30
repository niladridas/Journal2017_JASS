% Author: Niladri Das
% Affiliation: UQLab, Aerospace Engineering, TAMU
% Date: 14 May 2017
% 




k = linspace(0,1000,100);
plot(k,fun1(k,0))


% plot(k,1./value)
% 
% function f = fun1(k,x)
% % This function returns the vales evaluated at the given k
% f1 = exp(k*cos(x));
% f2 = besselj(0,k);
% f3 = (k+1).^2;
% f = (f1.*f3)./(f2);
% end
% 
function f = fun1(k,x)
% This function returns the vales evaluated at the given k
f1 = exp(k*cos(x));
f2 = besseli(0,k);
f = f1./(f2);
end