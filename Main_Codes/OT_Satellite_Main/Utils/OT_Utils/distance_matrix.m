% This function gives the distance values in form of a matrix D
function D = distance_matrix(X_f)
% The are M number of samples each with dimension 2
M = size(X_f,2);
D = zeros(M,M);
for i = 1:M-1
    for j = i+1:M
        % Modified Cost
        % Specific to the Satellite Problem
        % The first 5 variables are defined on the real space
        % The last element is a circular variable
        D1 = norm(X_f(1:5,i)-X_f(1:5,j));
        D2 = norm(min(abs(X_f(6,i)-X_f(6,j)),(2*pi-abs(X_f(6,i)-X_f(6,j)))));
        D(i,j) = (D1^2+D2^2)^(0.5);
%         
%          D(i,j)= norm(X_f(1:5,i)-X_f(1:5,j))+...
%              norm(min(abs(X_f(6,i)-X_f(6,j)),(2*pi-abs(X_f(6,i)-X_f(6,j)))));
         D(j,i) = D(i,j);
    end
end
end