function P = Optimal_Transport(X_f,Y,cost,Aeq,Aeq_1,weight)
M = size(X_f,2);
D = cost(X_f);
W = weight(X_f,Y);
D1 = reshape(D,[M*M,1]);

% x0 = 0.01*ones(M*M,1);
% fun = @(T) dot(D1,T); % Creating ghost function

beq = (1/M)*ones(M,1);
% opts = optimset('Display','iter','Algorithm','sqp',...
%                   'MaxFunEval',inf,'MaxIter',Inf);
% x = fmincon(fun,x0,[],[],[Aeq;Aeq_1],[beq;W'],zeros(M*M,1),[],[],opts);

A = [Aeq;Aeq_1];
% options = optimoptions('linprog','Algorithm','dual-simplex');
options = optimset('display', 'off');
x = linprog(D1,[],[],A,[beq;W'],zeros(M*M,1),[],[],options);

P = M.*x;
P = reshape(P,[M,M]);
% keyboard
% P = zeros(1,1);
end