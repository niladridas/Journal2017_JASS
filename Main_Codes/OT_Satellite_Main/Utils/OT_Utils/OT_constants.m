
function [Aeq,Aeq_1] = OT_constants(M)
Aeq = [];
for i = 1:M
    Aeq(i,(i-1)*M+1:i*M) = ones(M,1);  
end

Aeq_1 = zeros(M,M*M);
temp_Aeq_1 = zeros(1,M*M);
for i = 1:M
    temp_Aeq_1(1,(i-1)*M+1)= 1;
end
for i = 1:M
    Aeq_1(i,:) = circshift(temp_Aeq_1,[0,i-1]);
end
end