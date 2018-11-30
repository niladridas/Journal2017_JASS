
function out = mee2eci(input)
% input can be matrix with rows denoting one sample and column one variable
out = zeros(size(input,1),6);
for i = 1:size(input,1)
    out(i,:) = coe2eci(mee2coe(input(i,:)));
end
end