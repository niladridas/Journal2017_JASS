clear all;clc;
rng(1,'twister');
s = {};
for i = 1:200
    s{i} = rng;
    randi(1);
end
save('./data/s.mat','s');