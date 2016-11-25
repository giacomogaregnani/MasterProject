clear all

x = dlmread('jacHires.txt');
n = 8;
nTimesteps = size(x, 1) / 8;

eVal = zeros(nTimesteps, 8);
count = 1;
for i = 1 : nTimesteps
   J = x(count : count + 7, :);
   eVal(i, :) = eig(J)';
   count = count + 8;
end

plot(1 : nTimesteps, eVal);