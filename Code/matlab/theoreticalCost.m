clc
clear
close all
% desired accuracy
epsilon = 10 .^ [-0.001 : -1: -5]';
% order of the method
q = 1 : 2 : 5;
% number of levels
L = zeros(length(epsilon),length(q));
cost = L;
for i = 1 : length(q)
    L(:, i) = abs(log2(epsilon.^(1/q(i))));
    cost(:, i) = L(:, i) .* epsilon.^(-2);
end
markers = 'o*x><';
colors = 'brkgc';
legenditems = {};
for i = 1 : length(q)
    loglog(epsilon, cost(:, i), colors(i), 'Marker', markers(i))
    legenditems{i} = ['q = ' num2str(q(i)) ', MLMC'];
    hold on
end

% Cost MC
for i = 1 : length(q)
   loglog(epsilon, epsilon.^(-2 - 1/q(i)), ['--' colors(i)], 'Marker', markers(i));
   legenditems{i + length(q)} = ['q = ' num2str(q(i)) ', MC'];
end
set(gca, 'xdir','reverse')
legend(legenditems, 'Location', 'NW')
xlabel('\epl')
ylabel('cost')


