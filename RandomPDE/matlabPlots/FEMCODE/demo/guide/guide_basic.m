%% GUIDE - Matlab skills
% Ondrej Budac, 2014

%% 1. Arrays
%
% * used for structured data

a = 1:5;    % row vector
b = (1:5)'; % column vector
c = zeros(4,5); % 4 by 5 matrix
d = ones(2,3,4); % a 3D array
e = false(2,3,4,5); % a 4D logical array

dot(a,b); % works despite nonmatching sizes