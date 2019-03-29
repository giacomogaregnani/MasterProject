function [r_1,r_2,r_3,r_4] = compute_distance(X)
% Compute DISTANCES from boundaries

r_1 = X(1) + 1;
r_2 = X(2) + 1;
r_3 = 1 - X(1);
r_4 = 1 -X(2);

end