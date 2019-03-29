function dist = ComputeDistance(X)
% Compute DISTANCE FROM WALLS

r1 = X(1) + 1;
r2 = X(2) + 1;
r3 = 1 - X(1);
r4 = 1 -X(2);

dist = min([r1, r2, r3, r4]);

end