function D = densDiff(x1, W1, x2, W2)

r = x1(1, 3:4);
x1 = [x1(:, 1:2)];
nx1 = size(x1, 1);
x2 = [x2(:, 1:2)];
nx2 = size(x2, 1);

corners = zeros(4, 2);
xMin = min(min(x1(:, 1)), min(x2(:, 1)));
xMax = max(max(x1(:, 1)), max(x2(:, 1))) + r(1);
yMin = min(min(x1(:, 2)), min(x2(:, 2)));
yMax = max(max(x1(:, 2)), max(x2(:, 2))) + r(2);
xGrid = xMin : r(1) : xMax;
yGrid = yMin : r(2) : yMax;

origin = [xMin, yMin];
GridDens1 = zeros(length(xGrid), length(yGrid));
GridDens2 = zeros(length(xGrid), length(yGrid));

for i = 1 : nx1
    indices = coordToInd(origin, r, x1(i, 1:2));
    GridDens1(indices) = W1(i);
end

for i = 1 : nx2
    indices = coordToInd(origin, r, x2(i, 1:2));
    GridDens2(indices) = W2(i);
end

D = sum(sum(abs(GridDens1 - GridDens2)));