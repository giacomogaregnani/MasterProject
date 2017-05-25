function plotBoxCollection(fileName)

V = dlmread([fileName, '.txt']);
nCubes = size(V, 1);

figure
hold on
if size(V, 2) == 4
    
    for i = 1 : size(V, 1)
        X = V(i, 1:2);
        R = V(i, 3:4);
        fill([X(1); X(1)+R(1); X(1)+R(1); X(1)     ], ...
            [X(2); X(2);      X(2)+R(2); X(2)+R(2)], ...
            'red');
    end
else
    for i = 1 : nCubes
        drawCube(V(i, 1:3), V(i, 4:6), 'red');
    end
    view(57, 15)
end

return