function errNum = L2Error(mesh, uRef, u)

N = size(mesh.elements, 2);
f = @(x, y) reshape((uRef.evaluate(x, y) - u.evaluate(x, y)).^2, size(x,1), size(x,2));
errNum = zeros(N, 1);

for i = 1 : N
    vert = mesh.vertices(:, mesh.elements(1:3, i))';
    errNum(i) = sqrt(intm2(f, vert));
end