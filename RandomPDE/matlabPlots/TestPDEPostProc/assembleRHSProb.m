function F = assembleRHSProb(f, mesh, RightBC)

F = zeros(mesh.N - 1, 1);

F(1) = (mesh.xInt(2) - mesh.xMin) * f(mesh.xInt(1)) / 2;
for i = 2 : mesh.N - 2
   F(i) = (mesh.xInt(i+1) - mesh.xInt(i-1)) * f(mesh.xInt(i)) / 2;
end
F(end) = (mesh.xMax - mesh.xInt(end-1)) * f(mesh.xInt(end)) / 2;

F = [F; RightBC / mesh.h];