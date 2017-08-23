function mesh = buildMesh(mesh)

mesh.h = (mesh.xMax - mesh.xMin) / mesh.N;
mesh.x = linspace(mesh.xMin, mesh.xMax, mesh.N+1);
mesh.xInt = mesh.x(2:end-1);