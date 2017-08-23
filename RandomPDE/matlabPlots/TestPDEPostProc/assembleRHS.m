function F = assembleRHS(f, mesh, RightBC)

F = mesh.h * f(mesh.xInt);

F = [F; RightBC / mesh.h];