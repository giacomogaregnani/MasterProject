function uPoint = pointEval(u, x, mesh)

uPoint = interp1(mesh.x, u, x)';