function A = assembleMatrixProb_Cont(kappa, mesh)

A = sparse(mesh.N, mesh.N);
A(1, 1) = -1 / (mesh.xInt(1) - mesh.xMin) * kappa(2) ...
          -1 / (mesh.xInt(2) - mesh.xInt(1)) * kappa(4);
A(1, 2) = 1 / (mesh.xInt(2) - mesh.xInt(1)) * kappa(4);

mesh.xInt = [mesh.xInt, mesh.xMax];
for i = 2 : mesh.N-1
    A(i, i) = -1 / (mesh.xInt(i) - mesh.xInt(i-1)) * kappa(2*i) ...
              -1 / (mesh.xInt(i+1) - mesh.xInt(i)) * kappa(2*(i+1));
    A(i, i-1) = 1 / (mesh.xInt(i) - mesh.xInt(i-1)) * kappa(2*i);
    A(i, i+1) = 1 / (mesh.xInt(i+1) - mesh.xInt(i)) * kappa(2*(i+1));
end
A(end) = -1 / mesh.h;
A = -A;
end