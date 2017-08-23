function A = assembleMatrix_Cont(kappa, mesh)

A = sparse(mesh.N, mesh.N);
A(1, 1) = - kappa(2) - kappa(4);
A(1, 2) = kappa(4);
for i = 2 : mesh.N-1
    A(i, i) = -kappa(2*i) - kappa(2*(i+1));
    A(i, i-1) = kappa(2*i);
    A(i, i+1) = kappa(2*(i+1));
end
A(end) = -1;
A = -A / mesh.h;