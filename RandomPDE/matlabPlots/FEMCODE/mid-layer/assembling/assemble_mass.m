function M = assemble_mass(mesh, femspace)
%ASSEMBLE_MASS Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
NT=size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
mquad = 2*femspace.deg;
[lambda, weight] = quadpts(dim, mquad);
NQ = numel(weight);


%% INIT
mesh.volume = simplex_volume(mesh);

%% MATRIX GENERATION
M=sparse(ND, ND);
for i=1:LD
for j=1:LD
    if j<i
        continue;
    end
    Aij = zeros(NT,1);
    for m=1:NQ
        Aij = Aij + weight(m) .* ...
            evalb(mesh, 'all', lambda(m,:), i, 0, femspace.elemtype) .* ...
            evalb(mesh, 'all', lambda(m,:), j, 0, femspace.elemtype);
    end
    Aij = Aij .* mesh.volume;
    M = M + sparse(...
        double(get_dof(mesh, 'all', femspace, i)), ...
        double(get_dof(mesh, 'all', femspace, j)), ...
        Aij, ND, ND);
    if i<j
        M = M + sparse(...
            double(get_dof(mesh, 'all', femspace, j)), ...
            double(get_dof(mesh, 'all', femspace, i)), ...
            Aij, ND, ND);
    end
end
end
end

