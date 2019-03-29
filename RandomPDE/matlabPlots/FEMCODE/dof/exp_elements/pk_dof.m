function femspace = pk_dof( mesh, k )
%PK_DOF Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
N  = size(mesh.node,1);
NT = size(mesh.elem,1);
M = min(d, k-1);

% vertices DOFs
femspace.ndof = N;
femspace.dofs(1) = 1;

ind = 0;
for i=1:M
    V = factorial(dim+1) / factorial(i+1) / factorial(dim-i);
    W = factorial(k-1) / factorial(i) / factorial(k-1-i);
    [elem2sim, sim] = auxstructure(mesh, ...
        sprintf('elem2sim%d',i), sprintf('sim%d',i));
    
    femspace.elem2dof(:,ind+1: ind+V*W) = [N + double(elem2edge), ...
        N + NE + double(elem2edge), ...
        N + 2*NE + double(elem2face)];
    

    if (M==dim)
        femspace.dofs(M+1) = femspace.ndof;
        femspace.ndof = femspace.ndof + W*NT;
    end
    
        
end


end

