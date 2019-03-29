function div = assemble_div(mesh, ufemspace, pfemspace, vp)
%ASSEMBLE_DIV Assembles DIV MATRIX (Stokes problem)

%% CONSTANTS
dim = size(mesh.elem,2) - 1;
NT=size(mesh.elem,1);
uND = ufemspace.ndof;
pND = pfemspace.ndof;
if ~isfield(vp,'divquad'), 
	vp.divquad = ufemspace.deg + pfemspace.deg - 1; 
end
[lambda, weight] = quadpts(dim, vp.divquad);

%% ASSEMBLE
div = sparse(pND, dim*uND);
for dd=1:dim
    for j=1:ufemspace.ldof
        for k=1:pfemspace.ldof
            Bij = zeros(NT,1);
            for i= 1 :numel(weight)
                Bij = Bij + weight(i) * ...
                    evalb(mesh, 'all', lambda(i,:), j, dd, ufemspace.elemtype) .* ...
                    evalb(mesh, 'all', lambda(i,:), k, 0,  pfemspace.elemtype);
            end
            Bij = Bij .* mesh.volume;
            div = div - sparse(... % SIGN NECESSARY!
                double(get_dof(mesh,'all',pfemspace,k)), ...
                double(get_dof(mesh,'all',ufemspace,j)) + (dd-1)*uND, ...
                Bij, pND, dim*uND);
        end
    end
end

end

