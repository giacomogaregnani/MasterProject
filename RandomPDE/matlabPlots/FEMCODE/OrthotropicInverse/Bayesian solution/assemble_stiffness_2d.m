function stiffness = assemble_stiffness_2d(mesh, femspace, vp)

dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.alambda;
weight = vp.aweight;
NQ = numel(weight);

if ~isnumeric(vp.a)
	xloc = cell([1, NQ]);
	for m=1:numel(weight)
		xloc{m} = get_rc(mesh,[],'all',[],lambda(m,:));
	end
end

stiffness = sparse(2*ND, 2*ND);
for i = 1:LD
    for j = 1:LD
        Bij11 = zeros(NT,1);
        Bij12 = zeros(NT,1);
        Bij21 = zeros(NT,1);
        Bij22 = zeros(NT,1);
        for m = 1:NQ
            if ~isnumeric(vp.a)
                vpa11 = vp.a(xloc{m},1,1); vpa12 = vp.a(xloc{m},1,2); vpa21 = vp.a(xloc{m},2,1); vpa22 = vp.a(xloc{m},2,2);
                vpa33 = vp.a(xloc{m},3,3);
            else
                vpa11 = vp.a(:,1,1,m); vpa12 = vp.a(:,1,2,m); vpa21 = vp.a(:,2,1,m); vpa22 = vp.a(:,2,2,m);
                vpa33 = vp.a(:,3,3,m);
            end
            Bij11 = Bij11 + (weight(m)*vpa11).*evalb(mesh, 'all', lambda(m,:), i, 1, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 1, femspace.elemtype) + ...
                (weight(m)*vpa33).*evalb(mesh, 'all', lambda(m,:), i, 2, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 2, femspace.elemtype);
            Bij12 = Bij12 + (weight(m)*vpa12).*evalb(mesh, 'all', lambda(m,:), i, 1, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 2, femspace.elemtype) + ...
                (weight(m)*vpa33).*evalb(mesh, 'all', lambda(m,:), i, 2, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 1, femspace.elemtype);
            Bij21 = Bij21 + (weight(m)*vpa21).*evalb(mesh, 'all', lambda(m,:), i, 2, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 1, femspace.elemtype) + ...
                (weight(m)*vpa33).*evalb(mesh, 'all', lambda(m,:), i, 1, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 2, femspace.elemtype);
            Bij22 = Bij22 + (weight(m)*vpa22).*evalb(mesh, 'all', lambda(m,:), i, 2, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 2, femspace.elemtype) + ...
                (weight(m)*vpa33).*evalb(mesh, 'all', lambda(m,:), i, 1, femspace.elemtype).*...
                evalb(mesh, 'all', lambda(m,:), j, 1, femspace.elemtype);
        end
        Bij11 = Bij11.*mesh.volume;
        Bij12 = Bij12.*mesh.volume;
        Bij21 = Bij21.*mesh.volume;
        Bij22 = Bij22.*mesh.volume;
        stiffness = stiffness + sparse(double(get_dof(mesh, 'all', femspace, i)),...
            double(get_dof(mesh, 'all', femspace, j)), Bij11, 2*ND, 2*ND);
        stiffness = stiffness + sparse(double(get_dof(mesh, 'all', femspace, i)),...
            double(get_dof(mesh, 'all', femspace, j))+ND, Bij12, 2*ND, 2*ND);
        stiffness = stiffness + sparse(double(get_dof(mesh, 'all', femspace, i))+ND,...
            double(get_dof(mesh, 'all', femspace, j)), Bij21, 2*ND, 2*ND);
        stiffness = stiffness + sparse(double(get_dof(mesh, 'all', femspace, i))+ND,...
            double(get_dof(mesh, 'all', femspace, j))+ND, Bij22, 2*ND, 2*ND);
    end
end