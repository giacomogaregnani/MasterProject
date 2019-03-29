function [usol, ufemspace] = reconstruct_velocity(mesh, vp, psol, pfemspace)

dim = size(mesh.elem, 2) - 1;
lambda = vp.alambda;
weight = vp.aweight;
NQ = numel(weight);

%% ASSEMBLE ALL
g = []; % INIT
for k=1:NQ
    if ~isnumeric(vp.a)
        xloc = get_rc(mesh,[],'all',[],lambda(k,:));
    end
    
    f = evalf(mesh, 'all', pfemspace, lambda(k,:), vp.f);
    gg = 0; % INIT
    for i=1:dim
        ggg = []; % INIT
        for j=1:dim
            if vp.fully_discrete % fully discrete tensor
                vpa = vp.a(:,i,j,k);
            elseif isnumeric(vp.a)
                vpa = vp.a(i,j);
            else
                vpa = vp.a(xloc,i,j);
            end
            gggadd = vpa .* ...
                (evalf(mesh, 'all', pfemspace, lambda(k,:), psol, i) - f(:,i+1,:));
            ggg = cat(2,ggg,gggadd);
        end
        gg = gg + ggg;
    end
    g = cat(4, g, gg);
end
g = -g; % INVERSE
clear gg ggg gggadd;

%% RESTORE VALUES FROM INTEGRATION POINTS
[ufemspace, usol] = restore(mesh, g, vp.elemtype);
end

