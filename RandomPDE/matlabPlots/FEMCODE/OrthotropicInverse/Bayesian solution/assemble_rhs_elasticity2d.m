function rhs = assemble_rhs_elasticity2d(mesh, femspace, vp)
%ASSEMBLE_RHS Assembles RIGHT HAND SIDE

if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end

%% CONSTANTS
dim = size(mesh.elem,2) - 1;
ND = femspace.ndof;
LD = femspace.ldof;
lambda =  vp.flambda;
weight = vp.fweight;
NQ = numel(weight);

%% ASSEMBLING
rhs=zeros([ND, vp.fdim, vp.NF]);
for j=1:LD
    arhs = 0;
    for i=1:NQ
        %% L2 PART OF FORCE FIELD
        if isnumeric(vp.f)
            if ~vp.fully_discrete_rhs
                eval_fun = vp.f;
            else
                eval_fun = vp.f(:,:,:,i);
            end
        else
            eval_fun = evalf(mesh, 'all', [], lambda(i,:), vp.f);
        end
        arhs = arhs + weight(i) * ...
            bsxfun(@times, eval_fun(:,1:vp.fdim,:), ...
            evalb(mesh,'all',lambda(i,:), j, 0, femspace.elemtype));
        
        %% H-1 PART OF FORCE FIELD (only for Poisson)
        if vp.derivatives
            if vp.fdim > 1 && strcmp(vp.type,'micro_linear_elasticity_orth')
                if ~isnumeric(vp.a)
                    xloc = get_rc(mesh, [], 'all', [], lambda(i,:));
                    vpa11 = vp.a(xloc,1,1); vpa21 = vp.a(xloc,2,1);
                    vpa12 = vp.a(xloc,1,2); vpa22 = vp.a(xloc,2,2); vpa33 = vp.a(xloc,3,3);
                else
                    vpa11 = vp.a(:,1,1,i); vpa21 = vp.a(:,2,1,i);
                    vpa12 = vp.a(:,1,2,i); vpa22 = vp.a(:,2,2,i); vpa33 = vp.a(:,3,3,i);
                end
                arhs(:,1,:) =  arhs(:,1,:) + weight(i) * bsxfun(@times, vpa11, ...
                    bsxfun(@times, eval_fun(:,3,:), ...
                    evalb(mesh, 'all', lambda(i,:), j, 1, ...
                    femspace.elemtype)));
                arhs(:,2,:) =  arhs(:,2,:) + weight(i) * bsxfun(@times, vpa12, ...
                    bsxfun(@times, eval_fun(:,3,:), ...
                    evalb(mesh, 'all', lambda(i,:), j, 2, ...
                    femspace.elemtype)));
                arhs(:,1,:) =  arhs(:,1,:) + weight(i) * bsxfun(@times, vpa21, ...
                    bsxfun(@times, eval_fun(:,4,:), ...
                    evalb(mesh, 'all', lambda(i,:), j, 1, ...
                    femspace.elemtype)));
                arhs(:,2,:) =  arhs(:,2,:) + weight(i) * bsxfun(@times, vpa22, ...
                    bsxfun(@times, eval_fun(:,4,:), ...
                    evalb(mesh, 'all', lambda(i,:), j, 2, ...
                    femspace.elemtype)));
                arhs(:,1,:) =  arhs(:,1,:) + weight(i) * bsxfun(@times, vpa33, ...
                    bsxfun(@times, eval_fun(:,5,:), ...
                    evalb(mesh, 'all', lambda(i,:), j, 2, ...
                    femspace.elemtype)));
                arhs(:,2,:) =  arhs(:,2,:) + weight(i) * bsxfun(@times, vpa33, ...
                    bsxfun(@times, eval_fun(:,5,:), ...
                    evalb(mesh, 'all', lambda(i,:), j, 1, ...
                    femspace.elemtype)));
            elseif vp.fdim == 1
                for m=1:dim
                    for k=1:dim
                        if vp.fully_discrete
                            vpa = vp.a(:,m,k,i);
                        elseif isnumeric(vp.a) % numel(vp.a) == dim^2
                            vpa = vp.a(m,k);
                        else
                            xloc = get_rc(mesh, [],'all',[],lambda(i,:));
                            vpa = vp.a(xloc,m,k);
                        end
                        arhs = arhs + weight(i) * bsxfun(@times, vpa, ...
                            bsxfun(@times, eval_fun(:,k+1,:), ...
                            evalb(mesh, 'all', lambda(i,:), j, m, ...
                            femspace.elemtype)));
                    end
                end
            end
        end
    end
    arhs = bsxfun(@times, arhs, mesh.volume);
    dof = get_dof(mesh, 'all', femspace, j);
    for i=1:vp.fdim*vp.NF
        rhs(:,i) = rhs(:,i) + accumarray(dof, arhs(:,i), [ND, 1]);
    end
%     for i = 1:vp.fdim
%         for k = 1:vp.NF
%             rhs(:,i,k) = rhs(:,i,k) + accumarray(dof, arhs(:,i,k), [ND, 1]);
%         end
%     end
end
end
