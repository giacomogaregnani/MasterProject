function [ ASYM ] = assemble_interface_sym( mesh, femspace, vp )
%% Initialization
NA = femspace{1}.ndof;
dim = size(mesh{1}.elem,2) - 1;
ASYM = sparse(dim*NA, dim*NA);
    % quadrature formula on edges
    [lambda, weight] = quadpts(dim-1, 2*femspace{1}.deg);
    NQ = numel(weight);
    % for each js = 1:dim we consider those interface elements that have the
    % boundary edge oppoiste to the node with local coordinate js. This
    % simplifies how we obtain the degrees of freedom.
    for js=1:dim+1
      ind = mesh{1}.bdflag{3}(:,2) == js;
      % in this pass we only work with interface elements with numbers swhe
      swhe = mesh{1}.bdflag{3}(ind, 1);
      if isempty(swhe),
        continue;
      end
      % compute normals and tangents at the considered elements. Here we assume
      % that the dimension is 2. One of the few parts that will need
      % generalization before coming to 3D
      normals = get_normals(mesh{1}, swhe, js);
      tangents = [-normals(:,2), normals(:,1)]; % works only in 2D
      % modify the quadrature formula (just the barycentric coordinates) to
      % simplify the integration: we integrate over elements.
      lam = [lambda(:,1:js-1), zeros(size(lambda,1),1), lambda(:,js:dim)];
      % which local degrees of freedom are of interest (those that lie on the
      % interface edge).
      ldof = find(femspace{1}.nodelambda(:,js)==0)';
      % volumes (lengths in 2D, surfaces in 3D) of the interface edges.
      ssvol = subsimplex_volume(mesh{1}, [], swhe, [1:js-1,js+1:dim+1]);
      for k=1:NQ
        % coef is the part alpha * vp{2}.n / sqrt(tau * K * tau), computed at the
        % quadrature points. We assume that we are in 2 dimensions so that
        % there is just one tau.
        coef =0;
        % if the tensor is a function handle, we will need coordinates of the
        % quadrature nodes
        if ~isnumeric(vp{2}.a)
          xloc = get_rc(mesh{1},[],swhe,[],lam(k,:));
        end
        % now we assemble the values tau_{t1} * K_{t1, t2} * rau_{t2} and sum
        % them up
        for t1 = 1:dim
          for t2 = 1:dim
            if isnumeric(vp{2}.a)
              if numel(vp{2}.a) == 1
                if t1==t2,
                  vpa = vp{2}.a;
                else
                  continue;
                end
              else
                vpa = vp{2}.a(t1, t2);
                if vpa == 0
                  continue;
                end
              end
            else
              vpa = vp{2}.a(xloc,t1,t2);
            end
            coef = coef + vpa .* tangents(:,t1) .* tangents(:,t2);
          end
        end
        % the assembled coef is now modified to its final form, where
        % alpha is the Beavers-Saffman constant determined experimentally.
        coef = 1./sqrt(coef);
        % we continue by assembling the required matrix. We loop over local
        % degrees of freedom (of interest) twice and then assemble the values
        % into a sparse matrix.
        for m1=ldof
          % degrees of freedom (dof1) and basis functions (evb1) evaluated at
          % the quadrature points
          dof1 = get_dof(mesh{1}, swhe, femspace{1}, m1);
          evb1 = evalb(mesh{1}, swhe, lam(k,:), m1, 0, femspace{1}.elemtype);
          for m2=ldof
            % degrees of freedom (dof2) and basis functions (evb2) evaluated at
            % the quadrature points
            dof2 = get_dof(mesh{1}, swhe, femspace{1}, m2);
            evb2 = evalb(mesh{1}, swhe, lam(k,:), m2, 0, femspace{1}.elemtype);
            for t1 = 1:dim
              for t2 = 1:dim
                Aij = evb1 .* tangents(:,t1) .* evb2 .* tangents(:,t2) .* ...
                  weight(k) .* ssvol .* coef;
                ASYM = ASYM + sparse(double(dof1+(t1-1)*NA), double(dof2+(t2-1)*NA), Aij, dim*NA, dim*NA);
              end
            end
          end
        end
      end
    end
  end