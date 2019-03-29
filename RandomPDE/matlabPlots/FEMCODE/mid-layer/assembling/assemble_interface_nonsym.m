function ANSYM = assemble_interface_nonsym( mesh, femspace, interface )

%% Initialization
NA = femspace{1}.ndof;
dim = size(mesh{1}.elem,2) - 1;
NAD = femspace{2}.ndof;

ANSYM = sparse(dim*NA, NAD);
% quadrature formula on edges
[lambda, weight] = quadpts(dim-1, femspace{1}.deg + femspace{2}.deg);
NQ = numel(weight);
dof = cell(2,1);
% We go through interface{1} and interface{2}, one at a time.
for r = 1:2
  % for each js = 1:dim we consider those interface elements in mesh{r}
  % that have the boundary edge oppoiste to the node with local
  % coordinate js. This simplifies how we obtain the degrees of freedom
  for js=1:dim+1
    ind = mesh{r}.bdflag{3}(:,2) == js;
    % in this pass we only work with interface elements with numbers
    % whe1 (in mesh{r})
    whe1 = mesh{r}.bdflag{3}(ind, 1);
    % we restrict to those elements that are present in
    % interface{r}(:,1)
    [~,i1,i2] = intersect(whe1, interface{r}(:,1));
    whe1 = whe1(i1);
    % corresponding interface elements from mesh{3-r} have element
    % numbers whe2
    whe2 = interface{r}(i2,2);
    if isempty(whe1), continue; end
    % compute normals
    normals = (-1)^(r+1)*get_normals(mesh{r}, whe1, js); % TODO: watch out for direction!
    %Antoine: 27/3, changed sign
    % modify the quadrature formula (just the barycentric coordinates) to
    % simplify the integration: we integrate over elements.
    lam1 = [lambda(:,1:js-1), zeros(size(lambda,1),1), lambda(:,js:dim)];
    % which local degrees of freedom are of interest (those that lie on the
    % interface edge).
    ldof1 = find(femspace{r}.nodelambda(:,js)==0)';
    ldof2 = 1:femspace{3-r}.ldof;
    % volumes (lengths in 2D, surfaces in 3D) of the interface edges.
    ssvol = subsimplex_volume(mesh{r}, [], whe1, [1:js-1,js+1:dim+1]);
    for k=1:NQ
      % cartesian coordinates of quadrature nodes
      xloc = get_rc(mesh{r},[],whe1,[],lam1(k,:));
      % barycentric coordinates for integration in mesh{3-r}
      lam2 = get_lc2(mesh{3-r}, whe2, xloc);
      % We loop over local
      % degrees of freedom (of interest) in mesh{r}
      for m1=ldof1
        % degrees of freedom (dof1) and basis functions (evb1) evaluated at
        % the quadrature points in mesh{r}
        dof{r} = get_dof(mesh{r}, whe1, femspace{r}, m1);
        evb1 = evalb(mesh{r}, whe1, lam1(k,:), m1, 0, femspace{r}.elemtype);
        % We loop over local
        % degrees of freedom (of interest) in mesh{3-r}
        for m2=ldof2
          % degrees of freedom (dof2) and basis functions (evb2) evaluated at
          % the quadrature points in mesh{3-r}
          dof{3-r} = get_dof(mesh{3-r}, whe2, femspace{3-r}, m2);
          evb2 = evalb(mesh{3-r}, whe2, lam2, m2, 0, femspace{3-r}.elemtype);
          for t = 1:dim
            Aij = evb1 .* evb2 .* normals(:,t) .* weight(k) .* ssvol;
            ANSYM = ANSYM + sparse(double(dof{1}+(t-1)*NA), double(dof{2}), Aij, dim*NA, NAD);
          end
        end
      end
    end
  end
end
end

