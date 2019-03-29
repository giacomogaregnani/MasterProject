%% GUIDE 4 - Micro problems and correctors
% Ondrej Budac, Timothee Pouchon, 2014

%% Definitions
d=2; dim = 2;
a = @guide4_a1;
dbtype guide4_a1.m

%% Mesh
mesh = structured_mesh(repmat([0,1],[1,d]), 200, struct('periodic',true)); 
mesh.bdflag = 'dirichlet'; % applies to all the boundary = zero boundary
mesh.volume = simplex_volume(mesh);
NT = size(mesh.elem,1);

%% Cell problems XI
rhs = permute([zeros(d,1), -eye(d)],[3,2,1]);
vp = struct('a', a, 'f', rhs, 'elemtype', 'p1');
[xi, femspace, mesh, vp] = poisson(mesh, vp);

%%

% Test solvability
display(sum(assemble_rhs(mesh, femspace, vp),1));

%%

for i=1:d
    figure; simpplot_sol(mesh, xi(:,:,i)); title(['\xi_' num2str(i)]); snapnow;
end

%% Compute Homogenized tensor
a0 = zeros(d,d);
[lambda, weight] = quadpts(dim, femspace.deg);
for i=1:d
  for j=1:d
    res = 0;
    for s=1:numel(weight)
      res = res + weight(s) * ...
        evalf(mesh, 'all', femspace, lambda(s,:), a, [], i, j);
      for m = 1:dim
        res = res + weight(s) * ...
          evalf(mesh, 'all', femspace, lambda(s,:), a, [], i, m) .* ...
          evalf(mesh, 'all', femspace, lambda(s,:), xi(:,:,j), m);
      end
    end
    a0(i,j) = dot(res, mesh.volume);
  end
end
display(a0);

%% First correctors: THETA
vp = struct('a', a, 'elemtype', 'p1');
[lambda, weight] = quadpts(dim, 1);
NQ = numel(weight);
rhs = zeros(NT, d+1, d^2, NQ);
enctheta = reshape(1:d^2, [d,d]);
for i=1:d
    for j=1:d        
        for s=1:numel(weight)
            rhs(:, 1, enctheta(i,j), s) = ...
                evalf(mesh, 'all', femspace, lambda(s,:), a, [], i, j) ...
                -a0(i, j);
            for m=1:d
                rhs(:, 1, enctheta(i,j), s) = rhs(:, 1, enctheta(i,j), s) + ...
                    evalf(mesh, 'all', femspace, lambda(s,:), a, [], i, m) .* ...
                    evalf(mesh, 'all', femspace, lambda(s,:), xi(:,:,j), m);
            end
            rhs(:, i+1, enctheta(i,j), s) = ...
                evalf(mesh, 'all', femspace, lambda(s,:), xi(:,:,j));
        end        
    end
end
vp.f = rhs;
[theta, femspace, mesh, vp] = poisson(mesh, vp);

%% 

% Test solvability
display(sum(assemble_rhs(mesh, femspace, vp),1));

%% 

% plot
for i=1:d
    for j=1:d
        clf; simpplot_sol(mesh, theta(:,:,enctheta(i,j))); 
        title(['\xi_{' num2str(i) num2str(j) '}']); snapnow;
    end
end

%% Second correctors: Kappa
vp = struct('a', a, 'elemtype', 'p1');
[lambda, weight] = quadpts(dim, 1);
NQ = numel(weight);
rhs = zeros(NT, d+1, d^3, NQ);
enckappa = reshape(1:d^3, [d,d,d]);
for i=1:d
    for j=1:d
        for k=1:d
            for s=1:numel(weight)
                rhs(:, 1, enckappa(i,j,k), s) = ...
                    (evalf(mesh, 'all', femspace, lambda(s,:), a, [], i, j) - a0(i, j)) .*...
                    evalf(mesh, 'all', femspace, lambda(s,:), xi(:,:,k));
                for m=1:d
                    rhs(:, 1, enckappa(i,j,k), s) = rhs(:, 1, enckappa(i,j,k), s) + ...
                        evalf(mesh, 'all', femspace, lambda(s,:), a, [], i, m) .* ...
                        evalf(mesh, 'all', femspace, lambda(s,:), theta(:,:,enctheta(j,k)), m);
                end
                rhs(:, i+1, enckappa(i,j,k), s) = ...
                    evalf(mesh, 'all', femspace, lambda(s,:), theta(:,:,enctheta(j,k)));
            end
        end
    end
end
vp.f = rhs;
[kappa, femspace, mesh, vp] = poisson(mesh, vp);

%% 

% Test solvability
display(sum(assemble_rhs(mesh, femspace, vp),1));

%%

% plot
for i=1:d
    for j=1:d
        for k=1:d
            clf; simpplot_sol(mesh, kappa(:,:,enckappa(i,j,k)));
            title(['\kappa_{' num2str(i) num2str(j) num2str(k) '}']); snapnow;
        end
    end
end