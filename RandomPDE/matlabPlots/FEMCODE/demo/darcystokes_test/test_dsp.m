%% GUIDE 5 - Darcy-Stokes coupling
% Ondrej Budac, Antoine Imboden, 2014
clear mesh interface vp
mesh = cell(2,1);
mesh{1} = structured_mesh([0,2,0,1],[2,1],struct('centre',true));
mesh{2} = structured_mesh([0,2,-1,0],[2,1],struct('centre',true));


mesh{1}.bdflag{1} = get_bdflag(mesh{1}, @(x)(x(:,2) == 1)); % dirichlet
mesh{1}.bdflag{2} = get_bdflag(mesh{1}, @(x)(x(:,1).*(2-x(:,1)) == 0)); % Neumann
mesh{1}.bdflag{3} = get_bdflag(mesh{1}, @(x)(x(:,1).*(2-x(:,1)).*(1-x(:,2)) ~= 0)); % interface

mesh{2}.bdflag{1} = zeros(0,2); % dirichlet
mesh{2}.bdflag{2} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)).*(x(:,2)+1) == 0)); % Neumann
mesh{2}.bdflag{3} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)).*(1+x(:,2)) ~= 0)); % interface


vp{1}.a = 1;
vp{1}.f = [1,0];
vp{1}.elemtype = 'p2';
vp{1}.pelemtype = 'p1';
vp{1}.bc = {[0,0], [0,0]};

vp{2}.a = eye(2);
vp{2}.f = @(x)(repmat([0, 0, 0], size(x,1),1));
vp{2}.elemtype = 'p1';
vp{2}.bc = 0;
%Volumetric porosity qdarcy = n*Vdarcy, where qdarcy is the discharge vector
vp{2}.n = 1;
%Experimental constant, correlated with the average pore size
vp{2}.alpha = 1;
% Gravitational constant
vp{2}.g = 1;
% Fluid volumic density
vp{2}.rho = 1;

% for i=1:3 % increase to get a finer mesh
%   for r=1:2
%     mesh{r} = uniformrefine(mesh{r});
%   end
% end

%% find itnerface
for r=1:2
  for i=1:size(mesh{r}.bdflag{3},1)
    el = mesh{r}.bdflag{3}(i,1);
    sub = mesh{r}.elem(el,:);
    sub(mesh{r}.bdflag{3}(i,2)) = [];
    bd{r}(i,1:2) = get_rc(mesh{r}, sub);
  end
end

[~,i1,i2] = intersect(bd{1}, bd{2}, 'rows');
interface{1} = [mesh{1}.bdflag{3}(i1,1), mesh{2}.bdflag{3}(i2,1)];
interface{2} = zeros(0,2);

[ustokes,ufemspace,pstokes,pfemspace,phidarcy,phifemspace,mesh,vp] = darcystokesv3(mesh, interface, vp);

%% compute the velocity field in the porous part
qdarcy = reconstruct_velocity(mesh{2}, vp{2}, phidarcy, phifemspace);
udarcy = qdarcy./vp{2}.n;

%% compute the pressure in the porous part
pdarcy = phidarcy.*vp{2}.rho.*vp{2}.g;


%% Visualisation
figure;
simpplot_sol(mesh{1},pstokes(:,1))
simpplot_sol(mesh{2},pdarcy(:,1))
figure;

bary = get_rc(mesh{2});
NS = size(mesh{1}.node,1);
N1 = [mesh{1}.node(:,1); bary(:,1)];
N2 = [mesh{1}.node(:,2); bary(:,2)];
Q1 = [ustokes(1:NS, 1); udarcy(:, 1)];
Q2 = [ustokes(1:NS, 2); udarcy(:, 2)];
quiver(N1,N2,Q1,Q2);
