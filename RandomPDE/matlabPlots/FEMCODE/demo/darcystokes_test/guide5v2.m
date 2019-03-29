%% GUIDE 5 - Darcy-Stokes coupling
% Ondrej Budac, Antoine Imboden, 2014
clear mesh interface vp

mesh = cell(2,1);

%% Parameters
%Volumetric coefficient q = n*Vdarcy, where q is the discharge vector
n = 2;
%Experimental constant, correlated with the average pore size
alpha = 0.08;

%% Stokes part

mesh{1} = structured_mesh([0,2,0,1],[4,2]);
mesh{1}.elem([3,5,6],:) = [];
mesh{1} = renumber(mesh{1});
mesh{1}.bdflag{1} = get_bdflag(mesh{1}, @(x)(x(:,2) == 1)); % dirichlet
mesh{1}.bdflag{2} = get_bdflag(mesh{1}, @(x)(x(:,1).*(2-x(:,1)) == 0)); % Neumann
mesh{1}.bdflag{3} = get_bdflag(mesh{1}, @(x)(x(:,1).*(2-x(:,1)).*(1-x(:,2)) ~= 0)); % interface

%figure;
%simpplot(mesh{1}, struct('nodenum',true,'elemnum',true))

%% Darcy part

mesh{2} = structured_mesh([0,2,-1,0.5],[4,3]);
mesh{2}.elem([17,18,20,23,24],:) = [];
mesh{2} = renumber(mesh{2});
mesh{2}.bdflag{1} = get_bdflag(mesh{2}, @(x)(x(:,2) == -1)); % dirichlet
mesh{2}.bdflag{2} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)) == 0)); % Neumann
%mesh{2}.bdflag{1} = zeros(0,2); % dirichlet
%mesh{2}.bdflag{2} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)).*(x(:,2)+1) == 0)); % Neumann
mesh{2}.bdflag{3} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)).*(1+x(:,2)) ~= 0)); % interface

%figure;
%simpplot(mesh{2}, struct('nodenum',true,'elemnum',true))

vp{1}.a = 1;
vp{1}.f = [0,0];
vp{1}.elemtype = 'p2';
vp{1}.pelemtype = 'p1';
vp{1}.bc = {[1,0], [0,0]};

vp{2}.a = 0.0001 * eye(2);
vp{2}.f = @(x)(repmat([0, 0, 0],size(x,1),1));
vp{2}.elemtype = 'p1';
vp{2}.bc = 0;

%interface{1} = [1, 10; 3, 17; 10, 19; 5, 18; 4, 16]; % stokes - darcy
%interface{2} = zeros(0,2); % darcy - stokes

% these values were read from images. We'll figure out how to do that
% efficiently and automatically
%interface{1} = [1, 10; 3, 17; 10, 19]; % stokes - darcy
%interface{2} = [18, 5; 16,4]; % darcy - stokes

%% we don't allow a triangle to touch the boundary with two sides
for r=1:2
  mesh{r} = bisect(mesh{r},'all');
end

for i=1:3 % increase to get a finer mesh
  for r=1:2
    mesh{r} = uniformrefine(mesh{r});
  end
end
%% find itnerface
for r=1:2
  for i=1:size(mesh{r}.bdflag{3},1)
    el = mesh{r}.bdflag{3}(i,1);
    sub = mesh{r}.elem(el,:);
    sub(mesh{r}.bdflag{3}(i,2)) = [];
    bd{r}(i,1:2) = get_rc(mesh{r}, sub);
  end
end
[~,i1,i2] = intersect(bd{1},bd{2},'rows');
interface{1} = [mesh{1}.bdflag{3}(i1,1), mesh{2}.bdflag{3}(i2,1)];
interface{2} = zeros(0,2);

[usol,ufemspace,psol,pfemspace,phisol,phifemspace,mesh,vp] = darcystokesv2(mesh, interface, vp,n,alpha);

%% compute the velocity field in the porous part
q = reconstruct_velocity(mesh{2}, vp{2}, phisol, phifemspace);


%% Visualisation
figure;
simpplot_sol(mesh{1},psol(:,1))
simpplot_sol(mesh{2},phisol(:,1))
figure;
NS = size(mesh{1}.node,1);
quiver(mesh{1}.node(:,1), mesh{1}.node(:,2), usol(1:NS, 1), usol(1:NS, 2));
hold on;
bary = get_rc(mesh{2});
quiver(bary(:,1), bary(:,2), q(:, 1), q(:, 2));


% figure;
% simpplot_sol(mesh{1},usol(:,1))
% figure;
% simpplot_sol(mesh{1},usol(:,2))