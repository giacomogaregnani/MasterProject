%% GUIDE 5 - Darcy-Stokes coupling
% Ondrej Budac, Antoine Imboden, 2014
clear mesh interface vp

mesh = cell(2,1);


%% Stokes part

mesh{1} = structured_mesh([0,2,0,1],[4,2]);
mesh{1}.elem([3,5,6],:) = [];
mesh{1} = renumber(mesh{1});
mesh{1}.bdflag{1} = get_bdflag(mesh{1}, @(x)(x(:,2) == 1)); % dirichlet
mesh{1}.bdflag{2} = get_bdflag(mesh{1}, @(x)(x(:,1).*(2-x(:,1)) == 0)); % Neumann
mesh{1}.bdflag{3} = get_bdflag(mesh{1}, @(x)(x(:,1).*(2-x(:,1)).*(1-x(:,2)) ~= 0)); % interface

%figure;
%simpplot(mesh{1}, struct('nodenum',true,'elemnum',true))


mesh{2} = structured_mesh([0,2,-1,0.5],[4,3]);
mesh{2}.elem([17,18,20,23,24],:) = [];
mesh{2} = renumber(mesh{2});
mesh{2}.bdflag{1} = get_bdflag(mesh{2}, @(x)(x(:,2) == -1)); % dirichlet
mesh{2}.bdflag{2} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)) == 0)); % Neumann
mesh{2}.bdflag{3} = get_bdflag(mesh{2}, @(x)(x(:,1).*(2-x(:,1)).*(1+x(:,2)) ~= 0)); % interface

%figure;
%simpplot(mesh{2}, struct('nodenum',true,'elemnum',true))

vp{1}.a = 1;
vp{1}.f = [1,0];
vp{1}.elemtype = 'p2';
vp{1}.pelemtype = 'p1';
vp{1}.bc = [0,0];

vp{2}.a = 0.1 * eye(2);
vp{2}.f = 0;
vp{2}.elemtype = 'p1';
vp{2}.bc = 0;

%interface{1} = [1, 10; 3, 17; 10, 19; 5, 18; 4, 16]; % stokes - darcy
%interface{2} = zeros(0,2); % darcy - stokes

% these values were read from images. We'll figure out how to do that
% efficiently and automatically
interface{1} = [1, 10; 3, 17; 10, 19]; % stokes - darcy
interface{2} = [18, 5; 16,4]; % darcy - stokes

darcystokes(mesh, interface, vp);