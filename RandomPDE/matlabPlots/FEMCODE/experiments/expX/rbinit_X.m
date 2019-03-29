function RB = rbinit_X(RB)
% initialization of the RB experiment
RB.param.ref = [0,0]; % sets the micro mesh to the rectangular L-shape
RB.param.min = [-0.2,-0.2];    
RB.param.max = [0.2,0.2];
RB.param.test = @(x)(true(size(x,1),1));

mu1 = sym('mu1');
mu2 = sym('mu2');

mesh = struct;
mesh.node = [-1/2,0; mu1, mu2; 0,-1/2; 1/2,-1/2; 1/2,1/2; -1/2,1/2];
mesh.elem =  [2,6,1; 5,6,2; 4,5,2; 4,2,3];

RB.rmesh = struct;
RB.rmesh.node = subsmu(mesh.node, RB.param.ref);
RB.rmesh.elem = mesh.elem;

[RB.C, RB.G, RB.R] = get_affine_transformation(mesh, RB.param.ref);
end

