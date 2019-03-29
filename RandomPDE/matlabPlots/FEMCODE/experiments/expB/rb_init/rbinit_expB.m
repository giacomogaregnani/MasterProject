function RB = rbinit_expB(RB)
% initialization of the RB experiment
RB.param.ref = [0,0]; % sets the micro mesh to the rectangular L-shape
RB.param.min = [-0.2,-0.2];    
RB.param.max = [ 0.2, 0.2];
RB.param.test = @(x)(true(size(x,1),1));

mu1 = sym('mu1');
mu2 = sym('mu2');
RB.mesh.elem = [2,6,1; 5,6,2; 4,5,2; 4,2,3];
RB.mesh.node = [-1/2,0; ...
  mu1, mu2; ...
  0,-1/2; 1/2,-1/2; 1/2,1/2; -1/2,1/2];
end

