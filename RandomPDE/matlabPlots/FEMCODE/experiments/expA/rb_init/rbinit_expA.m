function RB = rbinit_expA(RB)
% initialization of the RB experiment
RB.param.ref = [0,3/4]; % sets the micro mesh to the rectangular L-shape
RB.param.min = [0,0];    
RB.param.max = [2,3];
RB.param.test = @(x)(~((x(:,2)>1) & (x(:,2) < 2) & (x(:,1) > 1)));

mu1 = sym('mu1');
mu2 = sym('mu2');
RB.mesh.elem = [2,6,1; 5,6,2; 4,5,2; 4,2,3];
RB.mesh.node = [-1/2,0; ...
  cos(2*pi*(mu2-mu1)/3)/5, cos(2*pi*(mu2+mu1)/3)/5; ...
  0,-1/2; 1/2,-1/2; 1/2,1/2; -1/2,1/2];
end

