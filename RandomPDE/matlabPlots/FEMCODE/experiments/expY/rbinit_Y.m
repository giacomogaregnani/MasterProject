function RB = rbinit_Y(RB)
% initialization of the RB experiment
RB.param.ref = [0,0]; % sets the micro mesh to the rectangular L-shape
RB.param.min = [-1/12,-1/12];    
RB.param.max = [1/12,1/12];
RB.param.test = @(x)(true(size(x,1),1));

mesh = [];
load('micro_rmesh_Y.mat','mesh');
RB.rmesh = mesh;

s1ref = [-1/2, -1/6, 1/6, 1/2];
s2ref = [-1/2, -1/6, 1/6, 1/2];

mu1 = sym('mu1');
mu2 = sym('mu2');
s1    = [-1/2, -1/6-mu1, 1/6+mu1, 1/2];
s2    = [-1/2, -1/6-mu2, 1/6+mu2, 1/2];

NT = size(mesh.elem,1);
bary = get_rc(mesh);
[RB.G,RB.C] = deal(cell(9,1));
RB.R = zeros(NT,1);
for i=1:3
  for j=1:3
    ind = i+3*(j-1);
    wh = (bary(:,1) >= s1ref(i)) & (bary(:,1) <= s1ref(i+1)) & (bary(:,2) >= s2ref(j)) & (bary(:,2) <= s2ref(j+1));
    RB.R(wh) = ind;
    RB.G{ind} = sym(zeros(2,2));
    RB.G{ind}(1,1) = (s1(i+1)-s1(i))/(s1ref(i+1)-s1ref(i));
    RB.G{ind}(2,2) = (s2(j+1)-s2(j))/(s2ref(j+1)-s2ref(j));
    
    RB.C{ind} = sym(zeros(2,1));
    RB.C{ind}(1) = s1(i) - s1ref(i)*RB.G{ind}(1,1);
    RB.C{ind}(2) = s2(j) - s2ref(j)*RB.G{ind}(2,2);
  end
end
end

