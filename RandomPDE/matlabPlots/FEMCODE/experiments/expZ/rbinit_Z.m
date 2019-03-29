function RB = rbinit_Z(RB)
% initialization of the RB experiment
RB.param.ref = [0,0,0]; % sets the micro mesh to the rectangular L-shape
RB.param.min = [-1/18,-1/18,-1/18];    
RB.param.max = [ 1/18, 1/18, 1/18];
RB.param.test = @(x)(true(size(x,1),1));

mesh = [];
load('micro_fmesh_Z.mat','mesh');
RB.rmesh = mesh;

s1ref = [-1/2, -1/6, 1/6, 1/2];
s2ref = [-1/2, -1/6, 1/6, 1/2];
s3ref = [-1/2, -1/6, 1/6, 1/2];

mu1 = sym('mu1');
mu2 = sym('mu2');
mu3 = sym('mu3');
s1    = [-1/2, -1/6-mu1, 1/6+mu1, 1/2];
s2    = [-1/2, -1/6-mu2, 1/6+mu2, 1/2];
s3    = [-1/2, -1/6-mu3, 1/6+mu3, 1/2];

NT = size(mesh.elem,1);
bary = get_rc(mesh);
[RB.G,RB.C] = deal(cell(27,1));
RB.R = zeros(NT,1);
for i=1:3
  for j=1:3
    for k=1:3
    ind = 1+(i-1)+3*(j-1)+9*(k-1);
    wh = (bary(:,1) >= s1ref(i)) & (bary(:,1) <= s1ref(i+1)) & ...
      (bary(:,2) >= s2ref(j)) & (bary(:,2) <= s2ref(j+1)) & ...
      (bary(:,3) >= s3ref(k)) & (bary(:,3) <= s3ref(k+1));

    RB.R(wh) = ind;
    RB.G{ind} = sym(zeros(3,3));
    RB.G{ind}(1,1) = (s1(i+1)-s1(i))/(s1ref(i+1)-s1ref(i));
    RB.G{ind}(2,2) = (s2(j+1)-s2(j))/(s2ref(j+1)-s2ref(j));
    RB.G{ind}(3,3) = (s3(k+1)-s3(k))/(s3ref(k+1)-s3ref(k));
    
    RB.C{ind} = sym(zeros(3,1));
    RB.C{ind}(1) = s1(i) - s1ref(i)*RB.G{ind}(1,1);
    RB.C{ind}(2) = s2(j) - s2ref(j)*RB.G{ind}(2,2);
    RB.C{ind}(3) = s3(k) - s3ref(k)*RB.G{ind}(3,3);
    end
  end
end
end

