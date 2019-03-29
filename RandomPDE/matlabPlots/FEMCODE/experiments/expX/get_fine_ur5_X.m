function [mesh, father] = get_fine_ur5_X(RB)
%creates a fine mesh for experiment A using the uniform refinement 6 times
mesh = RB.rmesh;
nref = 5;
[mesh, father] = bisect(mesh, [2,3], RB.R);
[mesh, father] = uniformrefine(mesh, father);
for iref = 2:nref
    [mesh, father] = uniformrefine(mesh, father);
end
mesh = periodize(mesh,[-1/2,1/2,-1/2,1/2]);
end
