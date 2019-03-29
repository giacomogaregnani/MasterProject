function [mesh, father] = get_fine_expA_ur5(RB)
%creates a fine mesh for experiment A using the uniform refinement 6 times
mesh = struct;
mesh.elem = RB.mesh.elem;
mesh.node = subsmu(RB.mesh.node, RB.param.ref);
nref = 5;
[mesh, father] = bisect(mesh, [2,3]);
[mesh, father] = uniformrefine(mesh, father);
mesh = periodize(mesh,[-1/2,1/2,-1/2,1/2]);
for iref = 2:nref
    [mesh, father] = uniformrefine(mesh, father);
end
end