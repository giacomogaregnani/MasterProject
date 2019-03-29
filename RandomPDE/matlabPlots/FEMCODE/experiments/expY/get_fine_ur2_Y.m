function [mesh, father] = get_fine_ur2_Y(RB)
%creates a fine mesh for experiment A using the uniform refinement 6 times
mesh = RB.rmesh;
nref = 5;
[mesh, father] = uniformrefine(mesh, RB.R);
[mesh, father] = uniformrefine(mesh, father);
end
