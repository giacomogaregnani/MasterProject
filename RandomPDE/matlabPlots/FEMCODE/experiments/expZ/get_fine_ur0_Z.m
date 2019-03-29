function [mesh, father] = get_fine_ur0_Z(RB)
%creates a fine mesh for experiment A using the uniform refinement 6 times
mesh = RB.rmesh;
father = RB.R;
end
