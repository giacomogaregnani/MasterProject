function diam = simplex_diameter(mesh)
%SIMPLEX_DIAMETER returns the diameter of all elements
% diam = simplex_diameter(mesh) computes the longes edge for every element
% in a mesh.

dim = size(mesh.elem,2) - 1;
diam = max(subsimplex_volume(mesh, dim - 1), [], 2);
end