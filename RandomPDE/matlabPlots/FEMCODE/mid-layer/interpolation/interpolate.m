function res = interpolate(mesh, femspace, f)
%INTERPOLATE nodal interpolation of a function handle f.
xloc = get_loc(mesh, femspace);
res = f(xloc);
end

