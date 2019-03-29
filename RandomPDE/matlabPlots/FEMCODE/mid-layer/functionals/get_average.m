function avg = get_average(mesh, femspace, f, varargin)
%GET_AVERAGE Computes an average of a discrete or continuous function

if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end
int = get_integral(mesh, femspace, f, varargin{:});
avg = int / sum(mesh.volume);
end

