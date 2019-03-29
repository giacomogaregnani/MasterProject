function mesh = cleanfields(mesh)
%CLEANFIELDS cleans auxiliary fields in a mesh: volume, dlambda, bary, ...

fieldlist = {'volume', 'dlambda','bary'};
for i=1:numel(fieldlist)
	if isfield(mesh,fieldlist{i})
		mesh=rmfield(mesh,fieldlist{i});
	end
end
end

