function [proj, inface] = project(points, subspace)
%% Project points p on the subplane spanned by node
dim = size(subspace,1)-1;
NP= size(points,1);
A = zeros(dim,dim);
V = subspace(2:dim+1,:)-repmat(subspace(1,:),[dim,1]); 
for i=1:dim
        for j=1:i
                A(i,j) = dot(V(i,:), V(j,:));
                A(j,i) = A(i,j);
        end
end
proj = (points - repmat(subspace(1,:), [NP,1])) * V' * A^(-1);
if nargout >1
	inface = ((sum(proj >=0 ,2) == dim) .* (sum(proj,2) <= 1)) > 0;
end
proj = repmat(subspace(1,:), [NP,1]) + proj * V;
end