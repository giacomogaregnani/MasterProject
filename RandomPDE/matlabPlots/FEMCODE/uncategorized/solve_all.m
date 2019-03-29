function sol = solve_all(A,rhs)
%SOLVE_ALL is a vectorized way to solve a large number of small linear eqs
%
% sol = solve_all(A, rhs) gives solution to many same-size linear systems
% in a vectorized manner. For the solution we use Crammer's rule which can
% give bad results in some cases. Be aware of this imprecision.
%
% for any i we have system A(i,:,:) * sol(i,:) = rhs(i,:).

%% INITIALIZATION
NT = size(A,1);
dim = size(A,2);
mdet = detn(A);
tol = 1e-30;    % tolerance for singularity (good idea ?!)
nonsing = abs(mdet)>tol;

%% Crammer's rule
sol = zeros(NT,dim);
B = cat(3, A, reshape(rhs,[NT, dim, 1])); 
for i=1:dim
	sol(:,i) = detn(B, 1:dim, [1:i-1, dim+1, i+1:dim]);
end
sol(nonsing,:) = sol(nonsing,:)./repmat(mdet(nonsing), [1,dim]);
sol(~nonsing,:) = Inf;

end
