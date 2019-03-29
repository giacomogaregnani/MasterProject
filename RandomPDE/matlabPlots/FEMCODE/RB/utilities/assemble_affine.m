function A = assemble_affine(mu, Theta, X)
Q = numel(Theta);
if issparse(X{1})
    A = sparse(size(X{1},1), size(X{1},2));
else
    A = 0;
end
for q=1:Q
    A = A + evalmu(Theta{q}, mu) * X{q};
end
end