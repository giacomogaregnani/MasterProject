function father = find_father(cmesh, fmesh)

NT = size(fmesh.elem, 1);
bary = get_rc(fmesh);

father = zeros(NT, 1);
for i=1:size(cmesh.elem, 1)
  lambda2 = get_lc2(cmesh, i * ones(NT,1), bary);
  father(all(lambda2>=0, 2)) = i;
end
end