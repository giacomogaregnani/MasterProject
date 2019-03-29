function vp = get_micro_variational_problem(pde)

vp.elemtype = 'p1';
vp.a = pde.micTensor;
vp.f = permute([zeros(3,2), -eye(3)], [3, 2, 1]);
vp.type = 'micro_linear_elasticity_orth';
end