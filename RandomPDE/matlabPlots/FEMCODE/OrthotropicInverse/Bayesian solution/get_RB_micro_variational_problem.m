function vp = get_RB_micro_variational_problem(pde)

vp.elemtype = 'p1';
vp.a = pde.micTensorQ;
vp.a_full = pde.micTensor;
vp.ThetaField = pde.ThetaField;
vp.f = permute([zeros(3,2), -eye(3)], [3, 2, 1]);
vp.type = 'micro_linear_elasticity_orth';
end