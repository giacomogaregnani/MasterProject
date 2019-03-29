clc
clear
idx = 3;
type = 'multi_par';
for i = 1:3
        flux = compute_flux_div_bc_enhanced(i);
        filename = ['flux', num2str(i), type '.mat'];
        save (filename);
end