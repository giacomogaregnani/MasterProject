function plot_error( err )
%PLOT_ERROR Summary of this function goes here
%   Detailed explanation goes here

params={'interpreter','latex','FontSize',14};
figure;
loglog(...
	err.dof, err.uL2(:,1,1), '-x', ...
	err.dof, err.uH1(:,1,1), '-x', ...
	err.dof, err.pL2(:,1,1), '-x', ...
	err.dof, err.est, '-o', ...
	[10^4,10^6], [10^-1,10^-3], ...
	[10^4,10^6], [10^-3,10^(-3-8/3)]);
legend({...
	'$\| u-u^h\|_{L^2}$', ...
    '$| u-u^h|_{H^1}$', ...
    '$\| p-p^h\|_{L^2}$', ...
	'$\mu_\Omega=\left(\sum_{K\in \mathcal{T}_H} \mu_K^2\right)^{1/2}$' ...
	'$O(N_{dof}^{-3/2})$', ...
	'$O(N_{dof}^{-2})$'
	},params{:},'location','southwest');
xlabel('$N_{dof}$',params{:})
end
