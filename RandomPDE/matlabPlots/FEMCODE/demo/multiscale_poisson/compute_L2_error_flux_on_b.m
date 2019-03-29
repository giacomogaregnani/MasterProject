function err = compute_L2_error_flux_on_b(bound, f1, f2)

N1 = size(f1,1); % can do once
N2 = size(f2,1); % can do once
x1 = linspace(bound(1),bound(2),N1); % can do once
x2 = linspace(bound(1),bound(2),N2); % can do once
f21 = interp1(x2, f2(:,1), x1); f22 = interp1(x2, f2(:,2), x1);
f23 = interp1(x2, f2(:,3), x1); f24 = interp1(x2, f2(:,4), x1);
f2i = [f21', f22', f23', f24'];

err1 = ((trapz(x1, ((f1(:,1)-f2i(:,1))).^2)));
err2 = ((trapz(x1, (f1(:,2)-f2i(:,2)).^2)));
err3 = ((trapz(x1, (f1(:,3)-f2i(:,3)).^2)));
err4 = ((trapz(x1, ((f1(:,4)-f2i(:,4))).^2)));
err = sqrt(err1+err2+err3+err4);