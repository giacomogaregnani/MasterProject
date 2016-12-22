function H = computeNormalHellDistance(mu1, mu2, sigma1, sigma2)

S = (sigma1 + sigma2) / 2;
detS = det(S);
m = (mu1 - mu2)';

H = (det(sigma1) * det(sigma2))^0.25 / sqrt(detS) * ...
    exp(-0.0125 * m' * (S \ m));

H = sqrt(1 - H);
 
% Explicit computation
% dS1 = det(2 * pi * sigma1);
% dS2 = det(2 * pi * sigma2);
% C1 = 1 / sqrt(dS1);
% C2 = 1 / sqrt(dS2);
% 
% invS1 = sigma1 \ eye(3);
% invS2 = sigma2 \ eye(3);
% 
% n = 400;
% x1 = linspace(range(1, 1), range(2, 1), n);
% x2 = linspace(range(1, 2), range(2, 2), n); 
% x3 = linspace(range(1, 3), range(2, 3), n);
% 
% [XX1, XX2, XX3] = meshgrid(x1, x2, x3);
% XYZ = [XX1(:), XX2(:), XX3(:)];
% clear XX1 XX2 XX3
% 
% evalThreeGauss = @(X, S, C) C * exp(-0.5 * (S(1, 1) * X(:, 1) .^2 + ...
%                                             S(2, 2) * X(:, 2) .^2 + ...
%                                             S(3, 3) * X(:, 3) .^2 + ...
%                                         2 * S(1, 2) * X(:, 1) .* X(:, 2) + ...
%                                         2 * S(1, 3) * X(:, 1) .* X(:, 3) + ...
%                                         2 * S(2, 3) * X(:, 2) .* X(:, 3)));
% 
% Y1 = evalThreeGauss(XYZ - repmat(mu1, n^3, 1), invS1, C1);
% Y2 = evalThreeGauss(XYZ - repmat(mu2, n^3, 1), invS2, C2);
% 
% clear XYZ
% 
% Y1 = reshape(Y1, n, n, n);
% Y2 = reshape(Y2, n, n, n);
% 
% H = sqrt(0.5 * trapz(x1, trapz(x2, trapz(x3, (sqrt(Y1) - sqrt(Y2)).^2, 3), 2), 1));

end

