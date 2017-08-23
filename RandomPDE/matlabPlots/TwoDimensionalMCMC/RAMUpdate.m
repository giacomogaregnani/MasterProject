function S = RAMUpdate(w, alpha, count, targetAlpha, S)

WWT = w * w';
WTW = w' * w;
nParam = length(w);
gammaI = min(1, 2 * count^(-0.6));
diffAlpha = alpha - targetAlpha;
coeff = gammaI * diffAlpha / WTW;
C = eye(nParam) + WWT * coeff;
S = chol(S * C * S')';

return