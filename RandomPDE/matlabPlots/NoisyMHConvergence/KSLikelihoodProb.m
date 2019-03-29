function l = KSLikelihoodProb(w, xObs, uObs, obsNoise, A, F, x, m, sigma)

uVec = [0; 1 / exp(w) * (A \ F); 0];
u = interp1(x, uVec, xObs);

diff = uObs - u - m;
% tmp = uObs - u + m;

% if abs(tmp) < abs(diff)
%     diff = tmp;
% end

l = -0.5 / (obsNoise^2 + sigma^2) * (diff' * diff);
