function ThetaField = evalTheta(mu, Theta)
Q = numel(Theta);
NP = size(mu,1);
ThetaField = zeros(NP,Q);
for q=1:Q
    ThetaField(:,q) = evalmu(Theta{q}, mu);
end
end