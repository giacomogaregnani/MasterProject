function W = BrownianMotion2D(Time,N,M)
% BM in 2D

h = (Time(2)-Time(1)) / N;
W = zeros(2*M,N+1);
W(:,2:N+1) = cumsum(sqrt(h) * randn(2*M,N),2);

end