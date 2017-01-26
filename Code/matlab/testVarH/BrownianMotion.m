function W = BrownianMotion(Time, N, M, size)
% Compute BROWNIAN MOTION in 2D case

h = (Time(2)-Time(1)) / N;
W = zeros(size*M,N+1);
W(:,2:N+1) = cumsum(sqrt(h) * randn(size*M,N), 2);

end