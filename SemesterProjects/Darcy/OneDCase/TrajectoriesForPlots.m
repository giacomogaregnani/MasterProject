function TrajectoriesForPlots(X0, f, g, Bounds, W, Time)
% Plot trajectories obtained with DEM and PURE KILLING BC's

[M, N] = size(W);
h = (Time(2)-Time(1))/N;
tau = Time(2) * ones(M,1);

X = zeros(M,N);
X(:,1) = X0;

for j = 1:M
    w = W(j,:);
    for i = 2:N
        X(j,i) = EMOneStep(X(j,i-1),f,g,w(i)-w(i-1),h);
        if X(j,i) > Bounds(2) || X(j,i) < Bounds(1)
            X(j,i+1:end) = X(j,i) * ones(1,N-i);
            break
        end
    end
end

times = linspace(Time(1),Time(2),N);

plot(times,Bounds(1)*ones(size(times)),'k','LineWidth',2)
hold on
plot(times,Bounds(2)*ones(size(times)),'k','LineWidth',2)
plot(times,X')
axis([Time(1) Time(2) Bounds(1)-0.1 Bounds(2)+0.1])


end