function [p] = NaiveForPlots(X0,f,g,W,Time,p)
% Function used for PLOTTING trajectories


if X0(1) >= 1 || X0(1) <= -1 || X0(2) >= 1 || X0(2) <= -1
    return
end

[TwoM,N] = size(W);
M = TwoM/2;
h = (Time(2)-Time(1))/(N-1);
tau = Time(2) * ones(M,1);

for j = 1:M
    w = W(2*j-1:2*j,:);
    x(:,1) = X0;
    for i = 2:N
        x(:,i) = EMOneStep(x(:,i-1),f,g,w(:,i)-w(:,i-1),h);
%         if x(1,i) >= 1 || x(1,i) <= -1 || x(2,i) >= 1 || x(2,i) <= -1
%             plot(x(1,:),x(2,:));
%             break
%         end
    end
end
plot(x(1, :), x(2, :), 'LineWidth', 1.5)

end