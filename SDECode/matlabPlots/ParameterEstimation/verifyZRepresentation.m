% Numerical Verification of Lemma 3.6 (Representation of dZ)
clc; clear; close all
%% Create the process

T = 100;
N = 10000;
h = T / N;
timeVec = linspace(0, T, N+1);

dW = sqrt(h)*randn(1, N);

gradV = @(x) x;
sigma = 0.5;
X = zeros(N+1, 1); 
for i = 1 : N
    X(i+1) = X(i) - h * gradV(X(i)) + sqrt(2*sigma) * dW(i);
end

%% Representation formula

delta = 100*h;
k = @(t, s) 1/delta * exp(-(t-s)/delta);
dtk = @(t, s) -1/delta^2 * exp(-(t-s)/delta);

% Build the integral
Y = zeros(N+1, 1);
for i = 2 : N+1 
    dtkEval = dtk(timeVec(i), timeVec(1:i-1))';
    xDiff = X(1:i-1) - X(i);
    Y(i) = h * sum(dtkEval .* xDiff);
end

% Build the process
Z = zeros(N+1, 1);
Z(1) = X(1);
for i = 1 : N
    Z(i+1) = Z(i) + h * (k(timeVec(i), 0) * X(i) + Y(i));
end

%% Convolution

kEval = k(timeVec, 0); 
Z2 = conv(kEval, X) * h;
Z2 = Z2(1:N+1);


%% Plot

figure
plot(timeVec, X, 'color', 0.8*ones(1, 3))
hold on
plot(timeVec, Z, 'k')
plot(timeVec, Z2, 'k--')