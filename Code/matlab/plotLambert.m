d = 0;
gamma = d - 0.5;
M = @(acc) exp(1/gamma * lambertw(-1, gamma * acc));
h = @(q,acc) acc.^(1/q);
cost = @(q,acc) M(acc) .* (h(1,acc).^(-1));

accTest = [0.35 : -0.01 : 0.04]

N = [10 100 1000 10000 10000];
figure
loglog(N, 1./sqrt(N));
hold on
loglog(N, log(N) ./ sqrt(N));
loglog(M(accTest), accTest)