phiFP = @(N, error) log(N)./error;

eTest = 10.^[-1:-0.01:-10];
N0 = 1000 * ones(size(eTest));
N = N0;

for i = 1 : length(N0)
    for k = 1 : 1000
        N(i) = phiFP(N(i), eTest(i));  
    end
end

loglog(eTest, N)