function sigma = batchMeansEstimator(theta)
%% Every column of theta must be an observations (i.e., theta [nParam, nSamples]) 

[nParam, nSamples] = size(theta);

batchSize = floor((nSamples)^(1/2));
nBatch = floor(nSamples / batchSize);

batchMeans = zeros(nParam, nBatch);

for k = 1 : nBatch
    for i = 1 : batchSize
        batchMeans(:, k) = batchMeans(:, k) + theta(:, (k-1)*batchSize + i);        
    end
    batchMeans(:, k) = batchMeans(:, k) / batchSize;
end

totalMean = mean(theta, 2);
sigma = 0;
for k = 1 : nBatch
    sigma = sigma + (batchMeans(:, k) - totalMean).^2;
end
sigma = sigma * batchSize / (nBatch - 1);
end