function l = likFunction(x, y, sigmaNoise)

l = exp(-(x - y)^2 / (2 * sigmaNoise^2));

return