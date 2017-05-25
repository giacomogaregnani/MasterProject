function I = coordToInd(origin, r, p)

I = ones(1, length(r)) + round((p - origin) ./ r);
