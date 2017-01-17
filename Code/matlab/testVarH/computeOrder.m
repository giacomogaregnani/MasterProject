function orders = computeOrder(h, X)

orders = log(X(1 : end - 1) ./ X(2 : end)) ./ log(h(1 : end - 1) ./ h(2 : end)); 

