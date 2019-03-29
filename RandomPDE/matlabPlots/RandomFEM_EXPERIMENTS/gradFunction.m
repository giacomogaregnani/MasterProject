function grad = gradFunction(x, uGrad, xq)

if size(xq) == [1, 1]
    grad = uGrad(find(x > xq, 1)-1);
else
    Nq = length(xq);
    grad = zeros(1, Nq);
    for i = 1 : Nq
        if xq(i) <= x(1) || xq(i) >= x(end)
            grad(i) = 0;
        else
            grad(i) = uGrad(find(x > xq(i), 1)-1);
        end
    end
end
end