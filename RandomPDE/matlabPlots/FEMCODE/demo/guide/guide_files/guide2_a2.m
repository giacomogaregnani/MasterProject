function tensor = guide2_a2(x, k ,l)
if k~=l
    tensor = zeros(size(x,1),1);
else
    tensor = 1 + 0.45*(sin(x(:,1)*2*pi) + sin(x(:,2)*2*pi));
end
end