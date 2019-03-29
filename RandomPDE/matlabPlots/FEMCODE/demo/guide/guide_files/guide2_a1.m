function tensor = guide2_a1(x, k ,l)
if k~=l
    tensor = zeros(size(x,1),1);
else
    tensor = 0.1 + x(:,k).^2;
end
end