function y = rkcProb(f, N, y, h, s)

kOld = y;
kCurr = y;
coeff = h / s^2;
coeff2 = 2 * coeff;

for i = 1 : N
    H = h + h^1.5 * (2 * rand(1) - 1);
    kNew = y + f(H, y) * coeff;
    for j = 1 : s
        kCurr = f(H, kNew) * coeff2 + 2 * kNew - kOld;
        kOld = kNew;
        kNew = kCurr;
    end
    y = kCurr;
end

end