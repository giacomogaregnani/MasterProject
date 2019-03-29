function y = rkc(f, N, y, h, s)

kOld = y;
kCurr = y;
coeff = h / s^2;
coeff2 = 2 * coeff;

for i = 1 : N
    kNew = y + f(h, y) * coeff;
    for j = 1 : s
        kCurr = f(h, kNew) * coeff2 + 2 * kNew - kOld;
        kOld = kNew;
        kNew = kCurr;
    end
    y = kCurr;
end

end