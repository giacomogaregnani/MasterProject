function y = impEulerLinear(A, y0, tSpan, h)

t = tSpan(1) : h : tSpan(2);
y = y0;

itMatrix = eye(length(y0)) - h * A;
for i = 2 : length(t)
    y = itMatrix \ y;
end