clc; clear; close all
%%

b = 3;
d = 0.1;
Cb = 1/gamma((b+1)/b);
tVec = linspace(0, 1, 1000);
myigamma = @(a, x) gamma(a) - igamma(a, x);
k = @(t, s) Cb / (d^(1/b)) * exp(-(t-s).^b / d);

%% PHI

phi = @(t) 1 - myigamma(1/b, t.^b/d) / gamma(1/b);

checkPhi = zeros(size(tVec));
for i = 1 : length(tVec)
    checkPhi(i) = 1 - integral(@(s) k(tVec(i),s), 0, tVec(i));
end

plot(tVec, phi(tVec));
hold on
plot(tVec, checkPhi)

display(['error phi = ' , num2str(max(abs(phi(tVec)-checkPhi)))])

%% INTEGRAL CHECK

I1 = @(t) b/gamma(1/b)^2 * 1/(2*d)^(1/b) * myigamma(1/b, 2*t.^b/d);
I2 = @(t) b/gamma(1/b)^2 * exp(-t.^b/d) / d^(1/b) .* myigamma(1/b, t.^b/d);
I3 = @(t) b^2/gamma(1/b)^2 * t/d^(2/b) .* exp(-2*t.^b/d);

integrand = @(t, s) (k(t,s)-k(t,0)).^2;

iNumerics = zeros(size(tVec));
for i = 1 : length(tVec)
    iNumerics(i) = integral(@(s) integrand(tVec(i),s), 0, tVec(i));
end
iTheory = I1(tVec) - 2*I2(tVec) + I3(tVec);

figure
plot(tVec, iNumerics)
hold on
plot(tVec, iTheory)

display(['error integral = ' , num2str(max(abs(iTheory-iNumerics)))])
