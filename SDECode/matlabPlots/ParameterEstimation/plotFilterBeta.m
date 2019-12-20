clc; clear; close all

eps = 0.1;
alpha = 0.5;
d = eps^alpha;
tVec = linspace(0, 1, 1000);

figure
hold on
for b = 1 : 5 : 40
    Cb = 1/gamma((b+1)/b);
    k = @(t, s) Cb / d * exp(-abs(t-s).^b / d^b);
    plot(tVec, k(1, tVec))
end

%% 

eps = 0.01;
b = 2;
Cb = 1/gamma((b+1)/b);
tVec = linspace(0, 1, 10000);

% figure
% hold on
% for alpha = 1
%     d = eps^alpha;
%     k = @(t, s) Cb / d * exp(-abs(t-s).^b / d);
%     kTaylor = @(t, s) max(0, Cb / d - Cb / (d^(b+1)) * (t-s).^b);
%     plot(tVec, k(tVec(end), tVec))
%     plot(tVec, kTaylor(tVec(end), tVec))
%     trapz(tVec, kTaylor(tVec(end), tVec))
% end
 
figure
hold on
for alpha = 0:0.5:1
    d = eps^alpha;
    k = @(t, s) Cb / d^(1/b) * exp(-abs(t-s).^b / d);
    plot(tVec, k(tVec(end), tVec))
end
