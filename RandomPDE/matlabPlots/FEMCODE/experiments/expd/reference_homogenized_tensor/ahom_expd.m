function itensor = ahom_expd(x,k,l)
%A_EXP_B Summary of this function goes here
%   Detailed explanation goes here

load('ahom_expd.mat','R','cs');

%% RAW DATA - NEED AVERAGING
r = x2par_expd(x);
itensor = ppval(cs{k,l},r);
end

