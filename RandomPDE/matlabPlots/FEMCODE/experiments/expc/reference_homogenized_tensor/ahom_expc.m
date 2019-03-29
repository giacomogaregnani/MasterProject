function itensor = ahom_expc(x,k,l)
%A_EXP_B Summary of this function goes here
%   Detailed explanation goes here

load('ahom_expc.mat','R','cs');

%% RAW DATA - NEED AVERAGING
r = x2par_expc(x);
itensor = ppval(cs{k,l},r);
end

