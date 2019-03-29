function itensor = ahom_expe(x,k,l)
%A_EXP_B Summary of this function goes here
%   Detailed explanation goes here

load('ahom_expe.mat','cs');

%% RAW DATA - NEED AVERAGING
r = x2par_expe(x);
itensor = ppval(cs{k,l},r);
end

