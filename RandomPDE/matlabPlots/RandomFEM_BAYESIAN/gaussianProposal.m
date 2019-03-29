function p = gaussianProposal(p, sigma)

d = length(p);
p = p + sigma.* randn(d, 1);

end