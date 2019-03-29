function w = GaussianProposal(w, sigma)

w = w + sigma * randn(size(w));