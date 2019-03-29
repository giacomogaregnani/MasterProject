function X = EMOneStepDarcy(X,u,sigma,dW,h)
% EM in DARCY case

X = X + u*h + sigma*dW;

end